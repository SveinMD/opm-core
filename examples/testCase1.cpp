/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "config.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cassert>
#include <opm/core/grid.h>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/io/vtk/writeVtkData.hpp>
#include <opm/core/linalg/LinearSolverUmfpack.hpp>
#include <opm/core/pressure/IncompTpfa.hpp>
#include <opm/core/pressure/FlowBCManager.hpp>
#include <opm/core/props/IncompPropertiesBasic.hpp>

#include <opm/core/transport/reorder/TransportSolverTwophaseReorder.hpp>

#include <opm/core/simulator/TwophaseState.hpp>
#include <opm/core/simulator/WellState.hpp>

#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/core/utility/StopWatch.hpp>

using std::string;

/// The Darcy law gives
///   \f[u_\alpha= -\frac1{\mu_\alpha} K_\alpha\nabla p_\alpha\f]
/// where \f$\mu_\alpha\f$ and \f$K_\alpha\f$ represent the viscosity
/// and the permeability tensor for each phase \f$\alpha\f$. In the two phase
/// case, we have either \f$\alpha=w\f$ or \f$\alpha=o\f$.
/// In this tutorial, we do not take into account capillary pressure so that
/// \f$p=p_w=p_o\f$ and gravity
/// effects. We  denote by \f$K\f$ the absolute permeability tensor and each phase
/// permeability is defined through its relative permeability by the expression
/// \f[K_\alpha=k_{r\alpha}K.\f]
/// The phase mobility are defined as
///  \f[\lambda_\alpha=\frac{k_{r\alpha}}{\mu_\alpha}\f]
/// so that the Darcy law may be rewritten as
///  \f[u_\alpha= -\lambda_\alpha K\nabla p.\f]
/// The conservation of mass for each phase writes:
/// \f[\frac{\partial}{\partial t}(\phi\rho_\alpha s_\alpha)+\nabla\cdot(\rho_\alpha u_\alpha)=q_\alpha\f]
/// where \f$s_\alpha\f$ denotes the saturation of the phase \f$\alpha\f$ and \f$q_\alpha\f$ is a source term. Let
/// us consider a two phase flow with oil and water. We assume that the rock and both fluid phases are incompressible. Since
/// \f$s_w+s_o=1\f$, we may add the conservation equations to get
///  \f[ \nabla\cdot u=\frac{q_w}{\rho_w}+\frac{q_o}{\rho_o}.\f]
/// where we define
///   \f[u=u_w+u_o.\f]
/// Let the total mobility be equal to
/// \f[\lambda=\lambda_w+\lambda_o\f]
/// Then, we have
/// \f[u=-\lambda K\nabla p.\f]
/// The set of equations
/// \f[\nabla\cdot u=\frac{q_w}{\rho_w}+\frac{q_o}{\rho_o},\quad u=-\lambda K\nabla p.\f]
/// is referred to as the <strong>pressure equation</strong>. We introduce
/// the fractional flow \f$f_w\f$
/// as
/// \f[f_w=\frac{\lambda_w}{\lambda_w+\lambda_o}\f]
/// and obtain
/// \f[\phi\frac{\partial s_w}{\partial t}+\nabla\cdot(f_w u)=\frac{q_w}{\rho_w}\f]
/// which is referred to as the <strong>transport equation</strong>. The pressure and
/// transport equation are coupled. In this tutorial, we implement a splitting scheme,
/// where, at each time step, we decouple the two equations. We solve first
/// the pressure equation and then update the water saturation by solving
/// the transport equation assuming that \f$u\f$ is constant in time in the time step
/// interval we are considering.

string replaceStrChar(string str, const string & replace, char ch)
{
	size_t found = str.find_first_of(replace);
	while( found != string::npos )
	{
		str[found] = ch;
		found = str.find_first_of(replace, found+1);
	}
	
	return str;
}

void printIterationsFromVector(const Opm::TransportSolverTwophaseReorder & transport_solver, int i, int num_cells, const char solver_type, const double comp_length, const double time_step)
{
	std::vector<int> iterations = transport_solver.getReorderIterations();
	std::ostringstream iterfilename;
	
	string str_comp_length = replaceStrChar(std::to_string(comp_length), ".", '_');
	string str_time_step = replaceStrChar(std::to_string(time_step), ".", '_');
	
	// Set filename and open
	iterfilename.str(""); 
	iterfilename << "testCase1-iterations-s-" << solver_type << "-T-" << str_comp_length << "-t-" << str_time_step << "-" << std::setw(3) << std::setfill('0') << i << ".txt";
	std::ofstream file; file.open(iterfilename.str().c_str());
	for ( int i = 0; i < num_cells; i++)
	{
		file << i << "\t" << iterations[i] << "\n";
	}
	file.close();
}

int main (int argc, char ** argv)
try
{
	double time_step_days = 0.1;
	double comp_length_days = 2;
	char solver_type = 'r';
	bool printIterations = false;
	// Check parameters
	if(argc > 1)
	{
		// -n: Newton solver
		// -r: Regula Falsi
		for(int i = 1; i < argc; i++)
		{
			if(std::string(argv[i]) == "-s")
			{
				i++;
				if(std::string(argv[i]) == "n")
				{
					std::cout << "Newton solver chosen for single cell problem\n";
					solver_type = 'n';
				}
				else if(std::string(argv[i]) == "t")
				{
					std::cout << "Trust region solver chosen for single cell problem\n";
					solver_type = 't';
				}
				else if(std::string(argv[i]) == "i")
				{
					std::cout << "Ridder's solver chosen for single cell problem\n";
					solver_type = 'i';
				}
				else if(std::string(argv[i]) == "b")
				{
					std::cout << "Brent's solver chosen for single cell problem\n";
					solver_type = 'b';
				}
				else if(std::string(argv[i]) == "r")
				{
					std::cout << "Regula Falsi solver chosen for single cell problem\n";
					solver_type = 'r';
				}
				else
				{
					solver_type = 'r';
				}
			
			}
			else if(std::string(argv[i]) == "-p")
			{
				printIterations = true;
			}
			else if(std::string(argv[i]) == "-t")
			{
				i++;
				time_step_days = std::atof(argv[i]);
			}
			else if(std::string(argv[i]) == "-T")
			{
				i++;
				comp_length_days = std::atof(argv[i]);
			}
			else
				std::cerr << "Invalid argument " << argv[i] << " passed to " << argv[0] << "\n";
		}
	}
		
    /// We define the grid. A Cartesian grid with 400 cells,
    /// each being 10m along each side. Note that we treat the
    /// grid as 3-dimensional, but have a thickness of only one
    /// layer in the Z direction.
    ///
    /// The Opm::GridManager is responsible for creating and destroying the grid,
    /// the UnstructuredGrid data structure contains the actual grid topology
    /// and geometry.
    int nx = 20;
    int ny = 20;
    int nz = 1;
    double dx = 10.0;
    double dy = 10.0;
    double dz = 10.0;
    using namespace Opm;
    GridManager grid_manager(nx, ny, nz, dx, dy, dz);
    const UnstructuredGrid& grid = *grid_manager.c_grid();
    int num_cells = grid.number_of_cells;

    /// We define the properties of the fluid.\n
    /// Number of phases, phase densities, phase viscosities,
    /// rock porosity and permeability.
    ///
    /// We always use SI units in the simulator. Many units are
    /// available for use, however.  They are stored as constants in
    /// the Opm::unit namespace, while prefixes are in the Opm::prefix
    /// namespace. See Units.hpp for more.
    int num_phases = 2;
    using namespace Opm::unit;
    using namespace Opm::prefix;
    std::vector<double> density(num_phases, 1000.0);
    std::vector<double> viscosity(num_phases, 1.0*centi*Poise);
    double porosity = 0.5;
    double permeability = 10.0*milli*darcy;

    /// We define the relative permeability function. We use a basic fluid
    /// description and set this function to be quadratic. For more realistic fluid, the
    /// saturation function may be interpolated from experimental data.
    SaturationPropsBasic::RelPermFunc rel_perm_func = SaturationPropsBasic::Quadratic;

    /// We construct a basic fluid and rock property object
    /// with the properties we have defined above.  Each property is
    /// constant and hold for all cells.
    IncompPropertiesBasic props(num_phases, rel_perm_func, density, viscosity,
                                porosity, permeability, grid.dimensions, num_cells);

    /// Gravity parameters. Here, we set zero gravity.
    const double *grav = 0;
    std::vector<double> omega;
    
    /// We set up the source term. Positive numbers indicate that the cell is a source,
    /// while negative numbers indicate a sink.
    std::vector<double> src(num_cells, 0.0);
    src[0] = 1.;
    src[num_cells-1] = -1.;
    
    /// We set up the boundary conditions. Letting bcs be empty is equivalent
    /// to no-flow boundary conditions.
    FlowBCManager bcs;
    
    /// We may now set up the pressure solver. At this point,
    /// unchanging parameters such as transmissibility are computed
    /// and stored internally by the IncompTpfa class. The null pointer
    /// constructor argument is for wells, which are not used in this tutorial.
    LinearSolverUmfpack linsolver;
    IncompTpfa psolver(grid, props, linsolver, grav, NULL, src, bcs.c_bcs());

    /// We set up a state object for the wells. Here, there are
    /// no wells and we let it remain empty.
    WellState well_state;
    
    /// We compute the pore volume
    std::vector<double> porevol;
    Opm::computePorevolume(grid, props.porosity(), porevol);
    
    /// Set up the transport solver. This is a reordering implicit Euler transport solver.
    const double tolerance = 1e-9;
    const int max_iterations = 30;
    Opm::TransportSolverTwophaseReorder transport_solver(grid, props, NULL, tolerance, max_iterations, solver_type, false);
	
    /// Time integration parameters
    //const double dt = 0.1*day;
    //const int num_time_steps = 20;
    const double comp_length = comp_length_days*day;
    const double dt = time_step_days*day;
    const int num_time_steps = comp_length/dt;
    
    std::cout << "Time step length: " << dt << std::endl;
    
    /// We define a vector which contains all cell indexes. We use this
    /// vector to set up parameters on the whole domain.
    std::vector<int> allcells(num_cells);
    for (int cell = 0; cell < num_cells; ++cell) {
        allcells[cell] = cell;
    }

    /// We set up a two-phase state object, and
    /// initialize water saturation to minimum everywhere.
    TwophaseState state;
    state.init(grid, 2);
    state.setFirstSat(allcells, props, TwophaseState::MinSat);

    /// This string stream will be used to construct a new
    /// output filename at each timestep.
    std::ostringstream vtkfilename;
	
	time::StopWatch clock;
	clock.start();
    for (int i = 0; i < num_time_steps; ++i) {

        psolver.solve(dt, state, well_state);

        transport_solver.solve(&porevol[0], &src[0], dt, state);
        
        if(printIterations)
		{
			printIterationsFromVector(transport_solver, i, num_cells, solver_type, comp_length_days, time_step_days);
			
	        vtkfilename.str("");
	        vtkfilename << "testCase1-s-" << solver_type << "-T-" << replaceStrChar(std::to_string(comp_length),".",'_') << "-t-" << replaceStrChar(std::to_string(dt),".",'_') << "-" << std::setw(3) << std::setfill('0') << i << ".vtu";
	        std::ofstream vtkfile(vtkfilename.str().c_str());
	        Opm::DataMap dm;
	        dm["saturation"] = &state.saturation();
	        dm["pressure"] = &state.pressure();
	        Opm::writeVtkData(grid, dm, vtkfile);
		}
    }
    clock.stop();
    std::cout << "Problem solved in " << clock.secsSinceStart() << " seconds \n";
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}

/// pythonscript3 python script to generate figures:
/// generate_doc_figures.py tutorial3
