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
#include <opm/core/props/IncompPropertiesShadow.hpp>
#include <opm/core/props/IncompPropertiesInterface.hpp>

#include <opm/core/transport/reorder/TransportSolverTwophaseReorder.hpp>

#include <opm/core/simulator/TwophaseState.hpp>
#include <opm/core/simulator/WellState.hpp>

#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/CaseUtilities.hpp>

#include <opm/core/utility/StopWatch.hpp>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>

using std::string;

int main (int argc, char ** argv)
try
{
	int nx = 20; int ny = 20; int nz = 1;
	int layer = 0;
	const int NPRINT = 100;
	int nprint = NPRINT;
	
	double xpos = 0; double ypos = 0;
	double dx = 10.0; double dy = 10.0; double dz = 10.0;
	double perm_mD = 10;
	double muw = 1; double muo = 1;
	double time_step_days = 0.1;
	double comp_length_days = 2;
	double srcVol = 0.2;
	double sinkVol = -srcVol;
	double grav_x = 0;
	double grav_y = 0;
	double grav_z = 0;
	double tol = 1e-9;
	double ddummy = 0.0;
	
	bool verbose = false;
	bool printFluxData = false;
	bool useInitialGuessApproximation = false;
	bool printVTU = false;
	bool printIterations = false;
	bool printIterOrVtu = false;
	bool time_pressure_solver = true;
	bool bdummy;
	
	Opm::RootFinderType solver_type = Opm::RegulaFalsiType;
	
	string perm_file_name = "spe_perm.dat";
	string print_points_file_name = "print_points.dat";
	string execName = boost::filesystem::path(std::string(argv[0])).stem().string();
	
	using namespace Opm;
	
	if(argc > 1)
		parseArguments(argc, argv, muw, muo, verbose, time_step_days, comp_length_days, 
					   dx, dy, dz, nx, ny, nz, solver_type, time_pressure_solver, printVTU, printIterations, nprint, 
					   print_points_file_name, perm_file_name, layer, xpos, ypos, perm_mD, bdummy,
					   srcVol, sinkVol, grav_x, grav_y, grav_z, tol, bdummy, bdummy, useInitialGuessApproximation, printFluxData,ddummy,ddummy);
	
	printIterOrVtu = printIterations || printVTU;
	
	if(verbose)
		std::cout << "----------------- Initializing problem -------------------\n";
	
    GridManager grid_manager(nx, ny, nz, dx, dy, dz);
    const UnstructuredGrid& grid = *grid_manager.c_grid();
    int num_cells = grid.number_of_cells;
    
    //std::cout << grid.cell_volumes[0] << std::endl;
    
    int num_phases = 2;
    using namespace Opm::unit;
    using namespace Opm::prefix;
    std::vector<double> density(num_phases, 1000.0);
    density[1] = 800.0;
    double visc_arr[] = {muw*centi*Poise, muo*centi*Poise};
    std::vector<double> viscosity(visc_arr, visc_arr + sizeof(visc_arr)/sizeof(double));
    double porosity = 0.5;
    double permeability = perm_mD*milli*darcy;
    SaturationPropsBasic::RelPermFunc rel_perm_func = SaturationPropsBasic::Quadratic;
    
    IncompPropertiesBasic props(num_phases, rel_perm_func, density, viscosity, porosity, permeability, grid.dimensions, num_cells);
    IncompPropertiesShadow shadow_props(props);
    
    const double *grav = 0;
    std::vector<double> omega;
    
    std::vector<double> porevol;
    Opm::computePorevolume(grid, props.porosity(), porevol);

	double injectedFluidAbsolute = srcVol; // m^3
	double poreVolume = porevol[0];
	double injectedFluidPoreVol = injectedFluidAbsolute/poreVolume;
	
    std::vector<double> src(num_cells, 0.0);
    src[0] = injectedFluidPoreVol; //1.;
    src[num_cells-1] = -injectedFluidPoreVol; //-1.;

    FlowBCManager bcs;

    LinearSolverUmfpack linsolver;
    IncompPropertiesInterface * prop_pointer = (IncompPropertiesInterface *) &props;
	IncompTpfa psolver(grid, *prop_pointer, linsolver, grav, NULL, src, bcs.c_bcs());
	
    WellState well_state;
    
    const double tolerance = tol;
    const int max_iterations = 50;
	Opm::TransportSolverTwophaseReorder transport_solver(grid, *prop_pointer, grav, tolerance, max_iterations, solver_type, verbose, useInitialGuessApproximation, printFluxData);

    const double comp_length = comp_length_days*day;
    const double dt = time_step_days*day;
    const int num_time_steps = comp_length/dt;
    nprint = std::min(nprint,num_time_steps);
    std::cout << "Time step length: " << dt << "\n";
	
	TwophaseState state;
    state.init(grid, 2);
    
    std::vector<int> allcells(num_cells);
    for (int cell = 0; cell < num_cells; ++cell) {
        allcells[cell] = cell;
    }
	state.setFirstSat(allcells, *prop_pointer, TwophaseState::MinSat);

    std::ostringstream vtkfilename;
	
	std::vector<int> print_points;
	if(printVTU)
	{
		if(nprint == NPRINT)
			nprint = num_time_steps;
		initPrintPointVector(print_points, num_time_steps, nprint, print_points_file_name);
	}
	
	if(verbose)
	{		
		std::cout << "----------------- Solving " << num_time_steps << " time steps -------------------\n";
		std::cout << "Press ENTER to continue computation... " << std::flush;
		std::cin.ignore(std::numeric_limits<std::streamsize> ::max(), '\n');
	}		
	
	std::vector<int>::iterator it = print_points.begin();
	time::StopWatch clock;
	double cputime = 0;
	if(time_pressure_solver)
		clock.start();
    for (int i = 0; i < num_time_steps; ++i) {
		if(verbose)
			std::cout << "*** Solving step " << i+1 << " of " << num_time_steps << " ***\n";
			
		if(verbose)
			std::cout << "Solving pressure system ...\n";
        psolver.solve(dt, state, well_state);
        
        if(verbose)
			std::cout << "Solving transport system:\n";
		if(!time_pressure_solver)
			clock.start();
        transport_solver.solve(&porevol[0], &src[0], dt, state);
        if(!time_pressure_solver)
		{
			clock.stop();
			cputime += clock.secsSinceStart();
		}
		
        if(printIterOrVtu)
		{
			if(printIterations)
				printIterationsFromVector(execName, transport_solver, i, num_cells, solver_type, comp_length_days, time_step_days, viscosity[0]/viscosity[1]);
			if(printVTU && it != print_points.end() && *it == i)
			{
				it++;
				printStateDataToVTKFile(execName, vtkfilename, state, grid, solver_type, comp_length_days, time_step_days, i);
			}
		}
    }
    if(time_pressure_solver)
    {
		clock.stop();
		cputime = clock.secsSinceStart();
	}
    std::cout << "Problem solved in " << cputime << " seconds \n";
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}
