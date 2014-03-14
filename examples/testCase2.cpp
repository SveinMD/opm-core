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
	int nxperm = 60; int nyperm = 220;
	int nprint = 100;
	
	double xpos = 0; double ypos = 0;
	double dxperm = 365.76; double dyperm = 670.56;
	double dx = 10.0; double dy = 10.0; double dz = 10.0;
	double muw = 1; double muo = 1;
	double time_step_days = 0.1;
	double comp_length_days = 2;
	double srcVol = 0.2;
	double sinkVol = -srcVol;
	double grav_x = 0;
	double grav_y = 0;
	double grav_z = 0;
	
	bool verbose = false;
	bool printIterations = false;
	bool solver_flag = false;
	
	char solver_type = 'r';
	
	string perm_file_name = "spe_perm.dat";
	string print_points_file_name = "print_points.dat";
	string execName = boost::filesystem::path(std::string(argv[0])).stem().string();
	
	using namespace Opm;
	
	if(argc > 1)
		parseArguments(argc, argv, muw, muo, verbose, solver_flag, time_step_days, comp_length_days, 
					   dx, dy, nx, ny, solver_type, printIterations, nprint, 
					   print_points_file_name, perm_file_name, layer, xpos, ypos,
					   srcVol, sinkVol, grav_x, grav_y, grav_z);
	
	char extra_solver_char = '\0';
	if(solver_flag)
		extra_solver_char = 'a';
			
	if(verbose)
		std::cout << "----------------- Initializing problem -------------------\n";
	
	std::vector<double> perm;
	buildPermData(perm_file_name, perm, layer, xpos, ypos, dx, dy, nx, ny, 
				  dxperm, dyperm, nxperm, nyperm, verbose);
	
    GridManager grid_manager(nx, ny, nz, dx, dy, dz);
    const UnstructuredGrid& grid = *grid_manager.c_grid();
    int num_cells = grid.number_of_cells;
    
    int num_phases = 2;
    using namespace Opm::unit;
    using namespace Opm::prefix;
    std::vector<double> density(num_phases, 1000.0);
    density[1] = 800.0;
    double visc_arr[] = {muw*centi*Poise, muo*centi*Poise};
    std::vector<double> viscosity(visc_arr, visc_arr + sizeof(visc_arr)/sizeof(double));
    double porosity = 0.5;
    double permeability = 10.0*milli*darcy;
    SaturationPropsBasic::RelPermFunc rel_perm_func = SaturationPropsBasic::Quadratic;

    IncompPropertiesBasic props(num_phases, rel_perm_func, density, viscosity,
                                porosity, permeability, grid.dimensions, num_cells);
    IncompPropertiesShadow shadow_props(props);
    
    const double grav_arr [] = {grav_x, grav_y, grav_z};
    const double *grav = &grav_arr[0]; //0;
    std::vector<double> omega;

	double injectedFluidAbsolute = srcVol; // m^3
	double poreVolume = dz*dx*dy*porosity/(nx*ny);
	double injectedFluidPoreVol = injectedFluidAbsolute/poreVolume;
	
    std::vector<double> src(num_cells, 0.0);
    src[0] = injectedFluidPoreVol; //1.;
    src[num_cells-1] = -injectedFluidPoreVol; //-1.;

    FlowBCManager bcs;

    LinearSolverUmfpack linsolver;
    IncompTpfa psolver(grid, shadow_props.usePermeability(&perm[0]), linsolver, grav, NULL, src, bcs.c_bcs());
	
    WellState well_state;
    
    std::vector<double> porevol;
    Opm::computePorevolume(grid, props.porosity(), porevol);
    
    const double tolerance = 1e-9;
    const int max_iterations = 30;
    Opm::TransportSolverTwophaseReorder transport_solver(grid, shadow_props.usePermeability(&perm[0]), grav, tolerance, max_iterations, solver_type, verbose, solver_flag);

    const double comp_length = comp_length_days*day;
    const double dt = time_step_days*day;
    const int num_time_steps = comp_length/dt;
    nprint = std::min(nprint,num_time_steps);
    std::cout << "Time step length: " << dt << std::endl;
	
	TwophaseState state;
    state.init(grid, 2);
    
    std::vector<int> allcells(num_cells);
    for (int cell = 0; cell < num_cells; ++cell) {
        allcells[cell] = cell;
    }
    state.setFirstSat(allcells, shadow_props.usePermeability(&perm[0]), TwophaseState::MinSat);
    
    /*unsigned int nxl = floor(nx*0.5);
    unsigned int nxr = nx - nxl;
    std::vector<int> lefthalf(nxr*ny);
    for(unsigned int cell = 0; cell < nxl*ny; ++cell)
    {
		int iy = (int)(cell/nxl);
		lefthalf[cell] = iy*nx + cell - iy*nxl;
	}
    std::vector<int> righthalf(nxr*ny);
	for(unsigned int cell = 0; cell < nxr*ny; cell++)
	{
		int iy = (int)(cell/nxl);
		righthalf[cell] = iy*nx + cell - iy*nxl + nxl;
	}
	state.setFirstSat(lefthalf, shadow_props.usePermeability(&perm[0]), TwophaseState::MaxSat);
	state.setFirstSat(righthalf, shadow_props.usePermeability(&perm[0]), TwophaseState::MinSat);*/
	
    std::ostringstream vtkfilename;
	
	std::vector<int> print_points;
	if(printIterations)
		initPrintPointVector(print_points, num_time_steps, nprint, print_points_file_name);
	
	if(verbose)
	{		
		std::cout << "----------------- Solving " << num_time_steps << " time steps -------------------\n";
		std::cout << "Press ENTER to continue computation... " << std::flush;
		std::cin.ignore(std::numeric_limits<std::streamsize> ::max(), '\n');
	}		
	
	std::vector<int>::iterator it = print_points.begin();
	time::StopWatch clock;
	clock.start();
    for (int i = 0; i < num_time_steps; ++i) {
		if(verbose)
			std::cout << "Solving pressure system ...\n";
        psolver.solve(dt, state, well_state);
        if(verbose)
			std::cout << "Solving transport system:\n";
        transport_solver.solve(&porevol[0], &src[0], dt, state);
        //transport_solver.solveGravity(&porevol[0], dt, state);
        
        if(printIterations && it != print_points.end() && *it == i) //( (i % plotInterval) == 0 ) )
		{
			it++;
			printIterationsFromVector(execName, transport_solver, i, num_cells, solver_type, comp_length_days, time_step_days);
			printStateDataToVTKFile(execName, vtkfilename, state, grid, solver_flag, solver_type, extra_solver_char, comp_length_days, time_step_days, *it /*i , plotInterval*/);
		}
		if(verbose)
			std::cout << "* Solved step " << i+1 << " of " << num_time_steps << " *\n";
    }
    clock.stop();
    std::cout << "Problem solved in " << clock.secsSinceStart() << " seconds \n";
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}
