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
	const int NPRINT = 100;
	int nprint = NPRINT;
	int zpos = 0;
	
	double xpos_double = 0.0, ypos_double = 0.0;
	double dx = 10.0; double dy = 10.0; double dz = 10.0;
	double muw = 1; double muo = 1;
	double perm_mD = 10;
	double time_step_days = 0.1;
	double comp_length_days = 2;
	double srcVol = 0.2;
	double sinkVol = -srcVol;
	double grav_x = 0;
	double grav_y = 0;
	double grav_z = 0;
	double tol = 1e-9;
	double denswater = 1000.0, densoil = 800.0;
	
	bool verbose = false;
	bool solve_gravity_column = false;
	bool useInitialGuessApproximation = false;
	bool printFluxValues = false;
	bool printVTU = false;
	bool printIterations = false;
	bool printIterOrVtu = false;
	bool time_pressure_solver = true;
	bool initBottomTop = false;
	bool initLeftRight = false;
	
	Opm::RootFinderType solver_type = Opm::RegulaFalsiType;
	
	string perm_file_name = "";
	string print_points_file_name = "print_points.dat";
	string execName = boost::filesystem::path(std::string(argv[0])).stem().string();
	
	using namespace Opm;
	
	bool bdummy;
	if(argc > 1)
		parseArguments(argc, argv, muw, muo, verbose, time_step_days, comp_length_days, 
					   dx, dy, dz, nx, ny, nz, solver_type, time_pressure_solver, printVTU, printIterations, nprint, 
					   print_points_file_name, perm_file_name, zpos, xpos_double, ypos_double, perm_mD, bdummy,
					   srcVol, sinkVol, grav_x, grav_y, grav_z, tol, initBottomTop, initLeftRight, useInitialGuessApproximation, printFluxValues,denswater, densoil);
	
	printIterOrVtu = printIterations || printVTU;
	
	if(verbose)
		std::cout << "----------------- Initializing problem -------------------\n";
	
    GridManager grid_manager(nx, ny, nz, dx, dy, dz);
    const UnstructuredGrid& grid = *grid_manager.c_grid();
    int num_cells = grid.number_of_cells;

    int num_phases = 2;
    using namespace Opm::unit;
    using namespace Opm::prefix;
    std::vector<double> density(num_phases, denswater);
    density[1] = densoil;
    double visc_arr[] = {muw*centi*Poise, muo*centi*Poise};
    std::vector<double> viscosity(visc_arr, visc_arr + sizeof(visc_arr)/sizeof(double));
    double porosity = 0.5;
    double permeability = perm_mD*milli*darcy;
    SaturationPropsBasic::RelPermFunc rel_perm_func = SaturationPropsBasic::Quadratic;
    
    IncompPropertiesBasic props(num_phases, rel_perm_func, density, viscosity, porosity, permeability, grid.dimensions, num_cells);
    IncompPropertiesShadow shadow_props(props);
    
    const double grav_arr [] = {grav_x, grav_y, grav_z};
    const double *grav = &grav_arr[0];
    solve_gravity_column = ( fabs( density[1]-density[0] ) > 0.0 ) && ( fabs(grav_x)+fabs(grav_y)+fabs(grav_z) > 0.0 );
    std::vector<double> omega;
    
    std::vector<double> porevol;
    Opm::computePorevolume(grid, props.porosity(), porevol);

	double injectedFluidAbsolute = srcVol;
	double injectedFluidPoreVol = injectedFluidAbsolute/porevol[0];
    std::vector<double> src(num_cells, 0.0);
    src[0] = injectedFluidPoreVol;
    src[num_cells-1] = -injectedFluidPoreVol;
	/*for(int l = 0; l < nz; l++)
	{
		int cell_src = l*nx*ny;
		int cell_sink = cell_src + nx*ny - 1;
		src[cell_src] = injectedFluidPoreVol;
		src[cell_sink] = -injectedFluidPoreVol;
	}*/

    FlowBCManager bcs;
	
	LinearSolverUmfpack linsolver;
    IncompPropertiesInterface * prop_pointer = (IncompPropertiesInterface *) &props;
	
	IncompTpfa psolver(grid, *prop_pointer, linsolver, grav, NULL, src, bcs.c_bcs());
	
    WellState well_state;
    
    const double tolerance = tol;
    const int max_iterations = 50;
	Opm::TransportSolverTwophaseReorder transport_solver(grid, *prop_pointer, grav, tolerance, max_iterations, solver_type, verbose, useInitialGuessApproximation, printFluxValues);
	
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

	if(initBottomTop)
	{
		unsigned int nzt = floor(nz*0.5);
	    unsigned int nzb = nz - nzt;
	    std::vector<int> tophalf(nzt*nx*ny);
	    for(unsigned int cell = 0; cell < nzt*nx*ny; ++cell)
	    {
			tophalf[cell] = cell;
		}
	    std::vector<int> bottomhalf(nzb*nx*ny);
		for(unsigned int cell = 0; cell < nzb*nx*ny; cell++)
		{
			bottomhalf[cell] = nzt*nx*ny + cell;
		}
		state.setFirstSat(tophalf, *prop_pointer, TwophaseState::MaxSat);
		state.setFirstSat(bottomhalf, *prop_pointer, TwophaseState::MinSat);	
		std::cout << "Init. sat done" << std::endl;
	}
	else if(initLeftRight)
	{
	    unsigned int nxl = floor(nx*0.5);
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
		state.setFirstSat(lefthalf, *prop_pointer, TwophaseState::MaxSat);
		state.setFirstSat(righthalf, *prop_pointer, TwophaseState::MinSat);
	}

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
		
        if(solve_gravity_column)
        {
			if(verbose)
				std::cout << "Solving gravity transport system:\n";
			transport_solver.solveGravity(&porevol[0], dt, state);
		}
		
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
