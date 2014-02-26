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

#include <opm/core/utility/StopWatch.hpp>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using std::string;

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

void parseArguments(int argc, char ** argv, 
double & muw, double & muo, bool & verbose, bool & solver_flag, 
double & time_step_days, double & comp_length_days, int & xdim, int & ydim, 
char & solver_type, bool & printIterations, string & perm_file_name, int & layer, int & xstart, int & xnum, 
					   int & ystart, int & ynum)
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
				std::cout << "Newton Raphson Trust region solver chosen for single cell problem\n";
				solver_type = 't';
				if(i+1 < argc && !boost::starts_with(std::string(argv[i+1]),"-"))
				{ 
					i++;
					solver_flag = true;
				}
			}
			else if(std::string(argv[i]) == "u")
			{
				std::cout << "Regula Falsi Trust region solver chosen for single cell problem\n";
				solver_type = 'u';
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
		else if(std::string(argv[i]) == "-d")
		{
			i++;
			xdim = std::atof(argv[i]);
			ydim = xdim;
			if(i+1 < argc && !boost::starts_with(std::string(argv[i+1]),"-"))
			{ 
				i++;
			    ydim = std::atoi(argv[i]);
			}
			std::cout << "Using " << xdim << " and " << ydim << " cell(s) in the x- and y-direction, respectively\n";
		}
		else if(std::string(argv[i]) == "-m")
		{
			i++;
			muw = std::atof(argv[i]);
			muo = muw;
			if(i+1 < argc && !boost::starts_with(std::string(argv[i+1]),"-"))
			{ 
				i++;
			    muo = std::atoi(argv[i]);
			}
			std::cout << "Using viscosity " << muw << " and " << muo << " cP for water and oil, respectively\n";
		}
		else if(std::string(argv[i]) == "-v")
		{
			verbose = true;
		}
		else if(std::string(argv[i]) == "-f")
		{
			i++;
			perm_file_name = std::string(argv[i]);
		}
		else if(std::string(argv[i]) == "--perm")
		{
			layer = atof(argv[++i]);
			xstart = atof(argv[++i]);
			xnum = atof(argv[++i]);
			ystart = atof(argv[++i]);
			ynum = atof(argv[++i]);
		}
		else
			std::cerr << "Invalid argument " << argv[i] << " passed to " << argv[0] << "\n";
	}
}

bool readPermData(string perm_file_name, std::vector<double> & Kx, std::vector<double> & Ky, std::vector<double> & Kz)
{
	string line;
	std::ifstream infile(perm_file_name);
	int curr_read_dim = 0; // Counts the current spatial dimension to be read, where K_x == 0, K_y == 1, K_z == 2
	if(infile.is_open())
	{
		while( getline(infile,line) )
		{
			boost::trim(line);
			if(line.empty())
				curr_read_dim++;
			else
			{
				std::vector<string> numbers;
				boost::split(numbers,line,boost::is_any_of("\t "));
				for(unsigned int i = 0;i < numbers.size(); i++)
				{
					boost::trim(numbers[i]);
					if(!numbers[i].empty())
					{
						double num = atof(numbers[i].c_str());
						if(curr_read_dim == 0) // K_x
							Kx.push_back(num);
						else if(curr_read_dim == 1) // K_y
							Ky.push_back(num);
						else if (curr_read_dim == 2) // K_z
							Kz.push_back(num);
						else
						{
							std::cout << "Too many empty lines!\n";
							return false;
						}
					}
				}
			}
		}
		infile.close();
		return true;
	}
	else
	{
		std::cout << "Failed to open permeability file " << perm_file_name << "\n";
		return false;
	}
}

void buildPermMatrixForRegion(std::vector<double> & perm, std::vector<double> Kx, std::vector<double> Ky, int layer, int xstart, int xnum, int ystart, int ynum)
{
	int xdim = 60; int ydim = 220;
	int layerInd = layer*xdim*ydim; // Index of first value in the layer
	int regionRowInd = xdim*ystart; // Index of the first element in the first row containing elements in the region, relative to selected layer
	for(int i = 0; i < ynum; i++) // Iterate over number of cells in y-dir of region
	{
		int colInd = xdim*i + xstart; // Index of first element of current row in region, relative to starting index of first element in the first row with elements in selected region
		for(int j = 0; j < xnum; j++) // Iterate over number of cells in x-dir of region
		{
			int index = layerInd + regionRowInd + colInd + j;
			perm.push_back(Kx[index]);
			perm.push_back(0); perm.push_back(0);
			perm.push_back(Ky[index]);
			//std::cout << Kx[index] << " ";
		}
		//std::cout << "\n";
	}
}

int main (int argc, char ** argv)
try
{
	int layer = 0; 													 // Selected layer
	int xstart = 0; int xnum = 20; int ystart = 0; int ynum = 20; // Selected region in layer
	double muw = 1;
	double muo = 1;
	bool verbose = false;
	bool solver_flag = false;
	double time_step_days = 0.1;
	double comp_length_days = 2;
	int xdim = 20;
	int ydim = 20;
	char solver_type = 'r';
	bool printIterations = false;
	string perm_file_name = "spe_perm.dat";
	// Check parameters
	if(argc > 1)
		parseArguments(argc, argv, muw, muo, verbose, solver_flag, time_step_days, comp_length_days, 
					   xdim, ydim, solver_type, printIterations, perm_file_name, layer, xstart, xnum, 
					   ystart, ynum);
	
	if(verbose)
	{
		std::cout << "----------------- Initializing problem -------------------\n";
		std::cout << "Reading permeabilities from file ...\n";
	}
	
		
	std::vector<double> Kx,Ky,Kz;
	readPermData(perm_file_name,Kx,Ky,Kz);
	
	if(verbose)
	{
		std::cout << "Finished reading permeabilities\n";
		std::cout << "Building permeability table ...\n";
	}
	
	std::vector<double> perm;
	buildPermMatrixForRegion(perm,Kx,Ky,layer,xstart,xnum,ystart,ynum);
	
	if(verbose)
	{
		std::cout << "Finished buildling permeability table\n";
	}
	
	//return 0;
		
    /// The Opm::GridManager is responsible for creating and destroying the grid,
    /// the UnstructuredGrid data structure contains the actual grid topology
    /// and geometry.
    int nx = xdim; // = 20;
    int ny = ydim; // = 20;
    int nz = 1;
    double dx = 10.0;
    double dy = 10.0;
    double dz = 10.0;
    using namespace Opm;
    GridManager grid_manager(nx, ny, nz, dx, dy, dz);
    const UnstructuredGrid& grid = *grid_manager.c_grid();
    int num_cells = grid.number_of_cells;

    int num_phases = 2;
    using namespace Opm::unit;
    using namespace Opm::prefix;
    std::vector<double> density(num_phases, 1000.0);
    double visc_arr[] = {muw*centi*Poise, muo*centi*Poise};
    std::vector<double> viscosity(visc_arr, visc_arr + sizeof(visc_arr)/sizeof(double));
    double porosity = 0.5;
    double permeability = 10.0*milli*darcy;
    SaturationPropsBasic::RelPermFunc rel_perm_func = SaturationPropsBasic::Quadratic;

    IncompPropertiesBasic props(num_phases, rel_perm_func, density, viscosity,
                                porosity, permeability, grid.dimensions, num_cells);

    const double *grav = 0;
    std::vector<double> omega;

    std::vector<double> src(num_cells, 0.0);
    src[0] = 1.;
    src[num_cells-1] = -1.;

    FlowBCManager bcs;

    LinearSolverUmfpack linsolver;
    IncompTpfa psolver(grid, props, linsolver, grav, NULL, src, bcs.c_bcs());

    WellState well_state;
    
    std::vector<double> porevol;
    Opm::computePorevolume(grid, props.porosity(), porevol);
    
    const double tolerance = 1e-9;
    const int max_iterations = 30;
    IncompPropertiesShadow shadow_props(props);
    Opm::TransportSolverTwophaseReorder transport_solver(grid, shadow_props.usePermeability(&perm[0]) , NULL, tolerance, max_iterations, solver_type, verbose, solver_flag);

    const double comp_length = comp_length_days*day;
    const double dt = time_step_days*day;
    const int num_time_steps = comp_length/dt;
    
    std::cout << "Time step length: " << dt << std::endl;

    std::vector<int> allcells(num_cells);
    for (int cell = 0; cell < num_cells; ++cell) {
        allcells[cell] = cell;
    }

    TwophaseState state;
    state.init(grid, 2);
    state.setFirstSat(allcells, props, TwophaseState::MinSat);

    std::ostringstream vtkfilename;
	
	if(verbose)
	{		
		std::cout << "----------------- Solving " << num_time_steps << " time steps -------------------\n";
		std::cout << "Press ENTER to continue computation... " << std::flush;
		std::cin.ignore(std::numeric_limits<std::streamsize> ::max(), '\n');
	}		
	
	time::StopWatch clock;
	clock.start();
    for (int i = 0; i < num_time_steps; ++i) {

        psolver.solve(dt, state, well_state);

        transport_solver.solve(&porevol[0], &src[0], dt, state);
        
        if(printIterations)
		{
			printIterationsFromVector(transport_solver, i, num_cells, solver_type, comp_length_days, time_step_days);
			
	        vtkfilename.str("");
	        vtkfilename << "testCase1-s-" << solver_type << "-T-" << replaceStrChar(std::to_string(comp_length_days),".",'_') << "-t-" << replaceStrChar(std::to_string(time_step_days),".",'_') << "-" << std::setw(3) << std::setfill('0') << i << ".vtu";
	        std::ofstream vtkfile(vtkfilename.str().c_str());
	        Opm::DataMap dm;
	        dm["saturation"] = &state.saturation();
	        dm["pressure"] = &state.pressure();
	        Opm::writeVtkData(grid, dm, vtkfile);
		}
		if(verbose)
			std::cout << "Solved step " << i+1 << " of " << num_time_steps << "\n";
    }
    clock.stop();
    std::cout << "Problem solved in " << clock.secsSinceStart() << " seconds \n";
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}
