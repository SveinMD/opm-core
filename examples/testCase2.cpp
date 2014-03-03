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


string replaceStrChar(string str, const string & replace, char ch);
void printIterationsFromVector(const Opm::TransportSolverTwophaseReorder & transport_solver, 
							   int i, int num_cells, const char solver_type, 
							   const double comp_length, const double time_step);
void parseArguments(int argc, char ** argv, double & muw, double & muo, 
					bool & verbose, bool & solver_flag, 
					double & time_step_days, double & comp_length_days, 
					double & xsize, double & ysize, int & xdim, int & ydim, char & solver_type, 
					bool & printIterations, string & perm_file_name, 
					int & layer, double & xpos, double & ypos);
void constructCacheFileName(std::ostringstream & filename, int layer, 
							double xstart, double xsize, int xnum, 
							double ystart, double ysize, int ynum);
void constructCacheFileName(std::ostringstream & filename, int layer, 
							int xstart, int xnum, int ystart, int ynum);
bool readPermDataFromCache(std::vector<double> & perm, 
						   std::ifstream & cachefile, std::string cachefilename = "");
bool readPermDataFromRawFile(string perm_file_name, std::vector<double> & Kx,
							 std::vector<double> & Ky, std::vector<double> & Kz);
void buildPermMatrixForRegion(std::vector<double> & perm, std::vector<double> Kx, 
							  std::vector<double> Ky, std::vector<double> Kz, int layer, int xstart, 
							  int xnum, int ystart, int ynum, double buildCache);
void buildPermData(string perm_file_name, std::vector<double> & perm, 
				   int layer, int xstart, int xnum, int ystart, int ynum, bool verbose);
void buildPermData(string perm_file_name, std::vector<double> & perm, int layer,
				   double xpos, double ypos, double xsize, double ysize, int nxcells,
				   int nycells, double xsizeperm, double ysizeperm, int nxcellsperm,
				   int nycellsperm, bool verbose);
double interpolate(double fa,double a,double fb,double b,double target);
void interpolatePermData(std::vector<double> & perm, std::vector<double> Kx, std::vector<double> Ky,
						 int layer, double xpos, double ypos, double xsize, double ysize,
						 int nxcells, int nycells, double xsizeperm, double ysizeperm,
						 int nxcellsperm, int nycellsperm, double buildCache);

int main (int argc, char ** argv)
try
{
	int nx = 20; int ny = 20; int nz = 1;
	int layer = 0;
	int nxperm = 60; int nyperm = 220;
	
	double xpos = 0; double ypos = 0;
	double dxperm = 365.76; double dyperm = 670.56;
	double dx = 10.0; double dy = 10.0; double dz = 10.0;
	double muw = 1; double muo = 1;
	double time_step_days = 0.1;
	double comp_length_days = 2;
	
	bool verbose = false;
	bool printIterations = false;
	bool solver_flag = false;
	
	char solver_type = 'r';
	
	string perm_file_name = "spe_perm.dat";
	
	if(argc > 1)
		parseArguments(argc, argv, muw, muo, verbose, solver_flag, time_step_days, comp_length_days, 
					   dx, dy, nx, ny, solver_type, printIterations, perm_file_name, layer, xpos, ypos);
	
	char extra_solver_char = '\0';
	if(solver_flag)
		extra_solver_char = 'a';
			
	if(verbose)
		std::cout << "----------------- Initializing problem -------------------\n";
	
	std::vector<double> perm;
	buildPermData(perm_file_name, perm, layer, xpos, ypos, dx, dy, nx, ny, 
				  dxperm, dyperm, nxperm, nyperm, verbose);
	
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
    IncompPropertiesShadow shadow_props(props);
    
    const double *grav = 0;
    std::vector<double> omega;

    std::vector<double> src(num_cells, 0.0);
    src[0] = 1.;
    src[num_cells-1] = -1.;

    FlowBCManager bcs;

    LinearSolverUmfpack linsolver;
    IncompTpfa psolver(grid, shadow_props.usePermeability(&perm[0]), linsolver, grav, NULL, src, bcs.c_bcs());
	
    WellState well_state;
    
    std::vector<double> porevol;
    Opm::computePorevolume(grid, props.porosity(), porevol);
    
    const double tolerance = 1e-9;
    const int max_iterations = 30;
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
    //state.setFirstSat(allcells, props, TwophaseState::MinSat);
    state.setFirstSat(allcells, shadow_props.usePermeability(&perm[0]), TwophaseState::MinSat);

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
	        vtkfilename << "testCase2-s-" << solver_type;
	        if(solver_flag)
				vtkfilename << extra_solver_char;
	        vtkfilename << "-T-" << replaceStrChar(std::to_string(comp_length_days),".",'_') << "-t-" << replaceStrChar(std::to_string(time_step_days),".",'_') << "-" << std::setw(3) << std::setfill('0') << i << ".vtu";
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
	iterfilename << "testCase2-iterations-s-" << solver_type << "-T-" << str_comp_length << "-t-" << str_time_step << "-" << std::setw(3) << std::setfill('0') << i << ".txt";
	std::ofstream file; file.open(iterfilename.str().c_str());
	for ( int i = 0; i < num_cells; i++)
	{
		file << i << "\t" << iterations[i] << "\n";
	}
	file.close();
}

void parseArguments(int argc, char ** argv, 
double & muw, double & muo, bool & verbose, bool & solver_flag, 
double & time_step_days, double & comp_length_days, double & xsize, double & ysize, int & xdim, int & ydim,
char & solver_type, bool & printIterations, string & perm_file_name, int & layer, double & xpos, double & ypos)
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
			xpos = atof(argv[++i]);
			ypos = atof(argv[++i]);
		}
		else if(std::string(argv[i]) == "--dim")
		{
			xsize = std::atof(argv[++i]);
			ysize = std::atof(argv[++i]);
			xdim = std::atoi(argv[++i]);
			ydim = std::atoi(argv[++i]);
			
			std::cout << "Using " << xdim << "x" << ydim << " cells on a " << xsize << " m by " << ysize << " m domain.\n";
		}
		else
			std::cerr << "Invalid argument " << argv[i] << " passed to " << argv[0] << "\n";
	}
}

void constructCacheFileName(std::ostringstream & filename, int layer, double xstart, double xsize, int xnum, double ystart, double ysize, int ynum)
{
	filename << "permeabilities-z-" << layer << "-Kx-" << replaceStrChar(boost::lexical_cast<std::string>(xstart),"_",'.') << "-" << replaceStrChar(boost::lexical_cast<std::string>(xsize),"_",'.') << "-" << xnum << "-Ky-" << replaceStrChar(boost::lexical_cast<std::string>(ystart),"_",'.') << "-" << replaceStrChar(boost::lexical_cast<std::string>(ysize),"_",'.') << "-" << ynum << ".cache";
}

void constructCacheFileName(std::ostringstream & filename, int layer, int xstart, int xnum, int ystart, int ynum)
{
	filename << "permeabilities-z-" << layer << "-Kx-" << xstart << "-" << xnum << "-Ky-" << ystart << "-" << ynum << ".cache";
}

bool readPermDataFromCache(std::vector<double> & perm, std::ifstream & cachefile, std::string cachefilename)
{
	string line;
	if(cachefile.is_open())
	{
		while( getline(cachefile,line) )
		{
			boost::trim(line);
			if(!line.empty())
			{
				std::vector<string> numbers;
				boost::split(numbers,line,boost::is_any_of("\t "));
				for(unsigned int i = 0;i < numbers.size(); i++)
				{
					boost::trim(numbers[i]);
					if(!numbers[i].empty())
					{
						double num = atof(numbers[i].c_str());
						perm.push_back(num);
					}
				}
			}
		}
		cachefile.close();
		return true;
	}
	else
	{
		std::cout << "Failed to open permeability cache file " << cachefilename << "\n";
		return false;
	}
}

bool readPermDataFromRawFile(string perm_file_name, std::vector<double> & Kx, std::vector<double> & Ky, std::vector<double> & Kz)
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

void buildPermMatrixForRegion(std::vector<double> & perm, std::vector<double> Kx, std::vector<double> Ky, std::vector<double> Kz, int layer, int xstart, int xnum, int ystart, int ynum, double buildCache)
{
	std::ofstream file;
	if(buildCache)
	{
		std::ostringstream filename;
		constructCacheFileName(filename,layer,xstart,xnum,ystart,ynum);
		file.open(filename.str().c_str());
	}
	int xdim = 60; int ydim = 220;
	int layerInd = layer*xdim*ydim; // Index of first value in the layer
	int regionRowInd = xdim*ystart; // Index of the first element in the first row containing elements in the region, relative to selected layer
	for(int i = 0; i < ynum; i++) // Iterate over number of cells in y-dir of region
	{
		int colInd = xdim*i + xstart; // Index of first element of current row in region, relative to starting index of first element in the first row with elements in selected region
		for(int j = 0; j < xnum; j++) // Iterate over number of cells in x-dir of region
		{
			int index = layerInd + regionRowInd + colInd + j;
			
			//perm.push_back(Kx[index]);
			//perm.push_back(0); perm.push_back(0);
			//perm.push_back(Ky[index]);
			
			perm.push_back(Kx[index]); perm.push_back(0); 			perm.push_back(0);
			perm.push_back(0);		   perm.push_back(Ky[index]);	perm.push_back(0);
			perm.push_back(0);		   perm.push_back(0);			perm.push_back(Kz[index]);
			
			if(buildCache)
				file << Kx[index] << "\t" << 0 		   << "\t" << 0 		<< "\t"
				     << 0         << "\t" << Ky[index] << "\t" << 0 		<< "\t"
				     << 0         << "\t" << 0 		   << "\t" << Kz[index] << "\n";
				//file << Kx[index] << "\t" << Ky[index] << "\t" << Kz[index] << "\n";
				//file << Kx[index] << "\t" << 0 << "\t" << 0 << "\t" << Ky[index] << "\n";
		}
	}
	file.close();
}

void buildPermData(string perm_file_name, std::vector<double> & perm, 
				   int layer, int xstart, int xnum, int ystart, int ynum, bool verbose)
{
	if(verbose)
		std::cout << "Initializing permeability table ... \n";
	std::ostringstream cachefilename;
	cachefilename.str("");
	constructCacheFileName(cachefilename,layer,xstart,xnum,ystart,ynum);
	std::ifstream cachefile; cachefile.open(cachefilename.str().c_str());
	if(cachefile.is_open())
	{
		if(verbose)
			std::cout << "Reading permeabilities from cache " << cachefilename.str() << " ... \n";
		
		readPermDataFromCache(perm, cachefile, cachefilename.str());
	}
	else
	{
		if(verbose)
			std::cout << "Reading permeabilities from file ...\n";
			
		std::vector<double> Kx,Ky,Kz;
		readPermDataFromRawFile(perm_file_name,Kx,Ky,Kz);
		
		if(verbose)
		{
			std::cout << "Finished reading permeabilities\n";
			std::cout << "Building permeability table ...\n";
		}
		
		buildPermMatrixForRegion(perm,Kx,Ky,Kz,layer,xstart,xnum,ystart,ynum,true);
	}
	
	if(verbose)
		std::cout << "Permeability table initialized. \n";
}

void buildPermData(string perm_file_name, std::vector<double> & perm, int layer,
				   double xpos, double ypos, double xsize, double ysize, int nxcells,
				   int nycells, double xsizeperm, double ysizeperm, int nxcellsperm,
				   int nycellsperm, bool verbose)
{
	if(verbose)
		std::cout << "Initializing permeability table ... \n";
	std::ostringstream cachefilename;
	cachefilename.str("");
	constructCacheFileName(cachefilename,layer,xpos,xsize,nxcells,ypos,ysize,nycells);
	std::ifstream cachefile; cachefile.open(cachefilename.str().c_str());
	if(cachefile.is_open())
	{
		if(verbose)
			std::cout << "Reading permeabilities from cache " << cachefilename.str() << " ... \n";
		readPermDataFromCache(perm, cachefile, cachefilename.str());
	}
	else
	{
		if(verbose)
			std::cout << "Reading permeabilities from file ...\n";
			
		std::vector<double> Kx,Ky,Kz;
		readPermDataFromRawFile(perm_file_name,Kx,Ky,Kz);
		
		if(verbose)
		{
			std::cout << "Finished reading permeabilities\n";
			std::cout << "Building permeability table ...\n";
		}
		
		interpolatePermData(perm, Kx, Ky, layer, xpos, ypos, xsize, ysize, nxcells,
							nycells, xsizeperm, ysizeperm, nxcellsperm, nycellsperm, true);
	}
	
	if(verbose)
		std::cout << "Permeability table initialized. \n";
}	

void interpolatePermData(std::vector<double> & perm, std::vector<double> Kx, std::vector<double> Ky,
						 int layer, double xpos, double ypos, double xsize, double ysize,
						 int nxcells, int nycells, double xsizeperm, double ysizeperm,
						 int nxcellsperm, int nycellsperm, double buildCache)
{
	std::ofstream file;
	if(buildCache)
	{
		std::ostringstream filename;
		constructCacheFileName(filename,layer,xpos,xsize,nxcells,ypos,ysize,nycells);
		file.open(filename.str().c_str());
	}
	
	//double xpos,ypos,xsize,ysize; // Physical dimensions of computational domain
	//int nxcells, nycells;		  // Number of cells in computational domain
	double dxcell=xsize/nxcells;  // Size of computation cell, x direction
	double dycell=ysize/nycells;  // Size of computation cell, y direction
	
	//double xsizeperm, ysizeperm;  // Physical dimensions of permeability data domain
	//int nxcellsperm, nycellsperm; // Number of cells in permeability data
	double dxcellperm = xsizeperm/nxcellsperm; // Size of cell in permeability data domain, x-direction
	double dycellperm = ysizeperm/nycellsperm; // size of cell in permeability data domain, y direction
	
	int indyposraw = floor(ypos/dycellperm);
	int indxposraw = floor(xpos/dxcellperm);
	int layerIndRaw = layer*nxcellsperm*nycellsperm;
	int regionRowIndRaw = nxcellsperm*indyposraw;
	
	if(xsize > xsizeperm || ysize > ysizeperm)
		std::cout << "Size of computational domain exceeds size of permeability data set!\n"; // TODO: Throw some size error
	
	std::vector<double> permDataRaw;  // Raw data from the permeability data domain
	std::vector<double> permDataComp; // Permeability data for the computational domain
	for(int j = 0; j < nycells; j++)
	{
		double ycell = ypos + j*dycell + dycell*0.5; // y-coordinate of computational cell centre
		int yindpermcell = floor(ycell/dycellperm);
		int colIndRaw = nxcellsperm*yindpermcell + indxposraw;
		for(int i = 0; i < nxcells; i++)
		{
			double xcell = xpos + i*dxcell + dxcell*0.5; // x-coordinate of computational cell centre
			int xindpermcell = floor(xcell/dxcellperm); // Index of cell in permeability domain
			
			int cellIndexRaw = layerIndRaw + regionRowIndRaw + colIndRaw + xindpermcell;
			
			bool printStuff = false;
			
			double Kx_target; double Ky_target;
			if( xcell <= dxcellperm*0.5 || xcell >= (xsizeperm - dxcellperm*0.5) || ycell <= dycellperm*0.5 || ycell >= (ysizeperm - dycellperm*0.5) )
			{
				if(printStuff)
				{
					std::cout << "** Edge case **\n";
					std::cout << "(" << i << "," << j << "),(" << xcell << "," << ycell << "),(" << dxcell << "," << dycell << "),(" << dxcellperm << "," << dycellperm << "),(" << xindpermcell << "," << yindpermcell << ")\n";
				}
				Kx_target = Kx[cellIndexRaw];
				Ky_target = Ky[cellIndexRaw];
				
				if(printStuff)
				{
					std::cout << Kx_target << "\n";
					std::cout << "** End edge case **\n";
				}
			}
			else // Interpolate
			{
				if(printStuff)
					std::cout << "** Inner case **\n";
				
				double y1 = yindpermcell*dycellperm+0.5*dycellperm;
				int ysign = (ycell > y1) ? 1 : -1;
				
				double x1 = xindpermcell*dxcellperm+0.5*dxcellperm;
				int xsign = (xcell > x1) ? 1 : -1;
				
				double y2 = (yindpermcell + ysign)*dycellperm+0.5*dycellperm;
				double x2 = (xindpermcell + xsign)*dxcellperm+0.5*dxcellperm;
				
				int p1 = cellIndexRaw; // Base cell
				int p2 = p1 + xsign; // Next cell in x-direction
				int p3 = p1 + ysign*nxcellsperm; // Next cell in y-direction
				int p4 = p3 + xsign; // Next cell in both x- and y-direction
				
				if(printStuff)
				{
					std::cout << "(" << i << "," << j << "),(" << xcell << "," << ycell << "),(" << xindpermcell << "," << yindpermcell << "),(" << p1 << ", " << p2 << ", " << p3 << ", " << p4 << ")\n";
					std::cout << "(" << y1 << "," << y2 << "),(" << Kx[p1] << "," << Kx[p3] << ")\n";
				}
				double Kx_y1 = interpolate(Kx[p1],y1,Kx[p3],y2,ycell);
				
				if(printStuff)
					std::cout << "(" << y1 << "," << y2 << "),(" << Kx[p2] << "," << Kx[p4] << ")\n";
				double Kx_y2 = interpolate(Kx[p2],y1,Kx[p4],y2,ycell);
				
				if(printStuff)
					std::cout << "(" << x1 << "," << x2 << "),(" << Kx_y1 << "," << Kx_y2 << ")\n";
				Kx_target = interpolate(Kx_y1,x1,Kx_y2,x2,xcell);
				
				double Ky_y1 = interpolate(Ky[p1],y1,Ky[p3],y2,ycell);
				double Ky_y2 = interpolate(Ky[p2],y1,Ky[p4],y2,ycell);
				Ky_target = interpolate(Ky_y1,x1,Ky_y2,x2,xcell);
				
				if(printStuff)
				{
					std::cout << "(" << xcell << "," << Kx_target << ")\n";
					std::cout << "** End inner case **\n";
				}
			}
			
			perm.push_back(Kx_target); perm.push_back(0); 			perm.push_back(0);
			perm.push_back(0);		   perm.push_back(Ky_target);	perm.push_back(0);
			perm.push_back(0);		   perm.push_back(0);			perm.push_back(0);
			if(buildCache)
				file << Kx_target << "\t" << 0 		   << "\t" << 0 << "\t"
				     << 0         << "\t" << Ky_target << "\t" << 0 << "\t"
				     << 0         << "\t" << 0 		   << "\t" << 0 << "\n";
				//file << Kx_target << "\t" << 0 << "\t" << 0 << "\t" << Ky_target << "\n";
		}
	}
	file.close();
}

double interpolate(double fa,double a,double fb,double b,double target)
{
	if( target < std::min(a,b) || target > std::max(a,b) )
		std::cout << "Interpolation error!\n"; // TODO: Throw some range error
	if(b == a)
		return fa;
	double ftarget = fa + (fb-fa)/(b-a)*(target-a);
	if(ftarget < 0)
		std::cout << "Negative permeability using (" << fa << "," << a << "),(" << fb << "," << b << ") with target " << target << "\n";
	return ftarget;
}
