#include "CaseUtilities.hpp"

#include "config.h"

//#include <iostream>
#include <iomanip>
#include <fstream>
//#include <vector>
#include <cassert>
#include <opm/core/grid.h>
//#include <opm/core/grid/GridManager.hpp>
#include <opm/core/io/vtk/writeVtkData.hpp>
#include <opm/core/linalg/LinearSolverUmfpack.hpp>
#include <opm/core/pressure/IncompTpfa.hpp>
#include <opm/core/pressure/FlowBCManager.hpp>
#include <opm/core/props/IncompPropertiesBasic.hpp>
#include <opm/core/props/IncompPropertiesShadow.hpp>

//#include <opm/core/transport/reorder/TransportSolverTwophaseReorder.hpp>

//#include <opm/core/simulator/TwophaseState.hpp>
#include <opm/core/simulator/WellState.hpp>

#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/core/utility/StopWatch.hpp>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

namespace Opm
{

using std::string;

void parseArguments(int argc, char ** argv, 
double & muw, double & muo, bool & verbose, bool & solver_flag, 
double & time_step_days, double & comp_length_days, double & xsize, 
double & ysize, int & xdim, int & ydim, char & solver_type, 
bool & printIterations, int & nprint, string & print_points_file_name,
string & perm_file_name, int & layer, double & xpos, double & ypos, 
double & srcVol, double & sinkVol, double & grav_x, double & grav_y, double & grav_z)
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
		else if(std::string(argv[i]) == "-i")
		{
			srcVol = std::atof(argv[++i]);
			sinkVol = -srcVol;
			if(i+1 < argc && !boost::starts_with(std::string(argv[i+1]),"-"))
				sinkVol = std::atof(argv[++i]);
		}
		else if(std::string(argv[i]) == "-g")
		{
			grav_x = std::atof(argv[++i]);
			grav_y = std::atof(argv[++i]);
			grav_z = std::atof(argv[++i]);
		}
		else if(std::string(argv[i]) == "-p")
		{
			printIterations = true;
		}
		else if(std::string(argv[i]) == "--print")
		{
			printIterations = true;
			nprint = std::atoi(argv[++i]);
			
		}
		else
			std::cerr << "Invalid argument " << argv[i] << " passed to " << argv[0] << "\n";
	}
}

void initPrintPointVector(std::vector<int> & print_points, int num_time_steps, int nprint, string print_points_file_name)
{
	int plotInterval = floor((double)num_time_steps/nprint + 0.5);
	print_points.push_back(0);
	std::vector<int> extra_print_points;
	readPrintPointsFromFile(print_points_file_name, extra_print_points);
	int last_interval_point = 0;
	while(print_points.back() < num_time_steps)
	{
		if(!extra_print_points.empty() && extra_print_points.back() < last_interval_point + plotInterval)
		{
			if(extra_print_points.back() != print_points.back())
				print_points.push_back(extra_print_points.back());
			extra_print_points.pop_back();
		}
		else
		{
			last_interval_point += plotInterval;
			if(last_interval_point != print_points.back())
				print_points.push_back(last_interval_point);
		}
	}
	for(std::vector<int>::iterator it = print_points.begin(); it != print_points.end(); it++)
		std::cout << *it << " ";
	std::cout << "\n";
}
bool readPrintPointsFromFile(std::string filename, std::vector<int> & print_points)
{
	std::ifstream printpointfile(filename);
	string line;
	if(printpointfile.is_open())
	{
		while( getline(printpointfile,line) )
		{
			boost::trim(line);
			if(!line.empty() && !boost::starts_with(line,"//"))
			{
				std::vector<string> numbers;
				boost::split(numbers,line,boost::is_any_of("\t "));
				for(int i = numbers.size()-1; i >= 0 ; i--)
				{
					boost::trim(numbers[i]);
					if(!numbers[i].empty())
					{
						int num = atoi(numbers[i].c_str());
						if( (!print_points.empty() && print_points.back() != num) || print_points.empty() )
							print_points.push_back(num);
					}
				}
			}
		}
		printpointfile.close();
		return true;
	}
	else
	{
		std::cout << "Failed to open print point file " << filename << "\n";
		return false;
	}
}

 void printStateDataToVTKFile(string execName, std::ostringstream & vtkfilename, Opm::TwophaseState state, const UnstructuredGrid& grid, bool solver_flag, char solver_type, char extra_solver_char, double comp_length_days, double time_step_days, int i /*, int plotInterval*/)
{
	vtkfilename.str("");
	vtkfilename << execName << "-s-" << solver_type;
	if(solver_flag)
		vtkfilename << extra_solver_char;
	vtkfilename << "-T-" << replaceStrChar(std::to_string(comp_length_days),".",'_') << "-t-" << replaceStrChar(std::to_string(time_step_days),".",'_') << "-" << std::setw(3) << std::setfill('0') <<  i << ".vtu"; //(int)(i/plotInterval) << ".vtu";
	printStateDataToVTKFile(vtkfilename.str(), state, grid);
}

void printStateDataToVTKFile(std::string vtkfilename, Opm::TwophaseState state, const UnstructuredGrid& grid)
{
	std::ofstream vtkfile(vtkfilename.c_str());
	Opm::DataMap dm;
	dm["saturation"] = &state.saturation();
	dm["pressure"] = &state.pressure();
	Opm::writeVtkData(grid, dm, vtkfile);
}

void printIterationsFromVector(string execName, const Opm::TransportSolverTwophaseReorder & transport_solver, int i, int num_cells, const char solver_type, const double comp_length, const double time_step)
{
	std::vector<int> iterations = transport_solver.getReorderIterations();
	std::ostringstream iterfilename;
	
	string str_comp_length = replaceStrChar(std::to_string(comp_length), ".", '_');
	string str_time_step = replaceStrChar(std::to_string(time_step), ".", '_');
	
	// Set filename and open
	iterfilename.str(""); 
	iterfilename << execName << "-iterations-s-" << solver_type << "-T-" << str_comp_length << "-t-" << str_time_step << "-" << std::setw(3) << std::setfill('0') << i << ".txt";
	std::ofstream file; file.open(iterfilename.str().c_str());
	for ( int i = 0; i < num_cells; i++)
	{
		file << i << "\t" << iterations[i] << "\n";
	}
	file.close();
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
			using Opm::unit::darcy;
			using Opm::prefix::milli;
			perm.push_back(Kx_target*milli*darcy); perm.push_back(0); 			perm.push_back(0);
			perm.push_back(0);		   perm.push_back(Ky_target*milli*darcy);	perm.push_back(0);
			perm.push_back(0);		   perm.push_back(0);			perm.push_back(0);
			if(buildCache)
				file << Kx_target*milli*darcy 	<< "\t" << 0 		   				<< "\t" << 0 << "\t"
				     << 0         				<< "\t" << Ky_target*milli*darcy 	<< "\t" << 0 << "\t"
				     << 0         				<< "\t" << 0 		   				<< "\t" << 0 << "\n";
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

}
