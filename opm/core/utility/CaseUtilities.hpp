#include "config.h"

#include <iostream>
//#include <iomanip>
//#include <fstream>
#include <vector>
//#include <cassert>
#include <opm/core/grid.h>
//#include <opm/core/grid/GridManager.hpp>
//#include <opm/core/io/vtk/writeVtkData.hpp>
//#include <opm/core/linalg/LinearSolverUmfpack.hpp>
//#include <opm/core/pressure/IncompTpfa.hpp>
//#include <opm/core/pressure/FlowBCManager.hpp>
//#include <opm/core/props/IncompPropertiesBasic.hpp>
//#include <opm/core/props/IncompPropertiesShadow.hpp>

#include <opm/core/transport/reorder/TransportSolverTwophaseReorder.hpp>

#include <opm/core/simulator/TwophaseState.hpp>
//#include <opm/core/simulator/WellState.hpp>

//#include <opm/core/utility/miscUtilities.hpp>
//#include <opm/core/utility/Units.hpp>
//#include <opm/core/utility/parameters/ParameterGroup.hpp>

//#include <opm/core/utility/StopWatch.hpp>

//#include <boost/algorithm/string/predicate.hpp>
//#include <boost/algorithm/string.hpp>
//#include <boost/lexical_cast.hpp>

namespace Opm
{

using std::string;

void initPrintPointVector(std::vector<int> & print_points, int num_time_steps, int nprint, bool printIterations, string print_points_file_name);
void printStateDataToVTKFile(std::ostringstream & vtkfilename, Opm::TwophaseState state, 
						 const UnstructuredGrid& grid, bool solver_flag, char solver_type, char extra_solver_char, 
						 double comp_length_days, double time_step_days, int i /*, int plotInterval*/);
void printStateDataToVTKFile(string vtkfilename, Opm::TwophaseState state, const UnstructuredGrid& grid);						 
string replaceStrChar(string str, const string & replace, char ch);
void printIterationsFromVector(const Opm::TransportSolverTwophaseReorder & transport_solver, 
							   int i, int num_cells, const char solver_type, 
							   const double comp_length, const double time_step);
void parseArguments(int argc, char ** argv, double & muw, double & muo, 
					bool & verbose, bool & solver_flag, double & time_step_days, double & comp_length_days,
					double & xsize, double & ysize, int & xdim, int & ydim, char & solver_type, 
					bool & printIterations, int & nprint, string & print_points_file_name,
					string & perm_file_name, int & layer, double & xpos, double & ypos, 
					double & srcVol, double & sinkVol, double & grav_x, double & grav_y, double & grav_z);
void constructCacheFileName(std::ostringstream & filename, int layer, 
							double xstart, double xsize, int xnum, 
							double ystart, double ysize, int ynum);
void constructCacheFileName(std::ostringstream & filename, int layer, 
							int xstart, int xnum, int ystart, int ynum);
bool readPrintPointsFromFile(std::string filename, std::vector<int> & print_points);
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
						 
}
