#include "config.h"

#include <cstring>
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

#include <boost/bimap.hpp>

namespace Opm
{
typedef boost::bimap<RootFinderType,std::string> rfs_map; // (r)oot (f)inder to (s)tring bidirectional map
using std::string;

rfs_map getRootFinderToStringMap();
string getIdentifierFromSolverType(const RootFinderType solver_type);
RootFinderType getSolverTypeFromIdentifier(string solver_type);
rfs_map getRootFinderToNameMap();
string getNameFromSolverType(const RootFinderType solver_type);

void initPrintPointVector(std::vector<int> & print_points, int num_time_steps, int nprint, string print_points_file_name);
void printStateDataToVTKFile(string execName, std::ostringstream & vtkfilename, Opm::TwophaseState state, 
						 const UnstructuredGrid& grid, const RootFinderType solver_type, 
						 double comp_length_days, double time_step_days, int i /*, int plotInterval*/);
void printStateDataToVTKFile(string vtkfilename, Opm::TwophaseState state, const UnstructuredGrid& grid);						 
string replaceStrChar(string str, const string & replace, char ch);
string replaceDot(string str);
string replaceDot(double num);
void printIterationsFromVector(string execName, const Opm::TransportSolverTwophaseReorder & transport_solver, 
							   int i, int num_cells, const RootFinderType solver_type, 
							   const double comp_length, const double time_step, const double visc_ratio);
void parseArguments(int argc, char ** argv, double & muw, double & muo, 
					bool & verbose, double & time_step_days, double & comp_length_days,
					double & xsize, double & ysize, double & zsize, int & xdim, int & ydim, int & zdim, RootFinderType & solver_type, 
					bool & printIterations, int & nprint, string & print_points_file_name,
					string & perm_file_name, int & layer, double & xpos, double & ypos, double & perm, bool & is_inhom_perm, 
					double & srcVol, double & sinkVol, double & grav_x, double & grav_y, double & grav_z, bool & initBottomTop, bool & initLeftRight);
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
						
/*template <class Functor>
class PrintFunctor
{
	public:
	//template <class Functor> 
	static void printFunctorValues(const Functor& f, int n, const char * filename)
	{
		double dx = 1/(n-1);
		std::ofstream file; file.open(filename);
		file << "x \t y \n";
		for (int i = 0; i < n; i++)
		{
			file << dx*i << " \t" << f(dx*i) << "\n";
		}
		file.close();
	}
};*/
						 
}
