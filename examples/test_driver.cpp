#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <iostream>
#include <opm/core/utility/RootFinders.hpp>
#include <math.h>

using namespace Opm;

struct Function1
{
	double operator()(double x) const
	{
		return exp(x)-exp(1);
	}
	double ds(double x) const
	{
		return exp(x);
	}
};

struct Function2
{
	static constexpr double M = 10.0;
	double operator()(double x) const
	{
		double x2 = pow(x,2.0);
		return x2/(x2+M*pow((1-x),2.0)) - 0.5;
	}
	double ds(double x) const
	{
		return 2*M*x*(1-x)/pow( M*pow(x-1,2.0) + pow(x,2.0), 2.0);
	}
	double ds2(double x) const
	{
		double x2 = pow(x,2.0);
		double xm2 = pow(x-1,2.0);
		return 2*M*( M*(2*x+1)*xm2 + x2*(2*x-3) )/pow(M*xm2 + x2, 3.0);
	}
};

typedef Function2 Function;
typedef NewtonRaphson<ThrowOnError> RootFinder;

int main(int argc, char ** argv)
{
	// -i: Initial guess. Default 0
	// -t: Tolerance. Default 1e-6
	// -m: Max iterations. Default 20
	// -v: Verbose. Default false
	double tolerance = 1e-6;
	double init_guess = 0.1;
	int max_iter = 20;
	bool verbose = false;
	
	for(int i = 1; i < argc; i = i+2)
	{
		if(std::string(argv[i]) == "-i")
			init_guess = atof(argv[i+1]);
		else if(std::string(argv[i]) == "-t")
			tolerance = atof(argv[i+1]);
		else if(std::string(argv[i]) ==	"-m")
			max_iter = atoi(argv[i+1]);
		else if(std::string(argv[i]) == "-v")
		{
			verbose = true;
			--i;
		}
		else
			std::cerr << "Invalid argument passed to" << argv[0] << "\n";
	}
	
	Function f;
	int iters_used = 0.0;

	double root = RootFinder::solve(f, init_guess, max_iter, tolerance, verbose, iters_used);
	std::cout << "A root candidate x = " << root << " was found after " << iters_used << " iterations.\n";
	
	return 0;
}
