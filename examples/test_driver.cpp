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
		return exp(x)-exp(1)+0.5;
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

typedef struct ThrowOnError TheErrorPolicy;

/*typedef boost::variant<RegulaFalsi<TheErrorPolicy>*, 
					   TrustRegion<TheErrorPolicy>*, 
					   Ridder<TheErrorPolicy>*,
					   Brent<TheErrorPolicy>*, 
					   NewtonRaphson<TheErrorPolicy>*> solverVariant;

template <class ErrorPolicy = ThrowOnError>
class PickSolver
{
	public:
	static solverVariant selectSolver(const char solver_type)
	{
		if(solver_type == 'r')
		{
			RegulaFalsi<ErrorPolicy> * solver = new RegulaFalsi<ErrorPolicy>();
			return solver;
		}
		else if(solver_type == 't')
		{
			TrustRegion<ErrorPolicy> * solver = new TrustRegion<ErrorPolicy>();
			return solver;
		}
		else if(solver_type == 'i')
		{
			Ridder<ErrorPolicy> * solver = new Ridder<ErrorPolicy>();
			return solver;
		}
		else if(solver_type == 'b')
		{
			Brent<ErrorPolicy> * solver = new Brent<ErrorPolicy>();
			return solver;
		}
		else if(solver_type == 'n')
		{
			NewtonRaphson<ErrorPolicy> * solver = new NewtonRaphson<ErrorPolicy>();
			return solver;
		}
		else
		{
			RegulaFalsi<ErrorPolicy> * solver = new RegulaFalsi<ErrorPolicy>();
			return solver;
		}
	}
		
};*/

/*BaseRootFinder * selectSolver(char solver_type)
{
	BaseRootFinder * solver;
	if(solver_type == 'b')
	{
		solver = new Brent<ThrowOnError>(); 
		return solver;
	}
	else if(solver_type == 'i')
	{
		solver = new Ridder<ThrowOnError>();
		return solver;
	}
}*/

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
	char solver_type = 'r';
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
		else if(std::string(argv[i]) == "-s")
		{
			solver_type = *argv[i+1];
		}
		else
			std::cerr << "Invalid argument passed to " << argv[0] << "\n";
	}
	
	Function f;
	int iters_used = 0.0;
	double root = -1e100;
	
	//solverVariant solver = PickSolver<ThrowOnError>::selectSolver(solver_type); 
	//BaseRootFinder * solver = selectSolver(solver_type);
	//root = solver->solve(f, init_guess, 0.0, 1.0, max_iter, tolerance, verbose, iters_used);
	
	if(solver_type == 'b')
	{
		std::cout << "Brent: \n";
		root = Brent<TheErrorPolicy>::solve(f, init_guess, 0.0, 1.0, max_iter, tolerance, verbose, iters_used);
	}
	else if(solver_type == 't')
	{
		std::cout << "Trust Region: \n";
		root = TrustRegion<TheErrorPolicy>::solve(f, init_guess, max_iter, tolerance, verbose, iters_used);
	}
	else if(solver_type == 'r')
	{
		std::cout << "Regula Falsi: \n";
		root = RegulaFalsi<TheErrorPolicy>::solve(f, init_guess, 0.0, 1.0, max_iter, tolerance, iters_used);
	}
	else if(solver_type == 'n')
	{
		std::cout << "Newton Raphson: \n";
		root = NewtonRaphson<TheErrorPolicy>::solve(f, init_guess, 0.0, 1.0, max_iter, tolerance, verbose, iters_used);
	}
	else if(solver_type == 'i')
	{
		std::cout << "Ridder's method: \n";
		root = Ridder<TheErrorPolicy>::solve(f, init_guess, 0.0, 1.0, max_iter, tolerance, verbose, iters_used);
	}
	else
	{
		std::cout << "Rergula Falsi(default): \n";
		root = RegulaFalsi<TheErrorPolicy>::solve(f, init_guess, 0.0, 1.0, max_iter, tolerance, iters_used);
	}
	
	
	std::cout << "A root x = " << root << " giving f(x) = " << f(root) << " was found after " << iters_used << " iterations.\n";
	
	return 0;
}
