//===========================================================================
//
// File: RootFinders.hpp
//
// Created: Thu May  6 19:59:42 2010
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Jostein R Natvig    <jostein.r.natvig@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2010 Statoil ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_ROOTFINDERS_HEADER
#define OPM_ROOTFINDERS_HEADER

#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/StopWatch.hpp>

#include <algorithm>
#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstring>

namespace Opm
{

    struct ThrowOnError
    {
        static double handleBracketingFailure(const double x0, const double x1, const double f0, const double f1)
        {
            OPM_THROW(std::runtime_error, "Error in parameters, zero not bracketed: [a, b] = ["
                  << x0 << ", " << x1 << "]    f(a) = " << f0 << "   f(b) = " << f1);
            return -1e100; // Never reached.
        }
        static double handleTooManyIterations(const double x0, const double x1, const int maxiter)
        {
            OPM_THROW(std::runtime_error, "Maximum number of iterations exceeded: " << maxiter << "\n"
                  << "Current interval is [" << std::min(x0, x1) << ", "
                  << std::max(x0, x1) << "]");
            return -1e100; // Never reached.
        }
        static double handleTooManyIterationsNewton(const double x, const int maxiter, double funcval)
        {
			OPM_THROW(std::runtime_error, "Maximum number of iterations exceeded: " << maxiter << ". "
				<< "Current best root guess is " << x << " giving function value " << funcval << "\n");
            return -1e100; // Never reached.
		}
    };

    struct WarnAndContinueOnError
    {
        static double handleBracketingFailure(const double x0, const double x1, const double f0, const double f1)
        {
            OPM_MESSAGE("Error in parameters, zero not bracketed: [a, b] = ["
                    << x0 << ", " << x1 << "]    f(a) = " << f0 << "   f(b) = " << f1
                    << "");
            return std::fabs(f0) < std::fabs(f1) ? x0 : x1;
        }
        static double handleTooManyIterations(const double x0, const double x1, const int maxiter)
        {
            OPM_MESSAGE("Maximum number of iterations exceeded: " << maxiter
                    << ", current interval is [" << std::min(x0, x1) << ", "
                    << std::max(x0, x1) << "]");
            return 0.5*(x0 + x1);
        }
        static double handleTooManyIterationsNewton(const double x, const int maxiter, double funcval)
        {
			OPM_MESSAGE("Maximum number of iterations exceeded: " << maxiter << ". "
				<< "Current best root guess is " << x << " giving function value " << funcval << "\n");
            return x;
		}
    };

    struct ContinueOnError
    {
        static double handleBracketingFailure(const double x0, const double x1, const double f0, const double f1)
        {
            return std::fabs(f0) < std::fabs(f1) ? x0 : x1;
        }
        static double handleTooManyIterations(const double x0, const double x1, const int /*maxiter*/)
        {
            return 0.5*(x0 + x1);
        }
        static double handleTooManyIterationsNewton(const double x, const int /*maxiter*/, double /*funcval*/)
        {
            return x;
		}
    };
    
    /*class BaseRootFinder
    {
		template <class Functor>
		inline static virtual double solve(const Functor& f,
								   const double initial_guess,
								   const double a, const double b,
								   const int max_iter,
								   const double tolerance,
								   bool verbose,
								   int & iterations_used) = 0;
	};*/
    
    template <class ErrorPolicy = ThrowOnError>
    class Brent //: public virtual BaseRootFinder
    {
		public:
		
		template <class Functor>
		inline static double solve(const Functor& f,
								   const double initial_guess,
								   const double a, const double b,
								   const int max_iter,
								   const double tolerance,
								   bool verbose,
								   int & iterations_used)
        {
			double xa = a;    double xb = b;    double xc = xa; double xs = 1.0; double xd = 1.0;
			double fa = f(a); double fb = f(b); double fc = fa; double fs = f(initial_guess);
			
			if(fa*fs < 0)
			{
				xb = initial_guess;
				fb = fs;
			}
			else if(fb*fs < 0)
			{
				xa = initial_guess;
				fa = fs;
			}
			
			if(fa*fb >= 0)
			{
				if ( std::abs(fa) <= tolerance) return xa;
				else if ( std::abs(fb) <= tolerance) return xb;
				else return ErrorPolicy::handleBracketingFailure(a,b,fa,fb);
			}
			
			if(verbose)
				std::cout << "----------------------- Brent's Method iteration ---------------------------\n"
						<< "Initial guess: " << initial_guess << "\n"
						<< "Max iterations: " << max_iter << "\n"
						<< "Error tolerance: " << tolerance << "\n"
						<< "# iter.\tx\t\tf(x)\t\tf_x(x) \n";
						
			if (verbose) printf("%d\t%8.3e\t%8.2e\n",iterations_used,fb < fa ? xb : xa,f(fb < fa ? xb : xa));	
			
			bool flag = true;
			for (int i = 0; i < max_iter; i++)
			{
				++iterations_used;
				double fab = fa-fb; 
				if(fa != fc && fb != fc)
				{
					double fac = fa-fc; double fbc = fb-fc;
					xs = xa*fb*fc/(fab*fac) + xb*fa*fc/(-fab*fbc) + xc*fa*fb/(fac*fbc);
				}
				else
					xs = xb + fb*(xb-xa)/fab;
				
				bool sbbc = std::abs(xs-xb) >= std::abs(xb-xc)/2;
				bool sbcd = std::abs(xs-xb) >= std::abs(xc-xd)/2;
				bool bct = std::abs(xb-xc) < tolerance;
				bool cdt = std::abs(xc-xd) < tolerance;
				double tmp = (3*xa+xb)/4;
				if( !((xs < tmp && xs > xb) || (xs > tmp && xs < xb)) || (flag && sbbc) || (!flag && sbcd) || (flag && bct) || (!flag && cdt) )
				{
					xs = 0.5*(xa+xb); flag = true;
				}
				else
					flag = false;
					
				fs = f(xs);
				
				if ( std::abs(fs) <= tolerance ) 
				{
					if(verbose) std::cout << "----------------------- End Brent's Method iteration ---------------------------\n";
					return xs;
				}
				
				xd = xc; 
				xc = xb; fc = fb;
				
				if( fa*fs < 0)
				{
					xb = xs; fb = fs;
				}
				else
				{
					xa = xs; fa = fs;
				}
				if (std::abs(fa) < std::abs(fb))
				{
					double temp = xa;
					xa = xb;
					xb = temp;
					temp = fa;
					fa = fb;
					fb = temp;
				}
				
				if (verbose) printf("%d\t%8.3e\t%8.2e\n",iterations_used,fs < fb ? xs : xb,f(fs < fb ? xs : xb));
				
				if( std::abs(xa-xb) <= tolerance ) 
				{
					if(verbose) std::cout << "----------------------- End Brent's Method iteration ---------------------------\n";
					return fs < fb ? xs : xb;
				}
				
			}
			return ErrorPolicy::handleTooManyIterationsNewton(xs, max_iter, f(xs));
		}
		
		template <class Functor>
		inline static double solve(const Functor& f,
								   const double a, const double b,
								   const int max_iter,
								   const double tolerance,
								   int& iterations_used)
        {
			return solve(f,0.5*(a+b),a,b,max_iter,tolerance,false,iterations_used);
		}
		
		template <class Functor>
		inline static double solve(const Functor& f,
								   const double initial_guess,
								   const double a, const double b,
								   const int max_iter,
								   const double tolerance,
								   int& iterations_used)
        {
			return solve(f,initial_guess,a,b,max_iter,tolerance,false,iterations_used);
		}
	};

	template <class ErrorPolicy = ThrowOnError>
	class Ridder //: public virtual BaseRootFinder
	{
		public:
		
		template <class Functor>
		inline static double solve(const Functor& f,
								   const double initial_guess,
								   const double a, const double b,
								   const int max_iter,
								   const double tolerance,
								   bool verbose,
								   int& iterations_used)
		{
			double x = initial_guess;
			double xa = a;
			double xb = b;
			double fa = f(xa); double fb = f(xb);
			double fnew = f(x);
			
			if(fnew*fa < 0)
				xb = x;
			else if(fnew*fb < 0)
				xa = x;
			
			if(fa*fb >= 0)
			{
				if (std::abs(fa) <= tolerance) return xa;
				else if (std::abs(fb) <= tolerance) return xb;
				else return ErrorPolicy::handleBracketingFailure(a,b,fa,fb);
			}
			
			if(verbose)
				std::cout << "----------------------- Ridder's Method iteration ---------------------------\n"
						<< "Initial guess: " << initial_guess << "\n"
						<< "Max iterations: " << max_iter << "\n"
						<< "Error tolerance: " << tolerance << "\n"
						<< "# iter.\tx\t\tf(x)\t\tf_x(x) \n";
            
            if (verbose) printf("%d\t%8.3e\t%8.2e\n",iterations_used,x,f(x));
			
			x = -10; // A saturation value always outside bracket
			for(int i = 0; i < max_iter; i++)
			{
				++iterations_used;
				double xc = 0.5*(xa+xb);
				double fc = f(xc);
				double s = sqrt(fc*fc-fa*fb);
				if (s == 0)
				{
					if(verbose) std::cout << "----------------------- End Ridder's Method iteration ---------------------------\n";
					return x;
				}
				double xnew = xc + (xc-xa)*((fa >= fb ? 1.0 : -1.0)*fc/s);
				if( std::abs(xnew-x) <= tolerance ) 
				{
					if(verbose) std::cout << "----------------------- End Ridder's Method iteration ---------------------------\n";
					return x;
				}
				x = xnew;
				fnew = f(x);
				if( std::abs(fnew) <= tolerance )
				{
					if(verbose) std::cout << "----------------------- End Ridder's Method iteration ---------------------------\n";
					return x;
				}
				if(fc*fnew < 0)
				{
					xa = xc; fa = fc;
					xb = x; fb = fnew;
				}
				else if(fa*fnew < 0)
				{
					xb = x; fb = fnew;
				}
				else if(fb*fnew < 0)
				{
					xa = x; fa = fnew;
				}
				else
					std::cerr << "Bracket violation. This should never happen!\n";
					
				if(xb < xa)
				{
					double temp = xa;
					xa = xb;
					xb = temp;
					temp = fa;
					fa = fb;
					fb = temp;
				}
				
				if( std::abs(xb-xa) <= tolerance)
				{
					if(verbose) std::cout << "----------------------- End Ridder's Method iteration ---------------------------\n";
					return x;
				}
				
				if (verbose) printf("%d\t%8.3e\t%8.2e\t%8.2e \n",iterations_used,x,f(x),f.ds(x));
			}
			
			return ErrorPolicy::handleTooManyIterationsNewton(x, max_iter, f(x));
		}
		
		template <class Functor>
		inline static double solve(const Functor& f,
								   const double a, const double b,
								   const int max_iter,
								   const double tolerance,
								   bool verbose,
								   int& iterations_used)
		{
			return solve(f,0.5*(a+b),a,b,max_iter,tolerance,verbose,iterations_used);
		}
		
		template <class Functor>
		inline static double solve(const Functor& f,
								   const double a, const double b,
								   const int max_iter,
								   const double tolerance,
								   int& iterations_used)
		{
			return solve(f,0.5*(a+b),a,b,max_iter,tolerance,false,iterations_used);
		}
		
		
		template <class Functor>
		inline static double solve(const Functor& f,
								   const double initial_guess,
								   const double a, const double b,
								   const int max_iter,
								   const double tolerance,
								   int& iterations_used)
		{
			return solve(f,initial_guess,a,b,max_iter,tolerance,false,iterations_used);
		}
	};
	
	template <class ErrorPolicy = ThrowOnError>
    class RegulaFalsiTrustRegion
    {
    public:

        /// Implements a modified regula falsi method as described in
        /// "Improved algorithms of Illinois-type for the numerical
        /// solution of nonlinear equations"
        /// by J. A. Ford.
        /// Current variant is the 'Pegasus' method.
        template <class Functor>
        inline static double solve(const Functor& f,
                                   const double a,
                                   const double b,
                                   const double visc_ratio,
                                   const int max_iter,
                                   const double tolerance,
                                   const bool verbose,
                                   int& iterations_used)
        {
            using namespace std;
            const double macheps = numeric_limits<double>::epsilon();
            const double eps = tolerance + macheps*max(max(fabs(a), fabs(b)), 1.0);

            double x0 = a;
            double x1 = b;
            double f0 = f(x0);
            const double epsF = tolerance + macheps*max(fabs(f0), 1.0);
            if (fabs(f0) < epsF) {
                return x0;
            }
            double f1 = f(x1);
            if (fabs(f1) < epsF) {
                return x1;
            }
            if (f0*f1 > 0.0) {
                return ErrorPolicy::handleBracketingFailure(a, b, f0, f1);
            }
            
            return loop(f,x0,x1,f0,f1,visc_ratio,max_iter,eps,epsF,iterations_used);
        }
        
        /// Implements a modified regula falsi method as described in
        /// "Improved algorithms of Illinois-type for the numerical
        /// solution of nonlinear equations"
        /// by J. A. Ford.
        /// Current variant is the 'Pegasus' method.
        /// This version takes an extra parameter for the initial guess.
        template <class Functor>
        inline static double solve(const Functor& f,
                                   const double initial_guess,
                                   const double a,
                                   const double b,
                                   const double visc_ratio,
                                   const int max_iter,
                                   const double tolerance,
                                   const bool verbose,
                                   int& iterations_used)
        {
            using namespace std;
            const double macheps = numeric_limits<double>::epsilon();
            const double eps = tolerance + macheps*max(max(fabs(a), fabs(b)), 1.0);

            double f_initial = f(initial_guess);
            const double epsF = tolerance + macheps*max(fabs(f_initial), 1.0);
            if (fabs(f_initial) < epsF) {
                return initial_guess;
            }
            double x0 = a;
            double x1 = b;
            double f0 = f_initial;
            double f1 = f_initial;
            if (x0 != initial_guess) {
                f0 = f(x0);
                if (fabs(f0) < epsF) {
                    return x0;
                }
            }
            if (x1 != initial_guess) {
                f1 = f(x1);
                if (fabs(f1) < epsF) {
                    return x1;
                }
            }
            if (f0*f_initial < 0.0) {
                x1 = initial_guess;
                f1 = f_initial;
            } else {
                x0 = initial_guess;
                f0 = f_initial;
            }
            if (f0*f1 > 0.0) {
                return ErrorPolicy::handleBracketingFailure(a, b, f0, f1);
            }
            
            return loop(f,x0,x1,f0,f1,visc_ratio,max_iter,eps,epsF,iterations_used);
        }

		template <class Functor>
        inline static double loop(const Functor& f,
                                   double x0,
                                   double x1,
                                   double f0,
                                   double f1,
                                   const double visc_ratio,
                                   const int max_iter,
                                   const double eps,
                                   const double epsF,
                                   int& iterations_used)
        {
			iterations_used = 0;
            // In every iteraton, x1 is the last point computed,
            // and x0 is the last point computed that makes it a bracket.
            while (fabs(x1 - x0) >= 1e-9*eps) {
				double xold = x1;
                double xnew = regulaFalsiStep(x0, x1, f0, f1);
                
                // Trust region
                xnew = std::max(std::min(xnew,1.0),0.0);
				if(dfw2(xnew,visc_ratio)*dfw2(xold,visc_ratio) < 0.0)
					xnew = (xnew+xold)/2.0;
                
                double fnew = f(xnew);
                ++iterations_used;
                if (iterations_used > max_iter) {
                    return ErrorPolicy::handleTooManyIterations(x0, x1, max_iter);
                }
                if (fabs(fnew) < epsF) {
                    return xnew;
                }
                // Now we must check which point we must replace.
                if ((fnew > 0.0) == (f0 > 0.0)) {
                    // We must replace x0.
                    x0 = x1;
                    f0 = f1;
                } else {
                    // We must replace x1, this is the case where
                    // the modification to regula falsi kicks in,
                    // by modifying f0.
                    // The 'Pegasus' method
                    const double gamma = f1/(f1 + fnew);
                    f0 *= gamma;
                }
                x1 = xnew;
                f1 = fnew;
            }
            return 0.5*(x0 + x1);
            
					
					
					
		}

    private:
        inline static double regulaFalsiStep(const double a,
                                             const double b,
                                             const double fa,
                                             const double fb)
        {
            assert(fa*fb < 0.0);
            return (b*fa - a*fb)/(fa - fb);
        }
		inline static double dfw2(double x, double M)
		{
			double x2 = pow(x,2.0);
			double xm2 = pow(x-1,2.0);
			return 2*M*( M*(2*x+1)*xm2 + x2*(2*x-3) )/pow(M*xm2 + x2, 3.0);
		}

    };
	
	// Newton Raphson Trust region solver
    template <class ErrorPolicy = ThrowOnError>
    class NewtonRaphsonTrustRegion
    {
		public:
		
		// Trust region solver
		// Uses supplied first and second derivatives for trust region scheme
		template <class Functor>
		inline static double solve(const Functor& f,
								   const double initial_guess,
								   const int max_iter,
								   const double tolerance,
								   const bool verbose,
								   int& iterations_used)
        {
            double x = initial_guess;
            double xNew = initial_guess;
            double xCorr = initial_guess;
            
            if(verbose)
            std::cout << "----------------------- Newton Trust Region iteration ---------------------------\n"
					  << "Initial guess: " << initial_guess << "\n"
                      << "Max iterations: " << max_iter << "\n"
                      << "Error tolerance: " << tolerance << "\n"
                      << "# iter.\tx\t\txCorr\t\tf(x)\t\tf_x(x) \n";
            
            while (std::abs(f(x)) > tolerance)
            {
				++iterations_used;
				if (iterations_used > max_iter)
                    return ErrorPolicy::handleTooManyIterationsNewton(x, max_iter, f(x));
                
				xNew = x - f(x)/f.ds(x);
				
				xCorr = std::max(std::min(xNew,1.0),0.0);
				if(f.ds2(xCorr)*f.ds2(x) < 0.0)
					xCorr = (xCorr+x)/2.0;
				
				if (verbose) printf("%d\t%8.2e\t%8.3e\t%8.2e\t%8.2e \n",iterations_used,xNew,xCorr,f(x),f.ds(x));
				
				x = xCorr;
			}
			if(verbose) std::cout << "---------------------- End Newton Trust Region iteration ------------------------\n";
				
            return x;
		}
		
		// Darcy flow trust region solver
		// Uses supplied viscosity ratio for trust region scheme
		// Reduced number of function calls by approximating dfw2(xNew)
		template <class Functor>
		inline static double solveApprox(const Functor& f,
								   const double initial_guess,
								   const double visc_ratio,
								   const int max_iter,
								   const double tolerance,
								   const bool verbose,
								   int& iterations_used)
        {
            double x = initial_guess;
            double xNew = initial_guess;
            
            if(verbose)
            std::cout << "----------------------- Newton Trust Region iteration ---------------------------\n"
					  << "Initial guess: " << initial_guess << "\n"
                      << "Max iterations: " << max_iter << "\n"
                      << "Error tolerance: " << tolerance << "\n"
                      << "# iter.\tx\t\tf(x)\t\tf_x(x) \n";
            
            double dfw2x = dfw2(x,visc_ratio);
            double dfw2xNew = 0;
            while (std::abs(f(x)) > tolerance)
            {
				++iterations_used;
				if (iterations_used > max_iter)
                    return ErrorPolicy::handleTooManyIterationsNewton(x, max_iter, f(x));
                
				// Scalar Newton's method
				xNew = x - f(x)/f.ds(x);
				
				// Restrict to interval [0,1]
				xNew = std::max(std::min(xNew,1.0),0.0);
				
				dfw2xNew = dfw2(xNew,visc_ratio);
				if(dfw2xNew*dfw2x < 0.0)
				{
					xNew = (xNew+x)*0.5;
					dfw2x = (dfw2xNew+dfw2x)*0.5;
				}
				else
					dfw2x = dfw2xNew;
				
				if (verbose)
					printf("%d\t%8.2e\t%8.2e\t%8.2e \n",iterations_used,xNew,f(x),f.ds(x));
				
				x = xNew;
			}
			if(verbose)
				std::cout << "---------------------- End Newton iteration ------------------------\n";
			
            return x;
		}
		
		// Darcy flow trust region solver
		// Uses supplied viscosity ratio for trust region scheme
		template <class Functor>
		inline static double solve(const Functor& f,
								   const double initial_guess,
								   const double visc_ratio,
								   const int max_iter,
								   const double tolerance,
								   const bool verbose,
								   int& iterations_used)
        {
            double x = initial_guess;
            double xNew = initial_guess;
            
            /*double thePoint = 0.4;
            std::string line;
            std::ifstream infile("pointdata.txt");
            if(infile.is_open())
            {
				if( getline(infile,line) )
					thePoint = atof(line.c_str());
				infile.close();
			}
			double xval;
            double xmin = 0; double xmax = 1;
            int n_points = 150;
            std::ofstream file;
            file.open("cell_residual.txt");
            for ( int i = 0; i <= n_points; i++)
            {
				xval = (xmax-xmin)/n_points*i;
				file << xval << "\t" << f(xval) << "\t" << f.ds(xval) << "\t" << fw(xval,visc_ratio) << "\t" << dfw2(xval,visc_ratio) << "\t" << f(thePoint) + f.ds(thePoint)*(xval-thePoint) << "\n";
			}
			file.close();*/
            
            if(verbose)
            std::cout << "----------------------- Newton Trust Region iteration ---------------------------\n"
					  << "Initial guess: " << initial_guess << "\n"
                      << "Max iterations: " << max_iter << "\n"
                      << "Error tolerance: " << tolerance << "\n"
                      << "# iter.\tx\t\tf(x)\t\tf_x(x) \n";
            
            while (std::abs(f(x)) > tolerance)
            {
				++iterations_used;
				if (iterations_used > max_iter)
                    return ErrorPolicy::handleTooManyIterationsNewton(x, max_iter, f(x));
                
				// Scalar Newton's method
				xNew = x - f(x)/f.ds(x);
				
				// Restrict to interval [0,1]
				xNew = std::max(std::min(xNew,1.0),0.0);
				
				// Trust Region
				if(dfw2(xNew,visc_ratio)*dfw2(x,visc_ratio) < 0.0)
					xNew = (xNew+x)/2.0;
				
				if (verbose)
					printf("%d\t%8.2e\t%8.2e\t%8.2e \n",iterations_used,xNew,f(x),f.ds(x));
				
				x = xNew;
			}
			if(verbose)
				std::cout << "---------------------- End Newton iteration ------------------------\n";
			
            return x;
		}
		
		inline static double fw(double x, double M)
		{
			double x2 = pow(x,2.0);
			return x2/(x2+M*pow((1-x),2.0));
		}
		inline static double dfw(double x, double M)
		{
			return 2*M*x*(1-x)/pow( M*pow(x-1,2.0) + pow(x,2.0), 2.0);
		}
		inline static double dfw2(double x, double M)
		{
			double x2 = pow(x,2.0);
			double xm2 = pow(x-1,2.0);
			return 2*M*( M*(2*x+1)*xm2 + x2*(2*x-3) )/pow(M*xm2 + x2, 3.0);
		}
		
		template <class Functor>
		inline static double solve(const Functor& f,
								   const double a,
                                   const double b,
								   const int max_iter,
								   const double tolerance,
								   int& iterations_used)
		{
			double initial_guess = (a+b)/2.0;
			return solve(f,initial_guess,max_iter,tolerance,false,iterations_used);
		}
		
		template <class Functor>
		inline static double solve(const Functor& f,
								   const double initial_guess,
								   const double a,
                                   const double b,
								   const int max_iter,
								   const double tolerance,
								   int& iterations_used)
		{
			return solve(f,initial_guess,max_iter,tolerance,false,iterations_used);
		}
	};

	// Newton Raphson solver
    template <class ErrorPolicy = ThrowOnError>
    class NewtonRaphson
    {
		public:
		// Newton Raphson solver
		template <class Functor>
		inline static double solve(const Functor& f,
								   const double initial_guess,
								   const double a, const double b,
								   const int max_iter,
								   const double tolerance,
								   bool verbose,
								   int& iterations_used)
        {
            double x = initial_guess;
            double xNew = initial_guess;
            
            if(verbose)
            std::cout << "----------------------- Newton iteration ---------------------------\n"
					  << "Initial guess: " << initial_guess << "\n"
                      << "Max iterations: " << max_iter << "\n"
                      << "Error tolerance: " << tolerance << "\n"
                      << "# iter.\tx\t\tf(x)\t\tf_x(x) \n";
            
            while (std::abs(f(x)) > tolerance)
            {
				++iterations_used;
				if (iterations_used > max_iter)
                    return ErrorPolicy::handleTooManyIterationsNewton(x, max_iter, f(x));
                
				xNew = x - f(x)/f.ds(x);
				if (verbose) printf("%d\t%8.3e\t%8.2e\t%8.2e \n",iterations_used,xNew,f(x),f.ds(x));
				
				x = xNew;
			}
			if(verbose) std::cout << "---------------------- End Newton iteration ------------------------\n";
				
            return x;
		}
		
		template <class Functor>
		inline static double solve(const Functor& f,
								   const double a,
                                   const double b,
								   const int max_iter,
								   const double tolerance,
								   int& iterations_used)
		{
			double initial_guess = (a+b)/2.0;
			return solve(f,initial_guess,a,b,max_iter,tolerance,false,iterations_used);
		}
		
		template <class Functor>
		inline static double solve(const Functor& f,
								   const double initial_guess,
								   const double a,
                                   const double b,
								   const int max_iter,
								   const double tolerance,
								   int& iterations_used)
		{
			return solve(f,initial_guess,a,b,max_iter,tolerance,false,iterations_used);
		}
		
		template <class Functor>
		inline static double solve(const Functor& f,
								   const double a,
                                   const double b,
								   const int max_iter,
								   const double tolerance,
								   bool verbose,
								   int& iterations_used)
		{
			double initial_guess = (a+b)/2.0;
			return solve(f,initial_guess,a,b,max_iter,tolerance,verbose,iterations_used);
		}
	};


    template <class ErrorPolicy = ThrowOnError>
    class RegulaFalsi
    {
    public:


        /// Implements a modified regula falsi method as described in
        /// "Improved algorithms of Illinois-type for the numerical
        /// solution of nonlinear equations"
        /// by J. A. Ford.
        /// Current variant is the 'Pegasus' method.
        template <class Functor>
        inline static double solve(const Functor& f,
                                   const double a,
                                   const double b,
                                   const int max_iter,
                                   const double tolerance,
                                   int& iterations_used)
        {
            using namespace std;
            const double macheps = numeric_limits<double>::epsilon();
            const double eps = tolerance + macheps*max(max(fabs(a), fabs(b)), 1.0);

            double x0 = a;
            double x1 = b;
            double f0 = f(x0);
            const double epsF = tolerance + macheps*max(fabs(f0), 1.0);
            if (fabs(f0) < epsF) {
                return x0;
            }
            double f1 = f(x1);
            if (fabs(f1) < epsF) {
                return x1;
            }
            if (f0*f1 > 0.0) {
                return ErrorPolicy::handleBracketingFailure(a, b, f0, f1);
            }
            iterations_used = 0;
            // In every iteraton, x1 is the last point computed,
            // and x0 is the last point computed that makes it a bracket.
            while (fabs(x1 - x0) >= 1e-9*eps) {
                double xnew = regulaFalsiStep(x0, x1, f0, f1);
                double fnew = f(xnew);
// 		cout << "xnew = " << xnew << "    fnew = " << fnew << endl;
                ++iterations_used;
                if (iterations_used > max_iter) {
                    return ErrorPolicy::handleTooManyIterations(x0, x1, max_iter);
                }
                if (fabs(fnew) < epsF) {
                    return xnew;
                }
                // Now we must check which point we must replace.
                if ((fnew > 0.0) == (f0 > 0.0)) {
                    // We must replace x0.
                    x0 = x1;
                    f0 = f1;
                } else {
                    // We must replace x1, this is the case where
                    // the modification to regula falsi kicks in,
                    // by modifying f0.
                    // 1. The classic Illinois method
//                  const double gamma = 0.5;
                    // @afr: The next two methods do not work??!!?
                    // 2. The method called 'Method 3' in the paper.
//                  const double phi0 = f1/f0;
//                  const double phi1 = fnew/f1;
//                  const double gamma = 1.0 - phi1/(1.0 - phi0);
                    // 3. The method called 'Method 4' in the paper.
//                  const double phi0 = f1/f0;
//                  const double phi1 = fnew/f1;
//                  const double gamma = 1.0 - phi0 - phi1;
//                  cout << "phi0 = " << phi0 <<" phi1 = " << phi1 <<
//                  " gamma = " << gamma << endl;
                    // 4. The 'Pegasus' method
                    const double gamma = f1/(f1 + fnew);
                    f0 *= gamma;
                }
                x1 = xnew;
                f1 = fnew;
            }
            return 0.5*(x0 + x1);
        }


        /// Implements a modified regula falsi method as described in
        /// "Improved algorithms of Illinois-type for the numerical
        /// solution of nonlinear equations"
        /// by J. A. Ford.
        /// Current variant is the 'Pegasus' method.
        /// This version takes an extra parameter for the initial guess.
        template <class Functor>
        inline static double solve(const Functor& f,
                                   const double initial_guess,
                                   const double a,
                                   const double b,
                                   const int max_iter,
                                   const double tolerance,
                                   int& iterations_used)
        {
            using namespace std;
            const double macheps = numeric_limits<double>::epsilon();
            const double eps = tolerance + macheps*max(max(fabs(a), fabs(b)), 1.0);

            double f_initial = f(initial_guess);
            const double epsF = tolerance + macheps*max(fabs(f_initial), 1.0);
            if (fabs(f_initial) < epsF) {
                return initial_guess;
            }
            double x0 = a;
            double x1 = b;
            double f0 = f_initial;
            double f1 = f_initial;
            if (x0 != initial_guess) {
                f0 = f(x0);
                if (fabs(f0) < epsF) {
                    return x0;
                }
            }
            if (x1 != initial_guess) {
                f1 = f(x1);
                if (fabs(f1) < epsF) {
                    return x1;
                }
            }
            if (f0*f_initial < 0.0) {
                x1 = initial_guess;
                f1 = f_initial;
            } else {
                x0 = initial_guess;
                f0 = f_initial;
            }
            if (f0*f1 > 0.0) {
                return ErrorPolicy::handleBracketingFailure(a, b, f0, f1);
            }
            iterations_used = 0;
            // In every iteraton, x1 is the last point computed,
            // and x0 is the last point computed that makes it a bracket.
            while (fabs(x1 - x0) >= 1e-9*eps) {
                double xnew = regulaFalsiStep(x0, x1, f0, f1);
                double fnew = f(xnew);
// 		cout << "xnew = " << xnew << "    fnew = " << fnew << endl;
                ++iterations_used;
                if (iterations_used > max_iter) {
                    return ErrorPolicy::handleTooManyIterations(x0, x1, max_iter);
                }
                if (fabs(fnew) < epsF) {
                    return xnew;
                }
                // Now we must check which point we must replace.
                if ((fnew > 0.0) == (f0 > 0.0)) {
                    // We must replace x0.
                    x0 = x1;
                    f0 = f1;
                } else {
                    // We must replace x1, this is the case where
                    // the modification to regula falsi kicks in,
                    // by modifying f0.
                    // 1. The classic Illinois method
//      const double gamma = 0.5;
                    // @afr: The next two methods do not work??!!?
                    // 2. The method called 'Method 3' in the paper.
// 		    const double phi0 = f1/f0;
// 		    const double phi1 = fnew/f1;
// 		    const double gamma = 1.0 - phi1/(1.0 - phi0);
                    // 3. The method called 'Method 4' in the paper.
// 		    const double phi0 = f1/f0;
// 		    const double phi1 = fnew/f1;
// 		    const double gamma = 1.0 - phi0 - phi1;
// 		    cout << "phi0 = " << phi0 <<" phi1 = " << phi1 <<
// 		    " gamma = " << gamma << endl;
                    // 4. The 'Pegasus' method
                    const double gamma = f1/(f1 + fnew);
                    f0 *= gamma;
                }
                x1 = xnew;
                f1 = fnew;
            }
            return 0.5*(x0 + x1);
        }


    private:
        inline static double regulaFalsiStep(const double a,
                                             const double b,
                                             const double fa,
                                             const double fb)
        {
            assert(fa*fb < 0.0);
            return (b*fa - a*fb)/(fa - fb);
        }


    };


    /// Attempts to find an interval bracketing a zero by successive
    /// enlargement of search interval.
    template <class Functor>
    inline void bracketZero(const Functor& f,
                            const double x0,
                            const double dx,
                            double& a,
                            double& b)
    {
        const int max_iters = 100;
        double f0 = f(x0);
        double cur_dx = dx;
        int i = 0;
        for (; i < max_iters; ++i) {
            double x = x0 + cur_dx;
            double f_new = f(x);
            if (f0*f_new <= 0.0) {
                break;
            }
            cur_dx = -2.0*cur_dx;
        }
        if (i == max_iters) {
            OPM_THROW(std::runtime_error, "Could not bracket zero in " << max_iters << "iterations.");
        }
        if (cur_dx < 0.0) {
            a = x0 + cur_dx;
            b = i < 2 ? x0 : x0 + 0.25*cur_dx;
        } else {
            a = i < 2 ? x0 : x0 + 0.25*cur_dx;
            b = x0 + cur_dx;
        }
    }





} // namespace Opm




#endif // OPM_ROOTFINDERS_HEADER
