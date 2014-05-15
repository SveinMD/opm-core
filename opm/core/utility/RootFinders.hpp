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
#include <opm/core/utility/CaseUtilities.hpp>

#include <algorithm>
#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstring>
#include <utility>

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
    
    template <class Functor>
class PrintFunctor
{
	public:
	//template <class Functor> 
	static void printFunctorValues(const Functor& f, int n, const char * filename)
	{
		double dx = 1.0/(n-1.0);
		std::ofstream file; file.open(filename);
		file << "x, \t y \n";
		for (int i = 0; i < n; i++)
		{
			file << (double)dx*i << ", \t" << (double)f(dx*i) << "\n";
		}
		file.close();
	}
};
    
    void addPointToVector(double x, double y, std::vector<std::pair<double,double>> & vec)
    {
		std::pair<double,double> thepair(x,y);
		vec.push_back(thepair);
	}
    double sign(double a,double b) {
		return ( (b) >= 0.0 ? fabs(a) : -fabs(a) );
	}
	double newt(double guess, double s, int iter)
	{
		//std::cout << guess << "\n";
		double prev_guess = guess;
		guess = 2*guess/(guess*guess+s);
		//std::cout << guess << "\n";
		int i = 0;
		while(i++ < iter && fabs(prev_guess - guess) > 10e-4)
		{
			prev_guess = guess;
			guess = 0.5*guess*(3-s*guess*guess);
			//std::cout << guess << "\n";
		}
		return guess;
	}
    
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
								   const double x1, const double x2,
								   const int max_iter,
								   const double tol,
								   bool verbose,
								   int & iterations_used, 
								   bool isTestRun, 
								   std::vector<std::pair<double,double>> & solution_path)
        {
			double eps = 3.0e-8;
			int iter;
			double a=x1,b=x2,c=x2,d=a,e=a,min1,min2;
			double fa=f(a),fb=f(b),fc=f(initial_guess),p,q,r,s,tol1,xm;
			
			if(fa*fc < 0)
			{
				b = initial_guess; fb = fc;
			}
			else if(fb*fc < 0)
			{
				a = initial_guess; fa = fc;
			}
			if(fa*fb >= 0)
			{
				if ( std::abs(fa) <= tol) return a;
				else if ( std::abs(fb) <= tol) return b;
				else return ErrorPolicy::handleBracketingFailure(a,b,fa,fb);
			}
			
			fc = fb;
			
			if(isTestRun)	
				addPointToVector(initial_guess,f(initial_guess),solution_path); // REMOVE FOR SPEED TESTS
			
			for(iter = 1; iter <= max_iter; iter++)
			{
				++iterations_used;
				if( (fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) 
				{
					c = a; fc = fa;
					e = d = b-a;
				}
				if( fabs(fc) < fabs(fb) ) 
				{
					a=b; 	b=c; 	c=a;
					fa=fb; 	fb=fc; 	fc=fa;
				}
				tol1=2.0*eps*fabs(b)+0.5*tol;
				//std::cout << tol1 << std::endl;
				xm = 0.5*(c-b);
				//if(fabs(xm) <= tol1 || fb == 0.0) return b;
				if(fabs(xm) <= tol1 || fabs(fb) <= tol1) return b;
				if(fabs(e) >= tol1 && fabs(fa) > fabs(fb)) 
				{
					s = fb/fa;
					if(a == c)
					{
						p=2.0*xm*s; q = 1.0-s;
					}
					else
					{
						q = fa/fc;
						r = fb/fc;
						p = s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
						q=(q-1.0)*(r-1.0)*(s-1.0);
					}
					if(p > 0.0) q = -q;
					p = fabs(p);
					min1=3.0*xm*q-fabs(tol1*q);
					min2=fabs(e*q);
					if(2.0*p < (min1 < min2 ? min1 : min2)) 
					{
						e=d; d = p/q;
					}
					else
					{
						d=xm; e=d;
					}
				}
				else
				{
					d=xm; e=d;
				}
				a=b; fa=fb;
				if(fabs(d) > tol1)
					b += d;
				else
					b += sign(tol1,xm);
				fb = f(b);
				if(isTestRun)
					addPointToVector(b,fb,solution_path); // REMOVE FOR SPEED TESTS
				
				/*if(iter > 30)
				{
					std::string file_name = "functor_values.data";
					PrintFunctor<Functor>::printFunctorValues(f,100,file_name.c_str());
					return ErrorPolicy::handleTooManyIterationsNewton(b, max_iter, f(b));
				}*/
			}
			
			return ErrorPolicy::handleTooManyIterationsNewton(a,b,max_iter);
		}
		
		template <class Functor>
		inline static double solve(const Functor& f,
								   const double initial_guess,
								   const double a, const double b,
								   const int max_iter,
								   const double tolerance,
								   bool verbose,
								   int & iterations_used,
								   bool dummy)
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
								   double x1, double x2,
								   const int max_iter,
								   const double tolerance,
								   bool verbose,
								   int& iterations_used//,
								   //double dummy_variable
								   , bool isTestRun, std::vector<std::pair<double,double>> & solution_path)
		{
			double invalid_ans = -1;
			int j;
			double ans,fh,fl,fm,fnew,s,xh,xl,xm,xnew=initial_guess,f_init=f(initial_guess);
			
			if(isTestRun)
				addPointToVector(initial_guess,f(initial_guess),solution_path); // REMOVE FOR SPEED TESTS
			
			fh = f(x2);
			if(fh*f_init < 0)
			{
				fl = f_init;
				x1 = initial_guess;
			}
			else
			{
				fl = f(x1);
				if(fl*f_init < 0)
				{
					fh = f_init;
					x2 = initial_guess;
				}
			}
			
			if( fl*fh < 0.0 ) {
				xl = x1; xh = x2;
				ans = invalid_ans;
				for(j = 1; j <= max_iter; j++) {
					++iterations_used;
					xm = 0.5*(xl+xh);
					fm = f(xm);
					s = sqrt(fm*fm-fl*fh); //newt(0.5,fm*fm-fl*fh,3); // Note: The iterative method returns the reciprocal root!
					if(s == 0.0) return ans; //if(s == 0.0) return ans; // if(fabs(s) <= tolerance) return ans; // 
					xnew = xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s); // *s); //
					if(isTestRun)
						addPointToVector(xnew,f(xnew),solution_path); // REMOVE FOR SPEED TESTS
					if (fabs(xnew-ans) <= tolerance) return xnew;
					ans = xnew;
					fnew = f(ans);
					if(fabs(fnew) <= tolerance) return ans; // if(fnew == 0.0) return ans;
					if(sign(fm,fnew) != fm) {
						xl = xm;  fl = fm;
						xh = ans; fh = fnew;
					}
					else if(sign(fl,fnew) != fl) {
						xh = ans; fh = fnew;
					}
					else if(sign(fh,fnew) != fh) {
						xl = ans; fl = fnew;
					}
					else {
						ErrorPolicy::handleBracketingFailure(xl,xh,fl,fh);
					}
					if(fabs(xh-xl) <= tolerance)
						return ans;
				}
				return ErrorPolicy::handleTooManyIterationsNewton(xnew, max_iter, f(xnew));
			}
			else {
				if(isTestRun)
					addPointToVector(x1,f(x1),solution_path); // REMOVE FOR SPEED TESTS
				if(fabs(fl) <= tolerance) return x1; // if(fl == 0.0) return x1;
				if(isTestRun)
					addPointToVector(x2,f(x2),solution_path); // REMOVE FOR SPEED TESTS
				if(fabs(fh) <= tolerance) return x2; // if(fh == 0.0) return x2;
				return ErrorPolicy::handleBracketingFailure(x1,x2,fl,fh);
			}
			return 0.0; // Unreachable
		}
		
		
		template <class Functor>
		inline static double solve(const Functor& f,
								   const double initial_guess,
								   const double a, const double b,
								   const int max_iter,
								   const double tolerance,
								   bool verbose,
								   int& iterations_used,
								   bool dummy
								   )
		{
			double s = 0.0;
			double x = initial_guess;
			double xa = a;
			double xb = b;
			double fa = f(xa); double fb = f(xb);
			double fnew = f(x);
			double xnew,xc,fc;
			
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
			
			/*if(verbose)
				std::cout << "----------------------- Ridder's Method iteration ---------------------------\n"
						<< "Initial guess: " << initial_guess << "\n"
						<< "Max iterations: " << max_iter << "\n"
						<< "Error tolerance: " << tolerance << "\n"
						<< "# iter.\tx\t\tf(x)\t\tf_x(x) \n";
            
            if (verbose) printf("%d\t%8.3e\t%8.2e\n",iterations_used,x,f(x));
			*/
			x = -10; // A saturation value always outside bracket
			for(int i = 0; i < max_iter; i++)
			{
				++iterations_used;
				xc = 0.5*(xa+xb);
				fc = f(xc);
				s = sqrt(fc*fc-fa*fb); // newt(2.4,fc*fc-fa*fb,6); //
				//std::cout << "Root: " << s << ", " << sqrt(fc*fc-fa*fb) << /*"," << f(0) << "," << f(1) << "," << f(xc)  << */ "\n"; // Note: The iterative method returns the reciprocal root!
				//return ErrorPolicy::handleTooManyIterationsNewton(x, max_iter, f(x));
				if (s == 0.0)
				{
					//if(verbose) std::cout << "----------------------- End Ridder's Method iteration ---------------------------\n";
					return x;
				}
				xnew = xc + (xc-xa)*((fa >= fb ? 1.0 : -1.0)*fc/s); // *s); // Consider changing 1/s to a multiplication
				if( std::abs(xnew-x) <= tolerance ) 
				{
					//if(verbose) std::cout << "----------------------- End Ridder's Method iteration ---------------------------\n";
					return x;
				}
				x = xnew;
				fnew = f(x);
				if( std::abs(fnew) <= tolerance )
				{
					//if(verbose) std::cout << "----------------------- End Ridder's Method iteration ---------------------------\n";
					return x;
				}
				if(fc*fnew < 0.0)
				{
					xa = xc; fa = fc;
					xb = x; fb = fnew;
				}
				else if(fa*fnew < 0.0)
				{
					xb = x; fb = fnew;
				}
				else if(fb*fnew < 0.0)
				{
					xa = x; fa = fnew;
				}
				else
					ErrorPolicy::handleBracketingFailure(xa,xb,fc,fnew);
					
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
					//if(verbose) std::cout << "----------------------- End Ridder's Method iteration ---------------------------\n";
					return x;
				}
				
				//if (verbose) printf("%d\t%8.3e\t%8.2e\t%8.2e \n",iterations_used,x,f(x),f.ds(x));
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
		inline static double solve(const Functor& f,
								   const double initial_guess,
								   const double visc_ratio,
								   const double inflec,
								   const int max_iter,
								   const double tolerance,
								   const bool verbose,
								   int& iterations_used, bool isTestRun, std::vector<std::pair<double,double>> & solution_path)
        {
            double x = initial_guess + 1 + 2*tolerance;
            double xNew = initial_guess;
            double fx = f(xNew);
            //double d2fw_inflec = visc_ratio*(2*inflec+1)*(inflec-1)*(inflec-1)+inflec*inflec*(2*inflec-3);
            //std::cout << "d2fw_inflec: " << d2fw_inflec << "\n";
            
            if(verbose)
            std::cout << "----------------------- Newton Trust Region iteration ---------------------------\n"
					  << "Initial guess: " << initial_guess << "\n"
                      << "Max iterations: " << max_iter << "\n"
                      << "Error tolerance: " << tolerance << "\n"
                      << "Inflection point: " << inflec << "\n"
                      //<< "Func. val. at inflec: " << d2fw_inflec << "\n"
                      << "# iter.\tx\t\tf(x)\t\tf_x(x) \n";
            
            /*if(verbose)
            {
				std::cout << "Viscosity ratio: " << visc_ratio << "\n";
	            double thePoint = 0.4;
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
				file.close();
			}*/
            
            if(isTestRun)
				addPointToVector(xNew,fx,solution_path); // REMOVE FOR SPEED TESTS
            
            while (fabs(fx) > tolerance && fabs(x-xNew) > tolerance)
            {
				++iterations_used;
				x = xNew;
				if (iterations_used > max_iter)
                    return ErrorPolicy::handleTooManyIterationsNewton(x, max_iter, f(x));
                
				xNew = x - fx/f.ds(x); // Scalar Newton method
				xNew = std::max(std::min(xNew,1.0),0.0); // Restrict to interval [0,1]
				if( (inflec-xNew)*(inflec-x) < 0 ) // Trust Region
					xNew = inflec;
				
				if (verbose)
					printf("%d\t%8.2e\t%8.2e\t%8.2e \n",iterations_used,xNew,f(xNew),f.ds(xNew));
				
				fx = f(xNew);
				if(isTestRun)
					addPointToVector(xNew,fx,solution_path); // REMOVE FOR SPEED TESTS
			}
			if(verbose)
				std::cout << "---------------------- End Newton iteration ------------------------\n";
			
            return xNew;
		}
		
		// Darcy flow trust region solver
		// Uses supplied viscosity ratio for trust region scheme
		template <class Functor>
		inline static double solveApprox(const Functor& f,
								   const double initial_guess,
								   const double visc_ratio,
								   const int max_iter,
								   const double tolerance,
								   const bool verbose,
								   int& iterations_used, bool isTestRun, std::vector<std::pair<double,double>> & solution_path)
        {
            double x = initial_guess + 1 + 2*tolerance;
            double xNew = initial_guess;
            double fx = f(xNew);
            
            if(isTestRun)
				addPointToVector(xNew,fx,solution_path); // REMOVE WHEN TESTING SPEED
            
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
            
            while (std::abs(fx) > tolerance && fabs(x-xNew) > tolerance)
            {
				++iterations_used;
				x = xNew;
				if (iterations_used > max_iter)
                    return ErrorPolicy::handleTooManyIterationsNewton(x, max_iter, f(x));
                
				// Scalar Newton's method
				xNew = x - fx/f.ds(x);
				
				// Restrict to interval [0,1]
				xNew = std::max(std::min(xNew,1.0),0.0);
				
				// Trust Region
				if(dfw2(xNew,visc_ratio)*dfw2(x,visc_ratio) < 0.0)
					xNew = (xNew+x)/2.0;
				
				if (verbose)
					printf("%d\t%8.2e\t%8.2e\t%8.2e \n",iterations_used,xNew,f(x),f.ds(x));
				
				fx = f(xNew);
				if(isTestRun)
					addPointToVector(xNew,fx,solution_path); // REMOVE WHEN TESTING SPEED
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
                                   int& iterations_used, bool isTestRun, std::vector<std::pair<double,double>> & solution_path)
        {
			
			if(isTestRun)
				addPointToVector(initial_guess,f(initial_guess),solution_path); // REMOVE FOR SPEED TESTS
            using namespace std;
            const double macheps = numeric_limits<double>::epsilon();
            const double eps = tolerance + macheps*max(max(fabs(a), fabs(b)), 1.0);
            double f_initial = f(initial_guess);
            const double epsF = tolerance + macheps*max(fabs(f_initial), 1.0);
            if (fabs(f_initial) < epsF) {
				if(isTestRun)
					addPointToVector(initial_guess,f_initial,solution_path); // REMOVE FOR SPEED TESTS
                return initial_guess;
            }
            double x0 = a;
            double x1 = b;
            double f0 = f_initial;
            double f1 = f_initial;
            if (x0 != initial_guess) {
                f0 = f(x0);
                if (fabs(f0) < epsF) {
					if(isTestRun)
						addPointToVector(x0,f0,solution_path); // REMOVE FOR SPEED TESTS
                    return x0;
                }
            }
            if (x1 != initial_guess) {
                f1 = f(x1);
                if (fabs(f1) < epsF) {
					if(isTestRun)
						addPointToVector(x1,f1,solution_path); // REMOVE FOR SPEED TESTS
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
                if(isTestRun)
					addPointToVector(xnew,fnew,solution_path); // REMOVE FOR SPEED TESTS
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
