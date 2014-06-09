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
		return ( b >= 0.0 ? fabs(a) : -fabs(a) );
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
    
    template <class ErrorPolicy = ThrowOnError>
    class Brent
    {
		public:
		
		/// Brent's method for finding 
		/// roots of scalar equations as described in 
		/// Brent, 1973, "Algorithms for Minimization 
		/// Without Derivatives".
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
			int iter;
			double a=x1,b=x2,c=x2,d=a,e=a;
			double fa=f(a),fb=f(b),fc=f(initial_guess),tol_mod,m,eps = 1.0e-16;
			bool a_equals_c = false;
			
			if(verbose)
				std::cout << "----------------------- Brent's Method iteration ---------------------------\n"
						<< "Initial guess: " << initial_guess << "\n"
						<< "Max iterations: " << max_iter << "\n"
						<< "Error tolerance: " << tol << "\n"
						<< "# iter.\tx\t\tf(x) \n";
						
			if((fa <= 0.0) == (fb <= 0.0))
			{
				if (verbose) printf("%d\t%8.3e\t%8.3e \n",iterations_used,fabs(fa) <= tol ? a : b , fabs(fa) <= tol ? fa : fb );
				if ( fabs(fa) <= tol) return a;
				else if ( fabs(fb) <= tol) return b;
				else return ErrorPolicy::handleBracketingFailure(a,b,fa,fb);
			}
			if((fa < 0.0) == (fc < 0.0))
			{
				a = initial_guess; fa = fc;
			}
			else
			{
				b = initial_guess; fb = fc;
			}
			c = b;
			fc = fb;
			
			if(isTestRun)	
				addPointToVector(initial_guess,f(initial_guess),solution_path);
			
			for(iter = 1; iter <= max_iter; iter++)
			{
				++iterations_used;
				if( (fb > 0.0) == (fc > 0.0) ) 
				{
					a_equals_c = true;
					c = a; fc = fa;
					e = d = b-a;
				}
				if( fabs(fc) < fabs(fb) ) 
				{
					a_equals_c = true; 
					a=b; 	b=c; 	c=a;
					fa=fb; 	fb=fc; 	fc=fa;
				}
				tol_mod = 2.0*eps*fabs(b)+0.5*tol;
				m = 0.5*(c-b);
				if (verbose) printf("%d\t%8.3e\t%8.3e\n",iterations_used,b,f(b));
				if(fabs(m) <= tol_mod || fabs(fb) <= tol_mod) return b;
				if(fabs(e) < tol_mod || fabs(fa) <= fabs(fb)) 
				{
					d=m; e=m;
				}
				else
				{
					double p,q,r,s=fb/fa;
					if(a_equals_c)
					{
						p=2.0*m*s; q = 1.0-s;
					}
					else
					{
						q = fa/fc;
						r = fb/fc;
						p = s*(2.0*m*q*(q-r)-(b-a)*(r-1.0));
						q=(q-1.0)*(r-1.0)*(s-1.0);
					}
					if(p > 0.0) q = -q;
					else p = fabs(p); //else p = -p;
					
					double criteria = 3.0*m*q-fabs(tol_mod*q);
					double criteria_alt = fabs(e*q);
					criteria = (criteria < criteria_alt ? criteria : criteria_alt);
					if(2.0*p < criteria) 
					{
						e=d; d=p/q;
					}
					else
					{
						d=m; e=m;
					}
				}
				a=b; fa=fb;
				a_equals_c = false;
				if(fabs(d) > tol_mod)
					b += d;
				else
					b += (m > 0 ? tol_mod : -tol_mod);
				fb = f(b);
				if(isTestRun)
					addPointToVector(b,fb,solution_path);
			}
			
			return ErrorPolicy::handleTooManyIterationsNewton(a,b,max_iter);
		}
		
		/// Brent's method for finding roots of scalar 
		/// equations as described in Brent, 1973, 
		/// "Algorithms for Minimization Without Derivatives".
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
			
			if((fa <= 0.0) == (fb <= 0))
			{
				if ( std::abs(fa) <= tolerance) return xa;
				else if ( std::abs(fb) <= tolerance) return xb;
				else return ErrorPolicy::handleBracketingFailure(a,b,fa,fb);
			}
			
			if((fa < 0.0) == (fs < 0.0))
			{
				xb = initial_guess;
				fb = fs;
			}
			else if((fb < 0.0) == (fs < 0.0))
			{
				xa = initial_guess;
				fa = fs;
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
	};

	template <class ErrorPolicy = ThrowOnError>
	class Ridder
	{
		public:
		
		/// Ridders' method for finding 
		/// roots of scalar equations as described in 
		/// Ridders, 1973, "A new Algorithm for Computing 
		/// a Single Root of a Real Continuous Function"
		template <class Functor>
		inline static double solve(const Functor& f,
								   const double initial_guess,
								   double x1, double x2,
								   const int max_iter,
								   const double tolerance,
								   bool verbose,
								   int& iterations_used, 
								   bool isTestRun, 
								   std::vector<std::pair<double,double>> & solution_path)
		{
			int j;
			double ans,fh,fl,fm,fnew,s,xh,xl,xm,xnew=initial_guess,f_init=f(initial_guess);
			
			if(isTestRun)
				addPointToVector(initial_guess,f(initial_guess),solution_path);
			
			if(verbose)
				std::cout << "----------------------- Ridder's Method iteration ---------------------------\n"
						<< "Initial guess: " << initial_guess << "\n"
						<< "Max iterations: " << max_iter << "\n"
						<< "Error tolerance: " << tolerance << "\n"
						<< "# iter.\tx\t\tf(x) \n";
            
            if (verbose) printf("%d\t%8.3e\t%8.2e\n",iterations_used,xnew,f(xnew));
            
			fh = f(x2);
			//std::cout << "x2 : " << x2 << ", fh: " << fh << std::endl;
			if( (fh < 0.0) == (f_init < 0.0) )
			{
				fl = f(x1);
				fh = f_init;
				x2 = initial_guess;
			}
			else
			{
				fl = f_init;
				x1 = initial_guess;
			}
			
			if( (fl < 0.0) == (fh > 0.0) ) {
				xl = x1; xh = x2;
				ans = (xl+xh)*0.5;
				for(j = 1; j <= max_iter; j++) {
					++iterations_used;
					xm = 0.5*(xl+xh);
					fm = f(xm);
					s = sqrt(fm*fm-fl*fh);
					if(s == 0.0) return ans; 
					xnew = xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s);
					if(isTestRun)
						addPointToVector(xnew,f(xnew),solution_path);
					if (verbose) printf("%d\t%8.3e\t%8.2e\n",iterations_used,xnew,f(xnew));
					if (fabs(xnew-ans) <= tolerance) return xnew;
					ans = xnew;
					fnew = f(ans);
					if(fabs(fnew) <= tolerance) return ans;
					if((fm > 0.0) == (fnew < 0.0)) {
						xl = xm;  fl = fm;
						xh = ans; fh = fnew;
					}
					else if((fl < 0.0) == (fnew < 0.0)) {
						xl = ans; fl = fnew;
					}
					else if((fh < 0.0) == (fnew < 0.0)) {
						xh = ans; fh = fnew;
					}
					else {
						ErrorPolicy::handleBracketingFailure(xl,xh,fl,fh);
					}
					if(fabs(xh-xl) <= tolerance)
						return ans;
				}
				if(verbose) std::cout << "----------------------- End Ridder's Method iteration ---------------------------\n";
				return ErrorPolicy::handleTooManyIterationsNewton(xnew, max_iter, f(xnew));
			}
			else {
				if(isTestRun)
					addPointToVector(x1,f(x1),solution_path);
				if(fabs(fl) <= tolerance) return x1;
				if(isTestRun)
					addPointToVector(x2,f(x2),solution_path);
				if(fabs(fh) <= tolerance) return x2;
				return ErrorPolicy::handleBracketingFailure(x1,x2,fl,fh);
			}
			return 0.0; // Unreachable
		}
		
		/// Ridders' method for finding 
		/// roots of scalar equations as described in 
		/// Ridders, 1973, "A new Algorithm for Computing 
		/// a Single Root of a Real Continuous Function"
		/// Returns the root value f(x) in the reference variable fnew
		template <class Functor>
		inline static double solve(const Functor& f,
								   const double initial_guess,
								   double x1, double x2,
								   const int max_iter,
								   const double tolerance,
								   bool verbose,
								   int& iterations_used, 
								   bool isTestRun, 
								   std::vector<std::pair<double,double>> & solution_path, double & fnew)
		{
			double invalid_ans = -1;
			double ans,fh,fl,fm,s,xh,xl,xm,xnew=initial_guess,f_init=f(initial_guess);
			
			if(isTestRun)
				addPointToVector(initial_guess,f(initial_guess),solution_path);
			
			if(verbose)
				std::cout << "----------------------- Ridder's Method iteration ---------------------------\n"
						<< "Initial guess: " << initial_guess << "\n"
						<< "Max iterations: " << max_iter << "\n"
						<< "Error tolerance: " << tolerance << "\n"
						<< "# iter.\tx\t\tf(x) \n";
            
            if (verbose) printf("%d\t%8.3e\t%8.2e\n",iterations_used,xnew,f(xnew));
            
			fh = f(x2);
			if( (fh < 0.0) == (f_init < 0.0) )
			{
				fl = f(x1);
				fh = f_init;
				x2 = initial_guess;
			}
			else
			{
				fl = f_init;
				x1 = initial_guess;
			}
			
			if( (fl < 0.0) == (fh > 0.0) ) {
				xl = x1; xh = x2;
				ans = invalid_ans;
				for(int j = 1; j <= max_iter; j++) {
					++iterations_used;
					xm = 0.5*(xl+xh);
					fm = f(xm);
					s = sqrt(fm*fm-fl*fh);
					if(s == 0.0) return ans; 
					xnew = xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s);
					if(isTestRun)
						addPointToVector(xnew,f(xnew),solution_path);
					if (verbose) printf("%d\t%8.3e\t%8.2e\n",iterations_used,xnew,f(xnew));
					if (fabs(xnew-ans) <= tolerance) return xnew;
					ans = xnew;
					fnew = f(ans);
					if(fabs(fnew) <= tolerance) return ans;
					if((fm > 0.0) == (fnew < 0.0)) {
						xl = xm;  fl = fm;
						xh = ans; fh = fnew;
					}
					else if((fl < 0.0) == (fnew < 0.0)) {
						xl = ans; fl = fnew;
					}
					else if((fh < 0.0) == (fnew < 0.0)) {
						xh = ans; fh = fnew;
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
					addPointToVector(x1,f(x1),solution_path);
				if(fabs(fl) <= tolerance) return x1;
				if(isTestRun)
					addPointToVector(x2,f(x2),solution_path);
				if(fabs(fh) <= tolerance) return x2;
				return ErrorPolicy::handleBracketingFailure(x1,x2,fl,fh);
			}
			return 0.0; // Unreachable
		}
	};
	
    template <class ErrorPolicy = ThrowOnError>
    class NewtonRaphsonTrustRegion
    {
		public:
		
		/// Trust region solver
		/// Assumes quadratic relative permeability 
		/// functions in s to compute derivatives.
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
		
		/// Trust region solver
		/// Uses supplied viscosity ratio for trust region scheme
		/// Requires a precomputed f_w inflection point as argument.
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
            double dfx = -1.0;
            double fx = f(xNew,dfx);
            
            if(verbose)
            std::cout << "----------------------- Newton Trust Region iteration ---------------------------\n"
					  << "Initial guess: " << initial_guess << "\n"
                      << "Max iterations: " << max_iter << "\n"
                      << "Error tolerance: " << tolerance << "\n"
                      << "Inflection point: " << inflec << "\n"
                      << "# iter.\tx\t\tf(x)\t\tf_x(x) \n";
            
            for(int i = 0; i < max_iter; i++)
            {
				x = xNew;
				
				xNew = x - fx/dfx; // Scalar Newton method
				++iterations_used;
				if(xNew > 1.0) xNew = 1.0; // Restrict to domain
				else if(xNew < 0.0) xNew = 0.0;
				if( (xNew < inflec) == (x > inflec) ) // Trust Region
					xNew = inflec;
				
				if (verbose)
					printf("%d\t%8.2e\t%8.2e\t%8.2e \n",iterations_used,xNew,f(xNew),f.ds(xNew));
				
				if(fabs(xNew-x) < tolerance)
					return xNew;
				
				fx = f(xNew,dfx);
				
				if(fabs(fx) < tolerance)
					return xNew;
			}
			if( fabs(fx) < tolerance ) return xNew;
			return ErrorPolicy::handleTooManyIterationsNewton(x, max_iter, f(x));
		}
		
		/*/// Trust region solver
		/// Uses supplied viscosity ratio for trust region scheme
		/// Requires a precomputed f_w inflection point as argument.
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
            double dfx = -1.0;
            double fx = f(xNew,dfx);
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
            
            //if(verbose)
            //{
				//std::cout << "Viscosity ratio: " << visc_ratio << "\n";
	            //double thePoint = 0.4;
	            //std::string line;
	            //std::ifstream infile("pointdata.txt");
	            //if(infile.is_open())
	            //{
					//if( getline(infile,line) )
						//thePoint = atof(line.c_str());
					//infile.close();
				//}
				//double xval;
	            //double xmin = 0; double xmax = 1;
	            //int n_points = 150;
	            //std::ofstream file;
	            //file.open("cell_residual.txt");
	            //for ( int i = 0; i <= n_points; i++)
	            //{
					//xval = (xmax-xmin)/n_points*i;
					//file << xval << "\t" << f(xval) << "\t" << f.ds(xval) << "\t" << fw(xval,visc_ratio) << "\t" << dfw2(xval,visc_ratio) << "\t" << f(thePoint) + f.ds(thePoint)*(xval-thePoint) << "\n";
				//}
				//file.close();
			//}
            
            if(isTestRun)
				addPointToVector(xNew,fx,solution_path); // REMOVE FOR SPEED TESTS
            
            while (fabs(fx) > tolerance && fabs(x-xNew) > tolerance)
            {
				++iterations_used;
				x = xNew;
				if (iterations_used > max_iter)
				{
					//PrintFunctor<Functor>::printFunctorValues(f, 100, "residual_fail.data");
                    return ErrorPolicy::handleTooManyIterationsNewton(x, max_iter, f(x));
				}
                
				xNew = x - fx/dfx; // Scalar Newton method
				//xNew = std::max(std::min(xNew,1.0),0.0); // Restrict to interval [0,1]
				if(xNew > 1.0) xNew = 1.0;
				else if(xNew < 0.0) xNew = 0.0;
				//if( (inflec-xNew)*(inflec-x) < 0 ) // Trust Region
				if( (xNew < inflec) == (x > inflec) ) // Trust Region
					xNew = inflec;
				
				if (verbose)
					printf("%d\t%8.2e\t%8.2e\t%8.2e \n",iterations_used,xNew,f(xNew),f.ds(xNew));
				
				fx = f(xNew,dfx);
				if(isTestRun)
					addPointToVector(xNew,fx,solution_path); // REMOVE FOR SPEED TESTS
			}
			if(verbose)
				std::cout << "---------------------- End Newton iteration ------------------------\n";
			
            return xNew;
		}*/
		
		/// Trust region solver
		/// Uses supplied viscosity ratio for trust region scheme
		/// Assumes quadratic relative permeability 
		/// functions in s to compute derivatives.
		template <class Functor>
		inline static double solveApprox(const Functor& f,
								   const double initial_guess,
								   const double visc_ratio,
								   const int max_iter,
								   const double tolerance,
								   const bool verbose,
								   int& iterations_used, bool isTestRun, std::vector<std::pair<double,double>> & solution_path)
        {
            double x = initial_guess + 10*tolerance;
            double xNew = initial_guess;
            double dfx = -1.0;
            double fx = f(xNew,dfx);
            
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
            
            while (fabs(fx) > tolerance && fabs(x-xNew) > tolerance)
            {
				++iterations_used;
				x = xNew;
				if (iterations_used > max_iter)
                    return ErrorPolicy::handleTooManyIterationsNewton(x, max_iter, f(x));
                
				// Scalar Newton's method
				xNew = x - fx/dfx;
				
				// Restrict to interval [0,1]
				xNew = std::max(std::min(xNew,1.0),0.0);
				
				// Trust Region
				if( (dfw2(xNew,visc_ratio) < 0.0) != (dfw2(x,visc_ratio) < 0.0) )
					xNew = (xNew+x)/2.0;
				
				if (verbose)
					printf("%d\t%8.2e\t%8.2e\t%8.2e \n",iterations_used,xNew,f(x),f.ds(x));
				
				fx = f(xNew,dfx);
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
	};

    template <class ErrorPolicy = ThrowOnError>
    class NewtonRaphson
    {
		public:
		/// Newton Raphson solver
		template <class Functor>
		inline static double solve(const Functor& f,
								   const double initial_guess,
								   const int max_iter,
								   const double tolerance,
								   bool verbose,
								   int& iterations_used)
        {
            double x = initial_guess, xNew;
            
            if(verbose)
            std::cout << "----------------------- Newton iteration ---------------------------\n"
					  << "Initial guess: " << initial_guess << "\n"
                      << "Max iterations: " << max_iter << "\n"
                      << "Error tolerance: " << tolerance << "\n"
                      << "# iter.\tx\t\tf(x)\t\tf_x(x) \n";
            double dfxds = 0.0;
            double fx = f(x,dfxds);
            for(int i = 0; i < max_iter; i++, iterations_used++)
            {
				//xNew = x - fx/f.ds(x);
				xNew = x - fx/dfxds;
				fx = f(xNew,dfxds);
				if (verbose) printf("%d\t%8.3e\t%8.2e\t%8.2e \n",iterations_used,xNew,fx,dfxds);
				if(fabs(x-xNew) < tolerance || fabs(fx) < tolerance)
				{
					if(verbose) std::cout << "---------------------- End Newton iteration ------------------------\n";
					return xNew;
				}
				x = xNew;
			}
			
			return ErrorPolicy::handleTooManyIterationsNewton(x, max_iter, fx);
		}
		
		/// Newton Raphson solver
		/// If isTestRun is true, the solution 
		/// updates are returned in a vector
		template <class Functor>
		inline static double solve(const Functor& f,
								   const double initial_guess,
								   const double f_init,
								   const int max_iter,
								   const double tolerance,
								   bool verbose,
								   int& iterations_used,
								   bool isTestRun, std::vector<std::pair<double,double>> & solution_path)
        {
            double x = initial_guess, xNew;
            
            if(verbose)
            std::cout << "----------------------- Newton iteration ---------------------------\n"
					  << "Initial guess: " << initial_guess << "\n"
                      << "Max iterations: " << max_iter << "\n"
                      << "Error tolerance: " << tolerance << "\n"
                      << "# iter.\tx\t\tf(x)\t\tf_x(x) \n";
            double fx = f_init; //f(x);
            double dfds = f.ds(x);
            for(int i = 0; i < max_iter; i++, iterations_used++)
            {
				xNew = x - fx/dfds;
				fx = f(xNew,dfds);
				if (verbose) printf("%d\t%8.3e\t%8.2e\t%8.2e \n",iterations_used,xNew,fx,f.ds(xNew));
				if(isTestRun)
					addPointToVector(xNew,fx,solution_path); // REMOVE FOR SPEED TESTS
				if(fabs(x-xNew) < tolerance || fabs(fx) < tolerance)
				{
					if(verbose) std::cout << "---------------------- End Newton iteration ------------------------\n";
						return xNew;
				}
				x = xNew;
			}
			
			return ErrorPolicy::handleTooManyIterationsNewton(x, max_iter, fx);
		}
		/*template <class Functor>
		inline static double solve(const Functor& f,
								   const double initial_guess,
								   const double f_init,
								   const int max_iter,
								   const double tolerance,
								   bool verbose,
								   int& iterations_used,
								   bool isTestRun, std::vector<std::pair<double,double>> & solution_path)
        {
            double x = initial_guess, xNew;
            
            if(verbose)
            std::cout << "----------------------- Newton iteration ---------------------------\n"
					  << "Initial guess: " << initial_guess << "\n"
                      << "Max iterations: " << max_iter << "\n"
                      << "Error tolerance: " << tolerance << "\n"
                      << "# iter.\tx\t\tf(x)\t\tf_x(x) \n";
            double fx = f_init; //f(x);
            for(int i = 0; i < max_iter; i++, iterations_used++)
            {
				xNew = x - fx/f.ds(x);
				fx = f(xNew);
				if (verbose) printf("%d\t%8.3e\t%8.2e\t%8.2e \n",iterations_used,xNew,fx,f.ds(x));
				if(isTestRun)
					addPointToVector(xNew,fx,solution_path); // REMOVE FOR SPEED TESTS
				if(fabs(x-xNew) < tolerance || fabs(fx) < tolerance)
				{
					if(verbose) std::cout << "---------------------- End Newton iteration ------------------------\n";
						return xNew;
				}
				x = xNew;
			}
			
			return ErrorPolicy::handleTooManyIterationsNewton(x, max_iter, fx);
		}*/
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
                                   int& iterations_used, bool verbose, bool isTestRun, std::vector<std::pair<double,double>> & solution_path)
        {
			if(verbose)
            std::cout << "----------------------- Regula Falsi iteration ---------------------------\n"
					  << "Initial guess: " << initial_guess << "\n"
					  << "Bracket: [" << a << "," << b << "]\n"
                      << "Max iterations: " << max_iter << "\n"
                      << "Error tolerance: " << tolerance << "\n"
                      << "# iter.\tx\t\tf(x)\n";

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
                if (verbose) printf("%d\t%8.3e\t%8.2e \n",iterations_used+1,xnew,fnew);
                if(isTestRun)
					addPointToVector(xnew,fnew,solution_path); // REMOVE FOR SPEED TESTS
// 		cout << "xnew = " << xnew << "    fnew = " << fnew << endl;
                ++iterations_used;
                if (iterations_used > max_iter) {
                    //return ErrorPolicy::handleTooManyIterations(x0, x1, max_iter);
                    return ErrorPolicy::handleTooManyIterationsNewton(xnew, x1, max_iter);
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
            if(verbose) std::cout << "----------------------- Regula Falsi iteration ---------------------------\n";
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

	template <class ErrorPolicy = ThrowOnError>
	class GlobalizedNewton
	{
		public:
		
		template <class Functor>
		inline static double solve(const Functor& f,
								   const double initial_guess,
								   double x1, double x2,
								   const int max_iter,
								   const double tolerance,
								   bool verbose,
								   int& iterations_used,
								   bool isTestRun,
								   std::vector<std::pair<double,double>> & solution_path)
		{
			double fnew = 0.0;
			// Run a few steps of the initializer, e.g. Ridders' method
			double sat = Ridder<ContinueOnError>::solve(f,initial_guess,x1,x2,2,tolerance,verbose,iterations_used,isTestRun,solution_path,fnew);
			//double sat = RegulaFalsi<ContinueOnError>::solve(f,initial_guess,x1,x2,2,tolerance,iterations_used,verbose,isTestRun,solution_path);
			if(fabs(fnew) < tolerance)
				return sat;
			// Refine answer from initializer with Newton-Raphson
			sat = NewtonRaphson<ThrowOnError>::solve(f,sat,fnew,max_iter,tolerance,verbose,iterations_used, isTestRun, solution_path);
			
			return sat;
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
