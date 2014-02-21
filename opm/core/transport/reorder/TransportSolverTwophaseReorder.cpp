/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

#include "config.h"
#include <opm/core/transport/reorder/TransportSolverTwophaseReorder.hpp>
#include <opm/core/props/IncompPropertiesInterface.hpp>
#include <opm/core/grid.h>
#include <opm/core/transport/reorder/reordersequence.h>
#include <opm/core/grid/ColumnExtract.hpp>
#include <opm/core/utility/RootFinders.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/pressure/tpfa/trans_tpfa.h>

#include <iostream>
#include <fstream>
#include <iterator>
#include <numeric>

#include <opm/core/utility/ErrorMacros.hpp>
#define _USE_MATH_DEFINES
#include <math.h>
//#include <complex>

#define EXPERIMENT_GAUSS_SEIDEL


namespace Opm
{

    // Choose error policy for scalar solves here.
    typedef RegulaFalsi<WarnAndContinueOnError> RootFinder;
    //typedef NewtonRaphson<WarnAndContinueOnError> RootFinder;
	
	TransportSolverTwophaseReorder::TransportSolverTwophaseReorder(const UnstructuredGrid& grid,
                                                                   const Opm::IncompPropertiesInterface& props,
                                                                   const double* gravity,
                                                                   const double tol,
                                                                   const int maxit,
                                                                   char solver_type,
                                                                   bool verbose,
                                                                   bool solver_flag)
        : grid_(grid),
          props_(props),
          tol_(tol),
          maxit_(maxit),
          solver_type_(solver_type),
          solver_flag_(solver_flag),
          darcyflux_(0),
          source_(0),
          dt_(0.0),
          saturation_(grid.number_of_cells, -1.0),
          fractionalflow_(grid.number_of_cells, -1.0),
          fractionalflowderivative_(grid.number_of_cells, -1.0),
          reorder_iterations_(grid.number_of_cells, 0),
          mob_(2*grid.number_of_cells, -1.0),
          dmob_(2*grid.number_of_cells, -1.0)
#ifdef EXPERIMENT_GAUSS_SEIDEL
        , ia_upw_(grid.number_of_cells + 1, -1),
          ja_upw_(grid.number_of_faces, -1),
          ia_downw_(grid.number_of_cells + 1, -1),
          ja_downw_(grid.number_of_faces, -1)
#endif
    {
        if (props.numPhases() != 2) {
            OPM_THROW(std::runtime_error, "Property object must have 2 phases");
        }
        visc_ = props.viscosity();
        int num_cells = props.numCells();
        smin_.resize(props.numPhases()*num_cells);
        smax_.resize(props.numPhases()*num_cells);
        std::vector<int> cells(num_cells);
        for (int i = 0; i < num_cells; ++i) {
            cells[i] = i;
        }
        props.satRange(props.numCells(), &cells[0], &smin_[0], &smax_[0]);
        if (gravity) {
            initGravity(gravity);
            initColumns();
        }
        setVerbose(verbose);
        initInflectionPoint(visc_[0]/visc_[1]);
        if(getVerbose())
			std::cout << "Inflection point " << flux_func_inflection_point_ << " initialized using viscosity ratio " << visc_[0]/visc_[1] << "\n";
    }

	TransportSolverTwophaseReorder::TransportSolverTwophaseReorder(const UnstructuredGrid& grid,
                                                                   const Opm::IncompPropertiesInterface& props,
                                                                   const double* gravity,
                                                                   const double tol,
                                                                   const int maxit,
                                                                   char solver_type,
                                                                   bool solver_flag)
        : grid_(grid),
          props_(props),
          tol_(tol),
          maxit_(maxit),
          solver_type_(solver_type),
          solver_flag_(solver_flag),
          darcyflux_(0),
          source_(0),
          dt_(0.0),
          saturation_(grid.number_of_cells, -1.0),
          fractionalflow_(grid.number_of_cells, -1.0),
          fractionalflowderivative_(grid.number_of_cells, -1.0),
          reorder_iterations_(grid.number_of_cells, 0),
          mob_(2*grid.number_of_cells, -1.0),
          dmob_(2*grid.number_of_cells, -1.0)
#ifdef EXPERIMENT_GAUSS_SEIDEL
        , ia_upw_(grid.number_of_cells + 1, -1),
          ja_upw_(grid.number_of_faces, -1),
          ia_downw_(grid.number_of_cells + 1, -1),
          ja_downw_(grid.number_of_faces, -1)
#endif
    {
        if (props.numPhases() != 2) {
            OPM_THROW(std::runtime_error, "Property object must have 2 phases");
        }
        visc_ = props.viscosity();
        int num_cells = props.numCells();
        smin_.resize(props.numPhases()*num_cells);
        smax_.resize(props.numPhases()*num_cells);
        std::vector<int> cells(num_cells);
        for (int i = 0; i < num_cells; ++i) {
            cells[i] = i;
        }
        props.satRange(props.numCells(), &cells[0], &smin_[0], &smax_[0]);
        if (gravity) {
            initGravity(gravity);
            initColumns();
        }
        initInflectionPoint(visc_[0]/visc_[1]);
    }

    TransportSolverTwophaseReorder::TransportSolverTwophaseReorder(const UnstructuredGrid& grid,
                                                                   const Opm::IncompPropertiesInterface& props,
                                                                   const double* gravity,
                                                                   const double tol,
                                                                   const int maxit)
        : grid_(grid),
          props_(props),
          tol_(tol),
          maxit_(maxit),
          darcyflux_(0),
          source_(0),
          dt_(0.0),
          saturation_(grid.number_of_cells, -1.0),
          fractionalflow_(grid.number_of_cells, -1.0),
          fractionalflowderivative_(grid.number_of_cells, -1.0),
          reorder_iterations_(grid.number_of_cells, 0),
          mob_(2*grid.number_of_cells, -1.0),
          dmob_(2*grid.number_of_cells, -1.0)
#ifdef EXPERIMENT_GAUSS_SEIDEL
        , ia_upw_(grid.number_of_cells + 1, -1),
          ja_upw_(grid.number_of_faces, -1),
          ia_downw_(grid.number_of_cells + 1, -1),
          ja_downw_(grid.number_of_faces, -1)
#endif
    {
        if (props.numPhases() != 2) {
            OPM_THROW(std::runtime_error, "Property object must have 2 phases");
        }
        visc_ = props.viscosity();
        int num_cells = props.numCells();
        smin_.resize(props.numPhases()*num_cells);
        smax_.resize(props.numPhases()*num_cells);
        std::vector<int> cells(num_cells);
        for (int i = 0; i < num_cells; ++i) {
            cells[i] = i;
        }
        props.satRange(props.numCells(), &cells[0], &smin_[0], &smax_[0]);
        if (gravity) {
            initGravity(gravity);
            initColumns();
        }
        solver_type_ = 'r';
        solver_flag_ = false;
    }


    TransportSolverTwophaseReorder::~TransportSolverTwophaseReorder()
    {
    }


    void TransportSolverTwophaseReorder::solve(const double* porevolume,
                                               const double* source,
                                               const double dt,
                                               TwophaseState& state)
    {
        darcyflux_ = &state.faceflux()[0];
        porevolume_ = porevolume;
        source_ = source;
        dt_ = dt;
        toWaterSat(state.saturation(), saturation_);

#ifdef EXPERIMENT_GAUSS_SEIDEL
        std::vector<int> seq(grid_.number_of_cells);
        std::vector<int> comp(grid_.number_of_cells + 1);
        int ncomp;
        compute_sequence_graph(&grid_, darcyflux_,
                               &seq[0], &comp[0], &ncomp,
                               &ia_upw_[0], &ja_upw_[0]);
        const int nf = grid_.number_of_faces;
        std::vector<double> neg_darcyflux(nf);
        std::transform(darcyflux_, darcyflux_ + nf, neg_darcyflux.begin(), std::negate<double>());
        compute_sequence_graph(&grid_, &neg_darcyflux[0],
                               &seq[0], &comp[0], &ncomp,
                               &ia_downw_[0], &ja_downw_[0]);
#endif
        std::fill(reorder_iterations_.begin(),reorder_iterations_.end(),0);
        reorderAndTransport(grid_, darcyflux_);
        toBothSat(saturation_, state.saturation());
    }


    const std::vector<int>& TransportSolverTwophaseReorder::getReorderIterations() const
    {
        return reorder_iterations_;
    }


    // Residual function r(s) for a single-cell implicit Euler transport
    //
    //     r(s) = s - s0 + dt/pv*( influx + outflux*f(s) )
    //
    // where influx is water influx, outflux is total outflux.
    // Influxes are negative, outfluxes positive.
    struct TransportSolverTwophaseReorder::Residual
    {
        int cell;
        double s0;
        double influx;    // sum_j min(v_ij, 0)*f(s_j) + q_w
        double dinflux;   // sum_j min(v_ij, 0)*dfds(s_j) + q_w
        double outflux;   // sum_j max(v_ij, 0) - q
        double comp_term; // q - sum_j v_ij
        double dtpv;    // dt/pv(i)
        const TransportSolverTwophaseReorder& tm;
        explicit Residual(const TransportSolverTwophaseReorder& tmodel, int cell_index)
            : tm(tmodel)
        {
            // Initialize [in|out]flux to include effect of transport source terms.
            cell    = cell_index;
            s0      = tm.saturation_[cell];
            double src_flux       = -tm.source_[cell];
            bool src_is_inflow = src_flux < 0.0;
            influx  =  src_is_inflow ? src_flux : 0.0;
            dinflux = 0.0;
            outflux = !src_is_inflow ? src_flux : 0.0;
            dtpv    = tm.dt_/tm.porevolume_[cell];

            // Compute fluxes over interior edges. Boundary flow is supposed to be
            // included in the transport source term, along with well sources.
            for (int i = tm.grid_.cell_facepos[cell]; i < tm.grid_.cell_facepos[cell+1]; ++i) {
                int f = tm.grid_.cell_faces[i];
                double flux;
                int other;
                // Compute cell flux
                if (cell == tm.grid_.face_cells[2*f]) {
                    flux  = tm.darcyflux_[f];
                    other = tm.grid_.face_cells[2*f+1];
                } else {
                    flux  =-tm.darcyflux_[f];
                    other = tm.grid_.face_cells[2*f];
                }
                // Add flux to influx or outflux, if interior.
                if (other != -1) {
                    if (flux < 0.0) {
                        influx  += flux*tm.fractionalflow_[other];
                        dinflux += flux*tm.fractionalflowderivative_[other];
                    } else {
                        outflux += flux;
                    }
                }
            }
            
            /*if(cell == 399)
            {
	            double xval;
	            double xmin = 0; double xmax = 1;
	            int n_points = 150;
	            std::ofstream file;
	            file.open("fractional_flow.txt");
	            for ( int i = 0; i <= n_points; i++)
	            {
					xval = (xmax-xmin)/n_points*i;
					file << xval << "\t" << tm.fracFlow(xval,cell) << "\t" << tm.fracFlowDerivative(xval,cell) << "\n";
				}
				file.close();
			}*/

        }
        double operator()(double s) const
        {
            return s - s0 + dtpv*(outflux*tm.fracFlow(s, cell) + influx);
        }
        double ds(double s) const
        {
			//return 1 + dtpv*(dinflux + outflux*tm.fracFlowDerivative(s,cell));
			return 1 + dtpv*(outflux*tm.fracFlowDerivative(s,cell));
		}
    };


    void TransportSolverTwophaseReorder::solveSingleCell(const int cell)
    {
        Residual res(*this, cell);
        // const double r0 = res(saturation_[cell]);
        // if (std::fabs(r0) < tol_) {
        //     return;
        // }
        int iters_used = 0;
        // saturation_[cell] = modifiedRegulaFalsi(res, smin_[2*cell], smax_[2*cell], maxit_, tol_, iters_used);
        
        //saturation_[cell] = RootFinder::solve(res, saturation_[cell], 0.0, 1.0, maxit_, tol_, iters_used); // Original. Commented 04.02.14 - Svein
        
		double inflec = getInflectionPoint();
        
        if(solver_type_ == 'n')
			saturation_[cell] = NewtonRaphson<ThrowOnError>::solve(res, saturation_[cell], 0.0, 1.0, maxit_, tol_, iters_used);
        else if(solver_type_ == 't')
        {
			double M = visc_[0]/visc_[1]; // Viscosity ratio, mu_w/mu_o for Trust Region scheme
			if(!solver_flag_)
			{
				if(getVerbose())
					std::cout << "Using precise Newton Raphson Trust Region method in cell " << cell << "\n";
				saturation_[cell] = NewtonRaphsonTrustRegion<ThrowOnError>::solve(res, saturation_[cell], M, inflec, maxit_, tol_, getVerbose(), iters_used);
			}
			else
			{
				if(getVerbose())
					std::cout << "Using approximate Newton Raphson Trust Region method in cell " << cell << "\n";
				/*if(cell == 25)
				{
					std::cout << "Using approximate Newton Raphson Trust Region method in cell " << cell << "\n";
					saturation_[cell] = NewtonRaphsonTrustRegion<ThrowOnError>::solveApprox(res, saturation_[cell], M, maxit_, tol_, true, iters_used);
				}
				else*/
				saturation_[cell] = NewtonRaphsonTrustRegion<ThrowOnError>::solveApprox(res, saturation_[cell], M, maxit_, tol_, getVerbose(), iters_used);
			}
		}
		else if(solver_type_ == 'u')
        {
			double M = visc_[0]/visc_[1]; // Viscosity ratio, mu_w/mu_o for Trust Region scheme
			saturation_[cell] = RegulaFalsiTrustRegion<ThrowOnError>::solve(res, saturation_[cell], 0.0, 1.0, M, maxit_, tol_, getVerbose(), iters_used);
		}
		else if(solver_type_ == 'i')
		{
			saturation_[cell] = Ridder<ThrowOnError>::solve(res, saturation_[cell], 0.0, 1.0, maxit_, tol_, iters_used);
		}
		else if(solver_type_ == 'b')
		{
			saturation_[cell] = Brent<ThrowOnError>::solve(res, saturation_[cell], 0.0, 1.0, maxit_, tol_, iters_used);
		}
        else
			saturation_[cell] = RootFinder::solve(res, saturation_[cell], 0.0, 1.0, maxit_, tol_, iters_used);
			
		/*if(iters_used == 0)
			std::cout << "No iterations used?! f(s):" << res(saturation_[cell]) << std::endl;
        else
			std::cout << "Iterations: " << iters_used << std::endl;*/
        
        // add if it is iteration on an out loop
        reorder_iterations_[cell] = reorder_iterations_[cell] + iters_used;
        fractionalflow_[cell] = fracFlow(saturation_[cell], cell);
        fractionalflowderivative_[cell] = fracFlowDerivative(saturation_[cell], cell);
    }

    void TransportSolverTwophaseReorder::solveMultiCell(const int num_cells, const int* cells)
    {

#ifdef EXPERIMENT_GAUSS_SEIDEL
        // Experiment: when a cell changes more than the tolerance,
        //             mark all downwind cells as needing updates. After
        //             computing a single update in each cell, use marks
        //             to guide further updating. Clear mark in cell when
        //             its solution gets updated.
        // Verdict: this is a good one! Approx. halved total time.
        std::vector<int> needs_update(num_cells, 1);
        // This one also needs the mapping from all cells to
        // the strongly connected subset to filter out connections
        std::vector<int> pos(grid_.number_of_cells, -1);
        for (int i = 0; i < num_cells; ++i) {
            const int cell = cells[i];
            pos[cell] = i;
        }

        // Note: partially copied from below.
        const double tol = 1e-9;
        const int max_iters = 300;
        // Must store s0 before we start.
        std::vector<double> s0(num_cells);
        // Must set initial fractional flows before we start.
        // Also, we compute the # of upstream neighbours.
        for (int i = 0; i < num_cells; ++i) {
            const int cell = cells[i];
            fractionalflow_[cell] = fracFlow(saturation_[cell], cell);
            fractionalflowderivative_[cell] = fracFlowDerivative(saturation_[cell], cell);
            s0[i] = saturation_[cell];
        }
        // Solve once in each cell.
        int num_iters = 0;
        int update_count = 0; // Change name/meaning to cells_updated?
        do {
            update_count = 0; // Must reset count for every iteration.
            for (int i = 0; i < num_cells; ++i) {
                if (!needs_update[i]) {
                    continue;
                }
                ++update_count;
                const int cell = cells[i];
                const double old_s = saturation_[cell];
                saturation_[cell] = s0[i];
                solveSingleCell(cell);
                const double s_change = std::fabs(saturation_[cell] - old_s);
                if (s_change > tol) {
                    // Mark downwind cells.
                    for (int j = ia_downw_[cell]; j < ia_downw_[cell+1]; ++j) {
                        const int downwind_cell = ja_downw_[j];
                        int ci = pos[downwind_cell];
                        if (ci != -1) {
                            needs_update[ci] = 1;
                        }
                    }
                }
                // Unmark this cell.
                needs_update[i] = 0;
            }
        } while (update_count > 0 && ++num_iters < max_iters);

        // Done with iterations, check if we succeeded.
        if (update_count > 0) {
            OPM_THROW(std::runtime_error, "In solveMultiCell(), we did not converge after "
                  << num_iters << " iterations. Remaining update count = " << update_count);
        }
        std::cout << "Solved " << num_cells << " cell multicell problem in "
                  << num_iters << " iterations." << std::endl;

#else
        double max_s_change = 0.0;
        const double tol = 1e-9;
        const int max_iters = 300;
        int num_iters = 0;
        // Must store s0 before we start.
        std::vector<double> s0(num_cells);
        // Must set initial fractional flows before we start.
        for (int i = 0; i < num_cells; ++i) {
            const int cell = cells[i];
            fractionalflow_[cell] = fracFlow(saturation_[cell], cell);
            fractionalflowderivative__[cell] = fracFlowDerivative(saturation_[cell], cell);
            s0[i] = saturation_[cell];
        }
        do {
            max_s_change = 0.0;
            for (int i = 0; i < num_cells; ++i) {
                const int cell = cells[i];
                const double old_s = saturation_[cell];
                saturation_[cell] = s0[i];
                solveSingleCell(cell);
                double s_change = std::fabs(saturation_[cell] - old_s);
                // std::cout << "cell = " << cell << "    delta s = " << s_change << std::endl;
                if (max_s_change < s_change) {
                    max_s_change = s_change;
                }
            }
            // std::cout << "Iter = " << num_iters << "    max_s_change = " << max_s_change
            //        << "    in cell " << max_change_cell << std::endl;
        } while (max_s_change > tol && ++num_iters < max_iters);
        if (max_s_change > tol) {
            OPM_THROW(std::runtime_error, "In solveMultiCell(), we did not converge after "
                  << num_iters << " iterations. Delta s = " << max_s_change);
        }
        std::cout << "Solved " << num_cells << " cell multicell problem in "
                  << num_iters << " iterations." << std::endl;
#endif // EXPERIMENT_GAUSS_SEIDEL
    }

    double TransportSolverTwophaseReorder::fracFlow(double s, int cell) const
    {
        double sat[2] = { s, 1.0 - s };
        double mob[2];
        props_.relperm(1, sat, &cell, mob, 0);
        mob[0] /= visc_[0];
        mob[1] /= visc_[1];
        return mob[0]/(mob[0] + mob[1]);
    }
    
    void TransportSolverTwophaseReorder::initInflectionPoint(const double M)
    {
		if(getVerbose())
			std::cout << "Computing fractional flow function inflection point using viscosity ratio " << M << " ...\n";
		
		std::vector<double> roots;
		computeRoots(roots,M);
		
		std::vector<double>::const_iterator iter = roots.begin()-1;
		if(getVerbose())
		{
			std::cout << "--- Found " << roots.size() << " roots ---\n";
			while(++iter != roots.end())
				std::cout << "x = " << *iter << "\n";
			std::cout << "--- End roots ---\n";
		}
		
		bool found = false;
		iter = roots.begin();
		while(!found && iter != roots.end())
		{
			if(checkRange(*iter))
			{
				flux_func_inflection_point_ = *iter;
				found = true;
			}
			iter++;
		}
		if(found == false)
			OPM_THROW(std::runtime_error,"Could not find the flow function inflection point.\n");
		if(getVerbose())
			std::cout << "Inflection point found at saturation " << flux_func_inflection_point_ << "\n";
	}
	void TransportSolverTwophaseReorder::computeRoots(std::vector<double> & roots, double M)
	{
		std::complex<double> up = computeU(M,true); 
		std::complex<double> um = computeU(M,false);
		//if(getVerbose())
		//	std::cout << "u_p: " << up << ", u_m: " << um << "\n";
		computeX(roots,up);
		computeX(roots,um);
	}
	std::complex<double> TransportSolverTwophaseReorder::computeU(double M, bool positive)
    {
		double t1,t2,t3,t4;
		t1 = 1.0/8.0; t2 = -0.25*M/(M+1); t3 = pow(0.5*M/(M+1)-0.25, 2)-1.0/16.0; t4 = 0.5*pow(fabs(t3), 0.5);
		if(getVerbose() && t3 > 0)
			std::cout << "Warning: Positive b^2-4ac, something smells fishy!\n";
		if(!positive)
			t4 = -t4;
		std::complex<double> var(t1+t2,t4);
		return var;
	}
	void TransportSolverTwophaseReorder::computeX(std::vector<double> & roots, std::complex<double> u)
	{
		if(getVerbose())
			std::cout << "Computing complex roots z from z^3 = " << u << "\n";
		double precision = 1e-8;
		
		double r = sqrt( u.real()*u.real() + u.imag()*u.imag() );
		double theta = M_PI;
		if(checkTarget(u.real(),precision))
			theta = M_PI*0.5;
		else
			theta = atan2(u.imag(),u.real());
			
		double sqrt3r = pow(r,1.0/3); 
		
		double innerTrig = computeInnerTrigArguments(theta,3,0);
		double real = sqrt3r*cos(innerTrig); double imag = sqrt3r*sin(innerTrig);
		std::complex<double> z1(real,imag);
		
		innerTrig = computeInnerTrigArguments(theta,3,1);
		real = sqrt3r*cos(innerTrig); imag = sqrt3r*sin(innerTrig);
		std::complex<double> z2(real,imag);
		
		innerTrig = computeInnerTrigArguments(theta,3,2);
		real = sqrt3r*cos(innerTrig); imag = sqrt3r*sin(innerTrig);
		std::complex<double> z3(real,imag);
		
		if(getVerbose())
			std::cout << "Found three complex roots z: " << z1 << ", " << z2 << ", " << z3 << "\n";
		
		z1 = computeY(z1); z2 = computeY(z2); z3 = computeY(z3); // Reusing z variables instead of making complex x and y vars
		
		if(getVerbose())
			std::cout << "Complex y-values from roots z: " << z1 << ", " << z2 << ", " << z3 << "\n";
		
		if(checkTarget(z1.imag(),precision))
			roots.push_back(computeX(z1.real()));
		if(checkTarget(z2.imag(),precision))
			roots.push_back(computeX(z2.real()));
		if(checkTarget(z3.imag(),precision))
			roots.push_back(computeX(z3.real()));
	}
	std::complex<double> TransportSolverTwophaseReorder::computeY(std::complex<double> z)
    {
		if(z.imag() == 0.0 && z.real() == 0.0)
		{
			OPM_THROW(std::runtime_error,"Can't compute y-value from a zero-valued u\n");
			return 1;
		}
		return z + 0.25/z;
	}
	// Utilities
	bool TransportSolverTwophaseReorder::checkTarget(double val, double target, double precision)
	{
		return val <= target + precision && val >= target - precision;
	}
	bool TransportSolverTwophaseReorder::checkTarget(double val, double precision)
	{
		return checkTarget(val, 0, precision);
	}
	bool TransportSolverTwophaseReorder::checkRange(double s)
    {
		return (s <= 1.0 && s >= 0.0);
	}
	double TransportSolverTwophaseReorder::computeInnerTrigArguments(double theta, double n, double k)
	{
		return theta/n + 2*M_PI*k/n;
	}
	// Real functions
	void TransportSolverTwophaseReorder::computeRootsFromSignCases(std::vector<double> & roots, double u)
	{
		if(u >= 0)
			roots.push_back(computeX(computeY(computeZ(u))));
		else
			computeX(roots,u);
	}
	void TransportSolverTwophaseReorder::computeX(std::vector<double> & roots, double u)
	{
		if(getVerbose())
			std::cout << "Computing complex roots from " << u << "\n";
		double precision = 1e-6;
		double theta = M_PI; double r = fabs(u);
		double sqrt3r = pow(r,1.0/3); 
		
		double innerTrig = computeInnerTrigArguments(theta,3,0);
		double real = sqrt3r*cos(innerTrig); double imag = sqrt3r*sin(innerTrig);
		std::complex<double> z1(real,imag);
		
		innerTrig = computeInnerTrigArguments(theta,3,1);
		real = sqrt3r*cos(innerTrig); imag = sqrt3r*sin(innerTrig);
		std::complex<double> z2(real,imag);
		
		innerTrig = computeInnerTrigArguments(theta,3,2);
		real = sqrt3r*cos(innerTrig); imag = sqrt3r*sin(innerTrig);
		std::complex<double> z3(real,imag);
		
		if(getVerbose())
			std::cout << "Complex roots: " << z1 << ", " << z2 << ", " << z3 << "\n";
		
		z1 = computeY(z1); z2 = computeY(z2); z3 = computeY(z3); // Reusing z variables instead of making complex x and y vars
		
		if(z1.imag() < precision && z1.imag() > -precision)
			roots.push_back(z1.real());
		if(z2.imag() < precision && z2.imag() > -precision)
			roots.push_back(z2.real());
		if(z2.imag() < precision && z3.imag() > -precision)
			roots.push_back(z3.real());
	}
	double TransportSolverTwophaseReorder::computeU(double M, bool positive,bool dummy)
    {
		double t1,t2,t3;
		t1 = 1.0/8.0; t2 = -0.25*M/(M+1); t3 = 0.5*pow(pow(0.5*M/(M+1)-0.25, 2)-1.0/16.0, 0.5);
		if(!positive)
			t3 = -t3;
		return t1+t2+t3;
	}
	double TransportSolverTwophaseReorder::computeZ(double u)
    {
		if(u < 0)
		{
			OPM_THROW(std::runtime_error,"Can't compute z-value from negative u\n");
			return 1;
		}
		return pow(u,1.0/3);
	}
	double TransportSolverTwophaseReorder::computeY(double z)
    {
		if(z == 0)
		{
			OPM_THROW(std::runtime_error,"Can't compute y-value from a zero-valued u\n");
			return 1;
		}
		return z + 0.25/z;
	}
	double TransportSolverTwophaseReorder::computeX(double y)
    {
		return y + 0.5;
	}
    
    void TransportSolverTwophaseReorder::initInflectionPoint_old(const double M)
    {
		/*double MP1 = M + 1;
		double M2,M3,M4,M5;
		M2 = M*M; M3 = M2*M; M4 = M3*M; M5 = M4*M;
		double inner_sqrt = 2*sqrt(-M5-4*M4-6*M3-4*M2-M);
		double alpha = 2*pow(-M3-M2+inner_sqrt+MP1,1.0/3.0); // Third root
		double beta = MP1/alpha;
		flux_func_inflection_point_ = 0.5*(beta + 1/beta + 1);*/
		
		double up,um;
		up = computeU(M,true,false); um = computeU(M,false,false);
		
		if(up >= 0 && um >= 0)
		{
			double zp = computeZ(up); double zm = computeZ(um);
			double xp,xm;
			if(zp >= 0 && zm >= 0)
			{
				xp = computeX(computeY(zp)); 
				xm = computeX(computeY(zm)); 
				if(checkRange(xp))
					flux_func_inflection_point_ = xp;
				else if(checkRange(xm))
					flux_func_inflection_point_ = xm;
				else
					OPM_THROW(std::runtime_error,"No inflection point found in range with two positive z-values\n");
			}
			else if(zp >= 0)
			{
				flux_func_inflection_point_ = computeX(computeY(zp));
				if(!checkRange(flux_func_inflection_point_))
					OPM_THROW(std::runtime_error,"Inflection point " << flux_func_inflection_point_ << " outside range using z_p\n");
			}
			else if(zm >= 0)
			{
				flux_func_inflection_point_ = computeX(computeY(zm));
				if(!checkRange(flux_func_inflection_point_))
					OPM_THROW(std::runtime_error,"Inflection point " << flux_func_inflection_point_ << " outside range using z_m\n");
			}
			else
				OPM_THROW(std::runtime_error, "No inflection point found due to negative z-values\n");
		}
		else if(up >= 0)
		{
			double z = computeZ(up);
			if(z >= 0)
			{
				flux_func_inflection_point_ = computeX(computeY(z));
				if(!checkRange(flux_func_inflection_point_))
					OPM_THROW(std::runtime_error,"Inflection point " << flux_func_inflection_point_ << " outside range using u_p\n");
			}
			else
				OPM_THROW(std::runtime_error,"No inflection point found due to negative z using u_p\n");
		}
		else if(um >= 0)
		{
			double z = computeZ(um);
			if(z >= 0)
			{
				flux_func_inflection_point_ = computeX(computeY(z));
				if(!checkRange(flux_func_inflection_point_))
					OPM_THROW(std::runtime_error,"Inflection point " << flux_func_inflection_point_ << " outside range using u_m\n");
			}
			else
				OPM_THROW(std::runtime_error,"No inflection point found due to negative z using u_m\n");
		}
		else
			OPM_THROW(std::runtime_error,"No inflection point found due to negative u-values This should never happen!\n");
		
		//double u = std::max(up,um);
		//double z = pow(u,1.0/3.0);
		//y = z + 1/(4*z);
		//flux_func_inflection_point_ = y + 1/2;
		std::cout << "Inflec. point fresh: " << flux_func_inflection_point_ << "\n";
	}
	
	double TransportSolverTwophaseReorder::getInflectionPoint() 
	{
		return flux_func_inflection_point_;
	}
    
    double TransportSolverTwophaseReorder::fracFlowDerivative(double s, int cell) const
    {
		double sat[2] = {s, 1.0 - s};
		double mob[2]; 
		double dmob[4];
		double smob, sdmob;
		
		props_.relperm(1, sat, &cell, mob, dmob);
		
		mob[0] /= visc_[0]; mob[1] /= visc_[1];
		smob =  mob[0] + mob[1];
		
		dmob[0] /= visc_[0]; dmob[3] /= -visc_[1];
		sdmob =  dmob[0] + dmob[3];
		//std::cout << "Mob.: " << mob[0] << ", " << mob[1] << "\ndMob: " << dmob[0] << ", " << dmob[3] << std::endl;
		return (dmob[0]*smob - mob[0]*sdmob)/(smob*smob);
	}

    // Residual function r(s) for a single-cell implicit Euler gravity segregation
    //
    //     r(s) = s - s0 + dt/pv*sum_{j adj i}( gravmod_ij * gf_ij ).
    //
    struct TransportSolverTwophaseReorder::GravityResidual
    {
        int cell;
        int nbcell[2];
        double s0;
        double dtpv;    // dt/pv(i)
        double gf[2];
        const TransportSolverTwophaseReorder& tm;
        explicit GravityResidual(const TransportSolverTwophaseReorder& tmodel,
                                 const std::vector<int>& cells,
                                 const int pos,
                                 const double* gravflux) // Always oriented towards next in column. Size = colsize - 1.
            : tm(tmodel)
        {
            cell = cells[pos];
            nbcell[0] = -1;
            gf[0] = 0.0;
            if (pos > 0) {
                nbcell[0] = cells[pos - 1];
                gf[0] = -gravflux[pos - 1];
            }
            nbcell[1] = -1;
            gf[1] = 0.0;
            if (pos < int(cells.size() - 1)) {
                nbcell[1] = cells[pos + 1];
                gf[1] = gravflux[pos];
            }
            s0      = tm.saturation_[cell];
            dtpv    = tm.dt_/tm.porevolume_[cell];

        }
        double operator()(double s) const
        {
            double res = s - s0;
            double mobcell[2];
            tm.mobility(s, cell, mobcell,0);
            for (int nb = 0; nb < 2; ++nb) {
                if (nbcell[nb] != -1) {
                    double m[2];
                    if (gf[nb] < 0.0) {
                        m[0] = mobcell[0];
                        m[1] = tm.mob_[2*nbcell[nb] + 1];
                    } else {
                        m[0] = tm.mob_[2*nbcell[nb]];
                        m[1] = mobcell[1];
                    }
                    if (m[0] + m[1] > 0.0) {
                        res += -dtpv*gf[nb]*m[0]*m[1]/(m[0] + m[1]);
                    }
                }
            }
            return res;
        }
        double ds(double s) const
        {
			double res = 1;
            double mobcell[2];
            double dmobcell[2];
            tm.mobility(s, cell, mobcell, dmobcell);
            for (int nb = 0; nb < 2; ++nb) {
                if (nbcell[nb] != -1) {
                    double m[2];
                    double dm[2];
                    if (gf[nb] < 0.0) {
                        m[0] = mobcell[0];
                        m[1] = tm.mob_[2*nbcell[nb] + 1];
                        dm[0] = -dmobcell[0]; // Minus sign handles the dkrds 'bug' in SaturationPropsBasic::relperm
                        dm[1] = tm.dmob_[2*nbcell[nb] + 1];
                    } else {
                        m[0] = tm.mob_[2*nbcell[nb]];
                        m[1] = mobcell[1];
                        dm[0] = tm.dmob_[2*nbcell[nb]];
                        dm[1] = dmobcell[1];
                    }
                    double msum = m[0] + m[1];
                    if (msum > 0.0) {
                        //res += -dtpv*gf[nb]*m[0]*m[1]/msum;
                        res += -dtpv*gf[nb]*( (dm[0]*m[0]+m[0]*dm[1])*msum - (dm[0]+dm[1])*m[0]*m[1] )/( msum*msum );
                    }
                }
            }
            return res;
		}
    };

    void TransportSolverTwophaseReorder::mobility(double s, int cell, double* mob, double* dmob) const
    {
        double sat[2] = { s, 1.0 - s };
        double dmobtemp[4];
        props_.relperm(1, sat, &cell, mob, dmobtemp);
        mob[0] /= visc_[0];
        mob[1] /= visc_[1];
        
        if(dmob != 0)
        {
			dmob[0] = dmobtemp[0]/visc_[0];
			dmob[1] = dmobtemp[3]/visc_[1];
		}
    }

    void TransportSolverTwophaseReorder::initGravity(const double* grav)
    {
        // Set up gravflux_ = T_ij g (rho_w - rho_o) (z_i - z_j)
        std::vector<double> htrans(grid_.cell_facepos[grid_.number_of_cells]);
        const int nf = grid_.number_of_faces;
        const int dim = grid_.dimensions;
        gravflux_.resize(nf);
        tpfa_htrans_compute(const_cast<UnstructuredGrid*>(&grid_), props_.permeability(), &htrans[0]);
        tpfa_trans_compute(const_cast<UnstructuredGrid*>(&grid_), &htrans[0], &gravflux_[0]);
        const double delta_rho = props_.density()[0] - props_.density()[1];
        for (int f = 0; f < nf; ++f) {
            const int* c = &grid_.face_cells[2*f];
            double gdz = 0.0;
            if (c[0] != -1 && c[1] != -1) {
                for (int d = 0; d < dim; ++d) {
                    gdz += grav[d]*(grid_.cell_centroids[dim*c[0] + d] - grid_.cell_centroids[dim*c[1] + d]);
                }
            }
            gravflux_[f] *= delta_rho*gdz;
        }
    }

    void TransportSolverTwophaseReorder::initColumns()
    {
        extractColumn(grid_, columns_);
    }

    void TransportSolverTwophaseReorder::solveSingleCellGravity(const std::vector<int>& cells,
                                                                const int pos,
                                                                const double* gravflux)
    {
        const int cell = cells[pos];
        GravityResidual res(*this, cells, pos, gravflux);
        if (std::fabs(res(saturation_[cell])) > tol_) {
            int iters_used = 0;
            saturation_[cell] = RootFinder::solve(res, smin_[2*cell], smax_[2*cell], maxit_, tol_, iters_used);
            reorder_iterations_[cell] = reorder_iterations_[cell] + iters_used;
        }
        saturation_[cell] = std::min(std::max(saturation_[cell], smin_[2*cell]), smax_[2*cell]);
        mobility(saturation_[cell], cell, &mob_[2*cell], &dmob_[2*cell]);
    }

    int TransportSolverTwophaseReorder::solveGravityColumn(const std::vector<int>& cells)
    {
        // Set up column gravflux.
        const int nc = cells.size();
        std::vector<double> col_gravflux(nc - 1);
        for (int ci = 0; ci < nc - 1; ++ci) {
            const int cell = cells[ci];
            const int next_cell = cells[ci + 1];
            for (int j = grid_.cell_facepos[cell]; j < grid_.cell_facepos[cell+1]; ++j) {
                const int face = grid_.cell_faces[j];
                const int c1 = grid_.face_cells[2*face + 0];
                const int c2 = grid_.face_cells[2*face + 1];
                if (c1 == next_cell || c2 == next_cell) {
                    const double gf = gravflux_[face];
                    col_gravflux[ci] = (c1 == cell) ? gf : -gf;
                }
            }
        }

        // Store initial saturation s0
        s0_.resize(nc);
        for (int ci = 0; ci < nc; ++ci) {
            s0_[ci] = saturation_[cells[ci]];
        }

        // Solve single cell problems, repeating if necessary.
        double max_s_change = 0.0;
        int num_iters = 0;
        do {
            max_s_change = 0.0;
            for (int ci = 0; ci < nc; ++ci) {
                const int ci2 = nc - ci - 1;
                double old_s[2] = { saturation_[cells[ci]],
                                    saturation_[cells[ci2]] };
                saturation_[cells[ci]] = s0_[ci];
                solveSingleCellGravity(cells, ci, &col_gravflux[0]);
                saturation_[cells[ci2]] = s0_[ci2];
                solveSingleCellGravity(cells, ci2, &col_gravflux[0]);
                max_s_change = std::max(max_s_change, std::max(std::fabs(saturation_[cells[ci]] - old_s[0]),
                                                               std::fabs(saturation_[cells[ci2]] - old_s[1])));
            }
            // std::cout << "Iter = " << num_iters << "    max_s_change = " << max_s_change << std::endl;
        } while (max_s_change > tol_ && ++num_iters < maxit_);

        if (max_s_change > tol_) {
            OPM_THROW(std::runtime_error, "In solveGravityColumn(), we did not converge after "
                  << num_iters << " iterations. Delta s = " << max_s_change);
        }
        return num_iters + 1;
    }

    void TransportSolverTwophaseReorder::solveGravity(const double* porevolume,
                                                      const double dt,
                                                      TwophaseState& state)
    {
        // Initialize mobilities.
        const int nc = grid_.number_of_cells;
        std::vector<int> cells(nc);
        for (int c = 0; c < nc; ++c) {
            cells[c] = c;
        }
        mob_.resize(2*nc);
        dmob_.resize(2*nc);
        double dmobtemp[4*nc];
        props_.relperm(cells.size(), &state.saturation()[0], &cells[0], &mob_[0], dmobtemp);
        const double* mu = props_.viscosity();
        for (int c = 0; c < nc; ++c) {
            mob_[2*c] /= mu[0];
            mob_[2*c + 1] /= mu[1];
            dmob_[2*c] = dmobtemp[4*c + 0]/mu[0];
            dmob_[2*c+1] = dmobtemp[4*c + 3]/mu[1];
        }

        // Set up other variables.
        porevolume_ = porevolume;
        dt_ = dt;
        toWaterSat(state.saturation(), saturation_);

        // Solve on all columns.
        int num_iters = 0;
        for (std::vector<std::vector<int> >::size_type i = 0; i < columns_.size(); i++) {
            // std::cout << "==== new column" << std::endl;
            num_iters += solveGravityColumn(columns_[i]);
        }
        std::cout << "Gauss-Seidel column solver average iterations: "
                  << double(num_iters)/double(columns_.size()) << std::endl;

        toBothSat(saturation_, state.saturation());
    }

} // namespace Opm



/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
