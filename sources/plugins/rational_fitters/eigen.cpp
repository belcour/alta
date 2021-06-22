/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2013, 2014, 2016 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

// Eigen
#include <Eigen/Dense>
#include <Eigen/SVD>

// STL
#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

// ALTA interfaces
#include <core/function.h>
#include <core/data.h>
#include <core/fitter.h>
#include <core/args.h>
#include <core/rational_function.h>
#include <core/vertical_segment.h>
#include <core/common.h>

using namespace alta;

/*! \brief A least square fitter for rational function using the library Eigen
 *  \ingroup plugins
 *  \ingroup fitters
 *  \todo :  WRITE MORE
 */
class rational_fitter_eigen : public fitter
{
	private: // data

		// min and Max usable np and nq values for the fitting
		int _np, _nq;

	public: //methods

		rational_fitter_eigen() 
		{
		}
		~rational_fitter_eigen() 
		{
		}

		bool fit_data(const ptr<data>& dat, ptr<function>& fit, const arguments&)
		{
			ptr<rational_function> r = dynamic_pointer_cast<rational_function>(fit) ;
			if(!r)
			{
				std::cerr << "<<ERROR>> not passing the correct function object to the fitter" << std::endl ;
				return false ;
			} 
			
			ptr<vertical_segment> d = dynamic_pointer_cast<vertical_segment>(dat) ;
			if(!d
			|| d->confidence_interval_kind() != vertical_segment::ASYMMETRICAL_CONFIDENCE_INTERVAL)
			{
				std::cerr << "<<ERROR>> not passing the correct data object to the fitter" << std::endl ;
				return false ;
			}

		// XXX: FIT and D may have different values of dimX() and dimY(), but
		// this is fine: we convert values as needed in operator().
			r->setMin(d->min());
			r->setMax(d->max());

			std::cout << "<<INFO>> np =" << _np << "& nq =" << _nq  << std::endl ;
			r->setSize(_np, _nq);

			timer time ;
			time.start() ;

			if(fit_data(d, _np, _nq, r))
			{
				time.stop();
				std::cout << "<<INFO>> got a fit" << std::endl ;
				std::cout << "<<INFO>> it took " << time << std::endl ;

				return true ;
			}

			std::cout << "<<INFO>> fit failed\r"  ;
			std::cout.flush() ;

			return false ;
		}

		void set_parameters(const arguments& args)
		{
			_np = args.get_float("np", 10) ;
			_nq = args.get_float("nq", 10) ;
		}
		
	protected: // methods

		bool fit_data(const ptr<vertical_segment>& d, int np, int nq, const ptr<rational_function>& r) 
		{
			// For each output dimension (color channel for BRDFs) perform
			// a separate fit on the y-1D rational function.
			for(int j=0; j<d->parametrization().dimY(); ++j)
			{
				rational_function_1d* rs = r->get(j);
				if(!fit_data(d, np, nq, j, rs))
				{
					return false ;
				}
			}

			return true ;
		}

		// dat is the data object, it contains all the points to fit
		// np and nq are the degree of the Rational Polynomial Function to fit to the data
		// y is the dimension to fit on the y-data (e.g. R, G or B for RGB signals)
		// the function returns a rational BRDF function and a boolean
		bool fit_data(const ptr<vertical_segment>& d, int np, int nq, int ny, rational_function_1d* r)
		{
			// Each constraint (fitting interval or point
			// add another dimension to the constraint
			// matrix
			Eigen::MatrixXd CI(np+nq, d->size()) ;
			for(int i=0; i<d->size(); ++i)	
			{		
				// A row of the constraint matrix has this 
				// form: [p_{0}(x_i), .., p_{np}(x_i), -f(x_i) q_{0}(x_i), .., -f(x_i) q_{nq}(x_i)]
				for(int j=0; j<np+nq; ++j)
				{
					// Filling the p part
					if(j<np)
					{
						const double pi = r->p(d->get(i), j) ;
						CI(j, i) =  pi ;
					}
					// Filling the q part
					else
					{
						vec v = d->get(i) ;
						const double y  = v[d->parametrization().dimX() + ny] ;
						const double qi = r->q(d->get(i), j-np) ;
						CI(j, i) = -y * qi ;
					}
				}
			}

		//	std::cout << CI << std::endl << std::endl ;

			// Perform the Eigen decomposition of CI CI'
			Eigen::MatrixXd M(np+nq, np+nq) ;
			M.noalias() = CI * CI.transpose();
		// Faster alternative for large np+nq (only the lower triangular half gets computed):
		// M.setZero();
		// M.selfadjointView<Lower>().rankUpdate(CI.transpose());
			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(M) ;
		/*
			std::cout << M.size() << std::endl << std::endl ;
			std::cout << M << std::endl << std::endl ;
		*/
			if(solver.info() == Eigen::Success)
			{
		/*
				std::cout << solver.eigenvalues() << std::endl ;
		*/
				// Calculate the minimum eigen value and
				// its position
		#ifdef NOT_WORKING
				int    min_id = 0;
				double min_val = solver.eigenvalues().minCoeff(&min_id);
		#else
				int    min_id = 0 ;
				double min_val = std::numeric_limits<double>::max() ;
				for(int i=0; i<solver.eigenvalues().size(); ++i)
				{
					double value = solver.eigenvalues()[i] ;
					if(value >= 0 && value < min_val)
					{
						min_id  = i ;
						min_val = value ;
					}
				}
				
				if(min_val == std::numeric_limits<double>::max())
				{
					return false;
				}
		#endif
				// Recopy the vector d
				vec p(np), q(nq);
				Eigen::VectorXd::Map(&p[0], np) = solver.eigenvectors().col(min_id).head(np);
				Eigen::VectorXd::Map(&q[0], nq) = solver.eigenvectors().col(min_id).tail(nq);
				
				r->update(p, q) ;
				std::cout << "<<INFO>> got solution " << *r << std::endl ;
				return true;
			}
			else
			{
				return false ;
			}
		}
};

ALTA_DLL_EXPORT fitter* provide_fitter()
{
	return new rational_fitter_eigen();
}
