/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2016 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include "rational_fitter_cgal.h"

#include <CGAL/basic.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>
#include <CGAL/MP_Float.h>
#include <Eigen/SVD>

#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

typedef CGAL::MP_Float ET ;
typedef CGAL::Quadratic_program<ET> Program ;
typedef CGAL::Quadratic_program_solution<ET> Solution ;
typedef CGAL::Quadratic_program_options Options ;

using namespace alta;

ALTA_DLL_EXPORT fitter* provide_fitter()
{
    return new rational_fitter_cgal();
}

rational_fitter_cgal::rational_fitter_cgal()
{
}
rational_fitter_cgal::~rational_fitter_cgal() 
{
}

bool rational_fitter_cgal::fit_data(const ptr<data>& dat, ptr<function>& fit, const arguments&)
{
	ptr<rational_function> r = dynamic_pointer_cast<rational_function>(fit) ;
	const ptr<vertical_segment> d = dynamic_pointer_cast<vertical_segment>(dat) ;
	if(!r || !d
     || d->confidence_interval_kind() != vertical_segment::ASYMMETRICAL_CONFIDENCE_INTERVAL)
	{
		std::cerr << "<<ERROR>> not passing the correct class to the fitter" << std::endl ;
		return false ;
	}

	// I need to set the dimension of the resulting function to be equal
	// to the dimension of my fitting problem
	r->setDimX(d->dimX()) ;
	r->setDimY(d->dimY()) ;
	r->setMin(d->min()) ;
	r->setMax(d->max()) ;

	std::cout << "<<INFO>> np in  [" << _min_np << ", " << _max_np 
	          << "] & nq in [" << _min_nq << ", " << _max_nq << "]" << std::endl ;


	int temp_np = _min_np, temp_nq = _min_nq ;
	while(temp_np <= _max_np || temp_nq <= _max_nq)
	{
		timer time ;
		time.start() ;

		r->setSize(temp_np, temp_nq);
		
		if(fit_data(d, temp_np, temp_nq, r))
		{
			time.stop();
			std::cout << "<<INFO>> got a fit using np = " << temp_np 
			          << " & nq =  " << temp_nq << "      " << std::endl
			          << "<<INFO>> it took " << time  << std::endl ;

			return true ;
		}

		std::cout << "<<INFO>> fit using np = " << temp_np << " & nq =  " << temp_nq << " failed\r"  ;
		std::cout.flush() ;

		if(temp_np <= _max_np)
		{
			++temp_np ;
		}
		if(temp_nq <= _max_nq)
		{
			++temp_nq ;
		}
	}

	return false ;
}

void rational_fitter_cgal::set_parameters(const arguments& args)
{
	_max_np = args.get_float("np", 10) ;
	_max_nq = args.get_float("nq", 10) ;
	_min_np = args.get_float("min-np", _max_np) ;
	_min_nq = args.get_float("min-nq", _max_nq) ;
}
		
bool rational_fitter_cgal::fit_data(const ptr<vertical_segment>& d, int np, int nq, const ptr<rational_function>& r) 
{
    // For each output dimension (color channel for BRDFs) perform
    // a separate fit on the y-1D rational function.
    for(int j=0; j<d->dimY(); ++j)
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
// np and nq are the degree of the RP to fit to the data
// y is the dimension to fit on the y-data (e.g. R, G or B for RGB signals)
// the function return a ration BRDF function and a boolean
bool rational_fitter_cgal::fit_data(const ptr<vertical_segment>& d, int np, int nq, int ny, rational_function_1d* r)
{
	// by default, we have a nonnegative QP with Ax - b >= 0
	Program qp (CGAL::LARGER, false, 0, false, 0) ; 

	// Select the size of the result vector to
	// be equal to the dimension of p + q
	for(int i=0; i<np+nq; ++i)
	{
		qp.set_d(i, i, 1.0) ; 
	}
	
	// Each constraint (fitting interval or point
	// add another dimension to the constraint
	// matrix
	Eigen::MatrixXd CI(np+nq, 2*d->size()) ;
	Eigen::VectorXd ci(2*d->size()) ;
	for(int i=0; i<d->size(); ++i)	
	{		
		// Norm of the row vector
		double a0_norm = 0.0 ;
		double a1_norm = 0.0 ;

		vec xi = d->get(i) ;
		
		// A row of the constraint matrix has this 
		// form: [p_{0}(x_i), .., p_{np}(x_i), -f(x_i) q_{0}(x_i), .., -f(x_i) q_{nq}(x_i)]
		// For the lower constraint and negated for 
		// the upper constraint
		for(int j=0; j<np+nq; ++j)
		{
			// Filling the p part
			if(j<np)
			{
				const double pi = r->p(xi, j) ;
				a0_norm += pi*pi ;
				a1_norm += pi*pi ;
				qp.set_a(j, i,  ET(pi)) ;
				qp.set_a(j, i+d->size(), ET(-pi)) ;
				CI(j, i) =  pi ;
				CI(j, i+d->size()) = -pi ;
			}
			// Filling the q part
			else
			{
				vec yl, yu ; 
				d->get(i, yl, yu) ;

				const double qi = r->q(xi, j-np) ;
				a0_norm += qi*qi * (yl[ny]*yl[ny]) ;
				qp.set_a(j, i, ET(-yl[ny] * qi)) ;
				CI(j, i) = -yl[ny] * qi ;
				
				a1_norm += qi*qi * (yu[ny]*yu[ny]) ;
				qp.set_a(j, i+d->size(),  ET(yu[ny] * qi)) ;
				CI(j, i+d->size()) = yu[ny] * qi ;
			}
		}
	
		// Set the c vector, will later be updated using the
		// delta parameter.
		ci(i) = sqrt(a0_norm) ;
		ci(i+d->size()) = sqrt(a1_norm) ;
	}
#ifdef DEBUG
	for(int j=0; j<d->size()*2; ++j)
	{
		for(int i=0; i<np+nq; ++i)
			std::cout << CI(i,j) << "\t";
		std::cout << std::endl; 
	}
	std::cout << std::endl ;
#endif
	// Update the ci column with the delta parameter
	// (See Celis et al. 2007 p.12)
	Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::HouseholderQRPreconditioner> svd(CI);
	const double sigma_m = svd.singularValues()(std::min(2*d->size(), np+nq)-1) ;
	const double sigma_M = svd.singularValues()(0) ;

#ifdef DEBUG
	std::cout << "<<DEBUG>> SVD = [ " ;
	for(int i=0; i<std::min(2*d->size(), np+nq); ++i)
	{
		std::cout << svd.singularValues()(i) << ", " ;
	}
	std::cout << " ]" << std::endl ;
#endif
	
	double delta = sigma_m / sigma_M ;
	if(std::isnan(delta) || (std::abs(delta) == std::numeric_limits<double>::infinity()))
	{
#ifdef DEBUG
		std::cerr << "<<ERROR>> delta factor is NaN of Inf" << std::endl ;
#endif
		return false ;
	}
	else if(delta == 0.0)
	{
		delta = 1.0 ;
	}

#ifdef DEBUG
	std::cout << "<<DEBUG>> delta factor: " << sigma_m << " / " << sigma_M << " = " << delta << std::endl ;
#endif
	for(int i=0; i<2*d->size(); ++i)	
	{		
		qp.set_b(i, ET(delta * ci(i))) ;
	}

#ifdef DEBUG
	// Export some informations on the problem to solve
	std::cout << "<<DEBUG>> " << qp.get_n() << " variables" << std::endl ;
	std::cout << "<<DEBUG>> " << qp.get_m() << " constraints" << std::endl ;
#endif

	if(qp.is_linear() && !qp.is_nonnegative())
	{
		std::cerr << "<<ERROR>> the problem should not be linear or negative!" << std::endl ;
		return false ;
	}

#ifdef EXTREM_DEBUG
	// Iterate over the rows
	std::cout << std::endl ;
	std::cout << "A = [" ;
	for(int j=0; j<qp.get_m(); ++j)
	{
		if(j == 0) std::cout << "   " ;

		// Iterate over the columns
		for(int i=0; i<qp.get_n(); ++i)
		{
			if(i > 0) std::cout << ",\t" ;
			std::cout << *(*(qp.get_a()+i) +j) ;
		}
		std::cout << ";" ;
		if(j < qp.get_n()-1) std::cout << std::endl ;
	}
	std::cout << "]" << std::endl << std::endl ;
	
	std::cout << std::endl ;
	std::cout << "D = [" ;
	for(int j=0; j<np+nq; ++j)
	{
		if(j == 0) std::cout << "   " ;

		// Iterate over the columns
		for(int i=0; i<np+nq; ++i)
		{
			if(i > 0) std::cout << ",\t" ;
			std::cout << *(*(qp.get_d()+i) +j) ;
		}
		std::cout << ";" ;
	}
	std::cout << "]" << std::endl << std::endl ;
#endif

	// solve the program, using ET as the exact type
	Options  o ;
	o.set_auto_validation(true) ;
	Solution s = CGAL::solve_quadratic_program(qp, ET(), o) ;

#ifdef DEBUG
	if(s.is_infeasible())
	{
		std::cout << "<<DEBUG>> the current program is infeasible" << std::endl ;
	}
#endif

	bool solves_qp = !s.is_infeasible() && s.solves_quadratic_program(qp)  ;
	for(int i=0; i<np+nq; ++i)
	{
		const double v = (double)CGAL::to_double(*(s.variable_numerators_begin()+i)) ;
		solves_qp = solves_qp && !std::isnan(v) && (v != std::numeric_limits<double>::infinity()) ;
	}

	if(solves_qp)
	{
		// Recopy the vector d
		vec p, q;
		p = vec::Zero(np) ; q = vec::Zero(nq) ;
		for(int i=0; i<np+nq; ++i)
		{
			const double v = CGAL::to_double(*(s.variable_numerators_begin()+i)) ;

			if(i < np)
			{
				p[i] = v ;
			}
			else
			{
				q[i-np] = v ;
			}
		}
 		r->update(p, q) ;
#ifdef DEBUG
        std::cout << "<<INFO>> got solution " << *r << std::endl ;
#endif
		
		return true;
	}
	else
	{
		return false; 
	}
}
