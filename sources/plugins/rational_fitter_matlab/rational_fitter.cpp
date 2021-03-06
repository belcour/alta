/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2016 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include "rational_fitter.h"

#include <Eigen/Dense>
#include <Eigen/SVD>

#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

#define BUFFER_SIZE 10000

using namespace alta;

ALTA_DLL_EXPORT fitter* provide_fitter()
{
	return new rational_fitter_matlab();
}

rational_fitter_matlab::rational_fitter_matlab() 
{
	// Create matlab engine
#ifdef WIN32
    if (!(ep = engOpen(NULL)))
#else
    if (!(ep = engOpen("matlab -nosplash")))
#endif
	{
		std::cerr << "<ERROR>> can't start MATLAB engine" << std::endl ;
	}
}
rational_fitter_matlab::~rational_fitter_matlab() 
{
	engClose(ep); 
}

bool rational_fitter_matlab::fit_data(const ptr<data>& dat, ptr<function>& fit, const arguments &args)
{
	ptr<rational_function> r = dynamic_pointer_cast<rational_function>(fit) ;
	const ptr<vertical_segment> d = dynamic_pointer_cast<vertical_segment>(dat) ;
	if(!r || !d || ep == NULL
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
            time.stop() ;
			std::cout << "<<INFO>> got a fit using np = " << temp_np << " & nq =  " << temp_nq << "      " << std::endl ;
            std::cout << "<<INFO>> it took " << time << std::endl ;

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

void rational_fitter_matlab::set_parameters(const arguments& args)
{
	_max_np = args.get_float("np", 10) ;
	_max_nq = args.get_float("nq", 10) ;
	_min_np = args.get_float("min-np", _max_np) ;
	_min_nq = args.get_float("min-nq", _max_nq) ;	

    _use_matlab = !args.is_defined("use-qpas");
}
		
bool rational_fitter_matlab::fit_data(const ptr<vertical_segment>& d, int np, int nq, const ptr<rational_function>& r) 
{
    // For each output dimension (color channel for BRDFs) perform
    // a separate fit on the y-1D rational function.
    for(int j=0; j<d->dimY(); ++j)
    {
        rational_function_1d* rs = r->get(j);
        rs->resize(np, nq);

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
bool rational_fitter_matlab::fit_data(const ptr<vertical_segment>& d, int np, int nq, int ny, rational_function_1d* r)
{
	// Size of the problem
	int N = np+nq ;
	int M = d->size() ;
	
	// Matrices of the problem
	Eigen::MatrixXd G (N, N) ;
	Eigen::VectorXd g (N) ;
	Eigen::MatrixXd CI(2*M, N) ;
	Eigen::VectorXd ci(2*M) ;

	// Select the size of the result vector to
	// be equal to the dimension of p + q
	for(int j=0; j<N; ++j)
	{
		for(int i=0; i<N; ++i)
		{
			if(i == j)
				G(i, j) = 1.0 ; 
			else
				G(i,j) = 0.0;
		}

		g(j) = 0.0 ;
	}
	
	// Each constraint (fitting interval or point
	// add another dimension to the constraint
	// matrix
	for(int i=0; i<M; ++i)	
	{		
		// Norm of the row vector
		double a0_norm = 0.0 ;
		double a1_norm = 0.0 ;
		double a0_sum  = 0.0 ;
		double a1_sum  = 0.0 ;

		vec xi = d->get(i) ;

		// A row of the constraint matrix has this 
		// form: [p_{0}(x_i), .., p_{np}(x_i), -f(x_i) q_{0}(x_i), .., -f(x_i) q_{nq}(x_i)]
		// For the lower constraint and negated for 
		// the upper constraint
		for(int j=0; j<N; ++j)
		{
			// Filling the p part
			if(j<np)
			{
				const double pi = r->p(xi, j) ;
				
				// Updating Eigen matrix
				CI(2*i+0, j) =  pi;
				CI(2*i+1, j) = -pi;
			}
			// Filling the q part
			else
			{
				vec yl, yu; 
				d->get(i, yl, yu);

				const double qi = r->q(xi, j-np);
				
				// Updating Eigen matrix
				CI(2*i+0, j) = -yu[ny] * qi;
				CI(2*i+1, j) =  yl[ny] * qi;
			}

			// Update the norm of the row
			a0_norm += CI(2*i+0, j)*CI(2*i+0, j);
			a1_norm += CI(2*i+1, j)*CI(2*i+1, j);
			a0_sum  += CI(2*i+0, j);
			a1_sum  += CI(2*i+1, j);
		}
	
		// Set the c vector, will later be updated using the
		// delta parameter.
		ci(2*i+0) = -sqrt(a0_norm) ;
		ci(2*i+1) = -sqrt(a1_norm) ;

#ifdef REMOVE
		if(std::abs(a0_sum) < 1.0E-2) { std::cout << "Constraint " << 2*i+0 << " has a low sum: " << a0_sum << std::endl; }
		if(std::abs(a1_sum) < 1.0E-2) { std::cout << "Constraint " << 2*i+1 << " has a low sum: " << a1_sum << std::endl; }
#endif
	}
	
	// Update the ci column with the delta parameter
	// (See Celis et al. 2007 p.12)
	Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::HouseholderQRPreconditioner> svd(CI);
	const double sigma_m = svd.singularValues()(std::min(2*M, N)-1) ;
	const double sigma_M = svd.singularValues()(0) ;
	double delta = sigma_m / sigma_M ;

#ifdef NOT_WORKING
	static double cond_number = std::numeric_limits<double>::max() ;
	cond_number = std::min(delta, cond_number) ;
#endif 

	if(std::isnan(delta) || (std::abs(delta) == std::numeric_limits<double>::infinity()))
	{
#ifdef DEBUG_MATRICES
		std::cerr << "<<ERROR>> delta factor is NaN of Inf" << std::endl ;
#endif
		return false ;
	}
	else if(delta < 1.0E-06)
	{
		delta = 1.0 ;
	}

#ifndef DEBUG_MATRICES
	std::cout << "<<DEBUG>> delta factor: " << sigma_m << " / " << sigma_M << " = " << delta << std::endl ;
#endif
	
	for(int i=0; i<2*M; ++i)	
	{		
		ci(i) = ci(i) * delta ; 
	}
	
#ifdef DEBUG_MATRICES
	std::cout << "CI = " << CI << std::endl << std::endl ;
#endif

	// Create the MATLAB defintion of objects
	// MATLAB defines a quad prog as
	//   1/2 x' H x + f' x with A x <= b 
	//
	mxArray *H, *f, *A, *b, *x, *flag;
	H = mxCreateDoubleMatrix(  N, N, mxREAL);
	f = mxCreateDoubleMatrix(  N, 1, mxREAL);
	A = mxCreateDoubleMatrix(2*M, N, mxREAL);
	b = mxCreateDoubleMatrix(2*M, 1, mxREAL);

	memcpy((void *)mxGetPr(H), (void *) G.data(),  N*N*sizeof(double));
	memcpy((void *)mxGetPr(f), (void *) g.data(),  N*sizeof(double));
	memcpy((void *)mxGetPr(A), (void *)CI.data(), (2*M)*N*sizeof(double));
	memcpy((void *)mxGetPr(b), (void *)ci.data(), (2*M)*sizeof(double));

	engPutVariable(ep, "H", H);
	engPutVariable(ep, "f", f);
	engPutVariable(ep, "A", A);
	engPutVariable(ep, "b", b);
	
	char* output = new char[BUFFER_SIZE+1];
	output[BUFFER_SIZE] = '\0';
	engOutputBuffer(ep, output, BUFFER_SIZE) ;
#ifdef DEBUG_MATRICES
	engEvalString(ep, "display(H)");
	std::cout << output << std::endl ;
	engEvalString(ep, "display(f)");
	std::cout << output << std::endl ;
	engEvalString(ep, "display(A)");
	std::cout << output << std::endl ;
	engEvalString(ep, "display(b)");
	std::cout << output << std::endl ;
#endif

    // Use Matlab quadratic programming solver
    if(_use_matlab)
    {
        engEvalString(ep, "options = optimset('Algorithm', 'interior-point-convex')");
        engEvalString(ep, "[x, fval, flag] = quadprog(H, f, A, b, [], [], [], [], [], options);");
#ifdef DEBUG
        std::cout << output << std::endl ;
#endif
    }
    else // Use QPAS solver
    {
        engEvalString(ep, "cd matlab;");
        engEvalString(ep, "[x, err] = qpas(H,f,A,b);");
#ifdef DEBUG
        std::cout << output << std::endl ;
#endif
        engEvalString(ep, "flag = err == 0.0;");
        engEvalString(ep, "cd ..;");
    }

#ifdef DEBUG_MATRICES
	engEvalString(ep, "display(x)");
	std::cout << output << std::endl ;
	engEvalString(ep, "display(flag)");
    std::cout << output << std::endl ;
#endif

	mxDestroyArray(H);
	mxDestroyArray(f);
	mxDestroyArray(A);
	mxDestroyArray(b);

	x    = engGetVariable(ep, "x") ;
	flag = engGetVariable(ep, "flag") ;
	if(x != NULL)
	{
		if(flag != NULL)
		{
			const double flag_val = mxGetScalar(flag);
			if(flag_val < 0)
			{
				mxDestroyArray(x);
				mxDestroyArray(flag);
			
#ifdef DEBUG
				std::cerr << "<<ERROR>> flag is an error flag with no solution" << std::endl ;
#endif
				return false ;
			}
			else if(flag_val != 0)
			{
				std::cout << "<<INFO>> the solution might not be the optimal solution, this might be" << std::endl;
				std::cout << "<<INFO>> caused by local min, or under precision steps" << std::endl;
			}

			double* val = (double*)mxGetData(x) ;
			vec a(np), b(nq);
			for(int i=0; i<N; ++i)
			{
				if(i < np)
				{
                    a[i] = val[i] ;
				}
				else
				{
                    b[i-np] = val[i] ;
				}
			}
			r->update(a, b) ;

			mxDestroyArray(x);
			mxDestroyArray(flag);
            return true ;
		}
		else
		{
#ifdef DEBUG
			std::cerr << "<<ERROR>> unable to gather result flag" << std::endl ;
#endif
			return false ;
		}
	}
	else 
	{
#ifdef DEBUG
		std::cerr << "<<ERROR>> unable to gather result x" << std::endl ;
#endif
		return false ;
	}
}
