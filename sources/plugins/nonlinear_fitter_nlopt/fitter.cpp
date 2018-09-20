/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2016 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include "fitter.h"

#include <nlopt.h>

#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>
#include <cassert>

#include <core/common.h>

using namespace alta;

ALTA_DLL_EXPORT fitter* provide_fitter()
{
	return new nonlinear_fitter_nlopt();
}

void print_nlopt_error(nlopt_result res, const std::string& string)
{
	if(res == NLOPT_FAILURE)
	{
		std::cerr << "<<ERROR>> generic failure for \"" << string << "\"" << std::endl;
	}
	else if(res == NLOPT_INVALID_ARGS)
	{
		std::cerr << "<<ERROR>> invalid arguments for \"" << string << "\"" << std::endl;
	}
	else if(res == NLOPT_OUT_OF_MEMORY)
	{
		std::cerr << "<<ERROR>> not enough memory for \"" << string << "\"" << std::endl;
	}
}

// The parameter of the function _f should be set prior to this function
// call. If not it will produce undesirable results.
void df(double* fjac, const nonlinear_function* f, const data* d)
{
	// Clean memory
	memset(fjac, 0.0, f->nbParameters()*sizeof(double));

	// Each constraint is of the form data point * color channel
	for(int s=0; s<d->size(); ++s)
	{
		// Get the data sample and extract the value part
		vec xi = d->get(s);
		vec _di = vec(d->parametrization().dimY());
		for(int i=0; i<d->parametrization().dimY(); ++i)
		{
			_di[i] = xi[d->parametrization().dimX() + i];
		}
		
		// Convert the sample point into the function space
		vec x(f->parametrization().dimX());
		params::convert(&xi[0],
                    d->parametrization().input_parametrization(),
                    f->parametrization().input_parametrization(),
                    &x[0]);
 
		// Get the jacobian of the function at position x_i for the current
		// set of parameters (set prior to function call)
		vec _jac = f->parametersJacobian(x);

		// Should add the resulting vector completely
		vec _y = (*f)(x) - _di;

		// Fill the columns of the matrix
		for(int j=0; j<f->nbParameters(); ++j)
		{
			// For each output channel, update the subpart of the
			// vector row
			for(int i=0; i<f->parametrization().dimY(); ++i)
			{
				fjac[j] += 2 * _y[i] * _jac[i*f->nbParameters() + j];
			}
		}
	}
}

double f(unsigned n, const double* x, double* dy, void* dat)
{
	nonlinear_function* _f = (nonlinear_function*)(((void**)dat)[0]);
	const data* _d = (const data*)(((void**)dat)[1]);

	// Update the parameters vector
	vec _p(_f->nbParameters());
	for(int i=0; i<_f->nbParameters(); ++i) { _p[i] = x[i]; }
	_f->setParameters(_p);

	// Create the result vector
	double y = 0.0;

	// Each constraint is of the form data point * color channel
	for(int s=0; s<_d->size(); ++s)
	{
		// Get the data sample and extract the value part
		vec xi = _d->get(s);
		vec _di = vec(_d->parametrization().dimY());
		for(int i=0; i<_d->parametrization().dimY(); ++i)
		{
			_di[i] = xi[_d->parametrization().dimX() + i];
		}

		// Convert the sample point into the function space
		vec x(_f->parametrization().dimX());
		params::convert(&xi[0],
                    _d->parametrization().input_parametrization(),
                    _f->parametrization().input_parametrization(),
                    &x[0]);

		// Should add the resulting vector completely
		vec _y = (*_f)(x) - _di;
		for(int i=0; i<_f->parametrization().dimY(); ++i)
		{
			y += pow(_y[i], 2);
		}
	}

	if(dy != NULL)
	{
		df(dy, _f, _d);
	}

	return y;
}

nonlinear_fitter_nlopt::nonlinear_fitter_nlopt() 
{
}
nonlinear_fitter_nlopt::~nonlinear_fitter_nlopt() 
{
}

bool nonlinear_fitter_nlopt::fit_data(const ptr<data>& d, ptr<function>& fit, const arguments &args)
{
  // XXX: FIT and D may have different values of dimX() and dimY(), but
  // this is fine: we convert values as needed in operator().
	fit->setMin(d->min());
	fit->setMax(d->max());

	// Convert the function and bootstrap it with the data
	ptr<nonlinear_function> nf = dynamic_pointer_cast<nonlinear_function>(fit);
	if(!nf)
	{
		std::cerr << "<<ERROR>> Nl-Opt Fitter. the function is not a non-linear function" << std::endl;
		return false;
	}
	nf->bootstrap(d, args);

	#ifdef DEBUG
	std::cout << "<<DEBUG>> Nl-Opt Fitter. number of parameters: " << nf->nbParameters() << std::endl;
	#endif
	
	if(nf->nbParameters() == 0)
	{
		return true;
	}

	// the following starting values provide a rough fit is the bootstrap flag is
	// enabled
	vec p = nf->parameters();

	// Check the optimizer to use 
	nlopt_algorithm algorithm;
	nlopt_result res;
	std::string optimizerName = args["nlopt-optimizer"];
	if(optimizerName == "local-COBYLA")
	{
		algorithm = NLOPT_LN_COBYLA;
	}
	else if(optimizerName == "local-BOBYQA")
	{
		algorithm = NLOPT_LN_BOBYQA;
	}
	else if(optimizerName == "local-NEWUOA")
	{
		algorithm = NLOPT_LN_NEWUOA_BOUND;
	}
	else if(optimizerName == "local-PRAXIS")
	{
		algorithm = NLOPT_LN_PRAXIS;
	}
    else if(optimizerName == "local-Nelder-Mead")
	{
		algorithm = NLOPT_LN_NELDERMEAD;
	}
	else if(optimizerName == "local-spblx")
	{
		algorithm = NLOPT_LN_SBPLX;
	}
	else if(optimizerName == "global-crs")
	{
		algorithm = NLOPT_GN_CRS2_LM;
	}
	else if(optimizerName == "global-stogo")
	{
		algorithm = NLOPT_GD_STOGO;
	}
	else if(optimizerName == "global-isres")
	{
		algorithm = NLOPT_GN_ISRES;
	}
	else if(optimizerName == "local-mma")
	{
		algorithm = NLOPT_LD_MMA; 
	}
	else if(optimizerName == "local-sqp")
	{
		algorithm = NLOPT_LD_SLSQP;
	}
	else
	{
		algorithm = NLOPT_LD_SLSQP;
	}

	// Create the optimizer
	nlopt_opt opt = nlopt_create(algorithm, nf->nbParameters());
	if(opt == NULL)
	{
		std::cerr << "<<ERROR>> unable to create the optimizer" << std::endl;
		return false;
	}


	// Set the bounds to the parameters
	vec p_min = nf->getParametersMin();
	vec p_max = nf->getParametersMax();
	nlopt_set_lower_bounds(opt, &p_min[0]);
	nlopt_set_upper_bounds(opt, &p_max[0]);

	// Set the stopping criterion to a maximum number of evaluation
	if(args.is_defined("nlop-max-num-iterations"))
	{
		res = nlopt_set_maxeval(opt, args.get_int("nlop-max-num-iterations", 10));
		if(res < 0) { print_nlopt_error(res, "nlopt_set_maxeval"); }
	}

	//Cf. here : http://ab-initio.mit.edu/wiki/index.php/NLopt_Introduction#Termination_conditions

	//Set by default 1e-4 as relative function tolerance
	nlopt_set_xtol_rel(opt, 1e-4);

	//Parse  other options related to stopping criterion
	std::string const FUNC_TOLERANCE = "nlop-relative-function-tolerance";
	if( args.is_defined(FUNC_TOLERANCE) )
	{
			double const new_func_tol = args.get_double(FUNC_TOLERANCE,1e-4);			
			res = nlopt_set_xtol_rel(opt, new_func_tol );		
			
			if( res < 0 ) { print_nlopt_error(res, "nlopt_set_xtol_rel"); }
	}

	//Absolute Tolerance Function on a value
	std::string const ABS_FUNC_TOLERANCE = "nlop-abs-function-tolerance";
	if( args.is_defined(ABS_FUNC_TOLERANCE))
	{
		double const new_func_tol = args.get_double(ABS_FUNC_TOLERANCE, 1e-6);
		res= nlopt_set_ftol_abs(opt, new_func_tol);
		
		if( res < 0) { print_nlopt_error(res, "nlopt_set_ftol_abs"); }
	}


	// Create the problem
	void* dat[2];
	dat[0] = (void*)nf.get();
	dat[1] = (void*)d.get();
	res = nlopt_set_min_objective(opt, f, dat);
	if(res < 0)
	{
		print_nlopt_error(res, "nlopt_set_min_objective");
		return false;
	}
	
	// Launch minimization
	double f_end;
	std::cout << p << std::endl;
	res = nlopt_optimize(opt, &p[0], &f_end);

	if(res > 0)
	{
		std::cout << "<<INFO>> optimized distance: " << f_end << std::endl;
		std::cout << "<<INFO>> found parameters: " << p << std::endl;
		nf->setParameters(p);
	}
	else
	{
		print_nlopt_error(res, "nlopt_optimize");
	}

	nlopt_destroy(opt);
	return res > 0;
}

void nonlinear_fitter_nlopt::set_parameters(const arguments& args)
{
}
