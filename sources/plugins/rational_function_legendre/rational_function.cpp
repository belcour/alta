/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2016 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include "rational_function.h"

#include <core/common.h>
#include <core/rational_function.h>

using namespace alta;

ALTA_DLL_EXPORT function* provide_function(const parameters& params)
{
    return new rational_function_legendre(params);
}

rational_function_legendre_1d::rational_function_legendre_1d()
{
}

rational_function_legendre_1d::rational_function_legendre_1d(int nX, int np, int nq) :
    rational_function_1d(nX, np, nq)
{
}

double rational_function_legendre_1d::legendre(double x, int i) const
{
	if(i == 0)
	{
		return 1;
	}
	else if(i == 1)
	{
		return x;
	}
	else
	{
		return ((2*i-1)*x*legendre(x, i-1) - (i-1)*legendre(x, i-2)) / (double)i ;
	}
}

// Get the p_i and q_j function
double rational_function_legendre_1d::p(const vec& x, int i) const
{
	std::vector<int> deg = index2degree(i);
	double res = 1.0;
	for(int k=0; k<_parameters.dimX(); ++k)
	{
		res *= legendre(2.0*((x[k] - _min[k]) / (_max[k]-_min[k]) - 0.5), deg[k]);
	}

	return res ;
}
double rational_function_legendre_1d::q(const vec& x, int i) const 
{
	return p(x, i);
}


rational_function_legendre::rational_function_legendre(const parameters& params)
    : rational_function(params)
{
}

#ifdef __GNUG__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

rational_function_legendre::rational_function_legendre()
{
}
#ifdef __GNUG__
#pragma GCC diagnostic pop
#endif
rational_function_legendre::~rational_function_legendre()
{
}

