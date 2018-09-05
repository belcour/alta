/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2016 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include "rational_function.h"

#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

#include <core/common.h>

using namespace alta;

ALTA_DLL_EXPORT function* provide_function(const alta::parameters& params)
{
    return new rational_function_chebychev(params);
}


rational_function_chebychev::rational_function_chebychev(const alta::parameters& params,
                                                         int np, int nq)
    : rational_function(params, np, nq)
{
}

#ifdef __GNUG__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

rational_function_chebychev::rational_function_chebychev() : rational_function() 
{
}

rational_function_chebychev::~rational_function_chebychev()
{
}


rational_function_chebychev_1d::rational_function_chebychev_1d()
{
}

rational_function_chebychev_1d::rational_function_chebychev_1d(int nX, int np, int nq) :
    rational_function_1d(nX, np, nq)
{
}

#ifdef __GNUG__
#pragma GCC diagnostic pop
#endif

rational_function_chebychev_1d::rational_function_chebychev_1d(const parameters& params,
                                                               int np, int nq)
    : rational_function_1d(params, np, nq)
{
}

double recursive_chebychev(double x, int i) {
   switch(i) {
      case 0:
         return 1.0;
      case 1:
         return x;
      case 2:
         return 2*x*x - 1;
      case 3:
         return (4*x*x - 3)*x;
      case 4:
         return 8.0*(x*x - 1.0)*x*x + 1.0;
      default:
         return 2.0*x*recursive_chebychev(x, i-1) - recursive_chebychev(x, i-2);
   }
}

double chebychev(double x, int i) {
#ifdef RECURSIVE_FORM
   recursive_chebychev(x, i);
#else
   if(std::abs(x) <= 1.0) {
      return cos(i * acos(x));
   } else {
      return recursive_chebychev(x, i);
   }
#endif
}

// Get the p_i and q_j function
double rational_function_chebychev_1d::p(const vec& x, int i) const
{
	std::vector<int> deg = index2degree(i);
	double res = 1.0;
	for(int k=0; k<_parameters.dimX(); ++k)
	{
		double xk = 2.0*((x[k] - _min[k]) / (_max[k]-_min[k]) - 0.5);
		res *= chebychev(xk, deg[k]);
	}

	return res ;
}
double rational_function_chebychev_1d::q(const vec& x, int i) const
{
	std::vector<int> deg = index2degree(i);
	double res = 1.0; 
	for(int k=0; k<_parameters.dimX(); ++k)
	{
		double xk = 2.0*((x[k] - _min[k]) / (_max[k]-_min[k]) - 0.5);
		res *= chebychev(xk, deg[k]);
	}

	return res ;
}
		
rational_function_1d* rational_function_chebychev::get(int i)
{
	// Check for consistency in the index of color channel
	if(i < _parameters.dimY())
	{
		if(rs[i] == NULL)
		{
			rs[i] = new rational_function_chebychev_1d(_parameters.dimX(), np, nq);

			// Test if the input domain is not empty. If one dimension of
			// the input domain is a point, I manually inflate this dimension
			// to avoid numerical issues.
			vec _min = min();
			vec _max = max();
			for(int k=0; k<_parameters.dimX(); ++k)
			{
				if(_min[k] == _max[k]) 
				{
					_min[k] -= 1.0;
					_max[k] += 1.0;
				}
			}

			rs[i]->setMin(_min) ;
			rs[i]->setMax(_max) ;
		}
		return rs[i];
	}
	else
	{
		std::cout << "<<ERROR>> tried to access out of bound 1D RF" << std::endl;
		return NULL;
	}
}
