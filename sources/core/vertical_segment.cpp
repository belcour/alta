/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2015, 2016 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include "vertical_segment.h"
#include "data_storage.h"

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cassert>

using namespace alta;

//#define RELATIVE_ERROR

// Note: For some reason returning an rvalue here doesn't work, so let's hope
// RVO comes into play.
static vec data_min(unsigned int size, const std::vector<vec>& data)
{
    vec min(size);

    for(int k = 0; k < size; ++k)
    {
        min[k] = std::numeric_limits<double>::max();
    }

    for (auto&& item: data)
    {
        for (int k = 0; k < size; ++k)
        {
            min[k] = std::min(min[k], item[k]);
        }
    }

    return min;
}

static vec data_max(unsigned int size, const std::vector<vec>& data)
{
    vec max(size);

    for(int k = 0; k < size; ++k)
    {
        max[k] = -std::numeric_limits<double>::max();
    }

    for (auto&& item: data)
    {
        for (int k = 0; k < size; ++k)
        {
            max[k] = std::max(max[k], item[k]);
        }
    }

    return max;
}


vertical_segment::vertical_segment(const parameters& params,
                                   std::vector<vec>&& input_data)
    : data(params, input_data.size(),
           data_min(params.dimX(), input_data),
           data_max(params.dimX(), input_data)),
      _data(input_data), _is_absolute(true), _dt(0.1)
{
}

vertical_segment::vertical_segment(const parameters& params, unsigned int size):
    data(params, size), _is_absolute(true), _dt(0.1)
{
    initializeToZero(size);
}


void
vertical_segment::initializeToZero( unsigned int number_of_data_elements )
{
	_data.clear();

  unsigned int const  size_of_one_element = _parameters.dimX() + 3*_parameters.dimY();
	
	for( unsigned int i=0; i < number_of_data_elements; i++)
	{
		vec initial_element = vec::Zero( size_of_one_element );
		_data.push_back( initial_element );
	}
}

void vertical_segment::get(int i, vec& x, vec& yl, vec& yu) const
{
#ifdef DEBUG
    assert(i >= 0 && i < _data.size());
#endif
    x.resize(_parameters.dimX()); yl.resize(_parameters.dimY()) ; yu.resize(_parameters.dimY()) ;
    for(int j=0; j<_parameters.dimX(); ++j)
    {
        x[j] = _data[i][j];
    }
    for(int j=0; j<_parameters.dimY(); ++j)
    {
        yl[j] = _data[i][_parameters.dimX() + 1*_parameters.dimY() + j];
        yu[j] = _data[i][_parameters.dimX() + 2*_parameters.dimY() + j];
    }
}
		
void vertical_segment::get(int i, vec& yl, vec& yu) const
{
	yl.resize(_parameters.dimY()) ; yu.resize(_parameters.dimY()) ;
	for(int j=0; j<_parameters.dimY(); ++j)
	{
		yl[j] = _data[i][_parameters.dimX() + _parameters.dimY() + j] ;
		yu[j] = _data[i][_parameters.dimX() + 2*_parameters.dimY() + j] ;
	}
}

vec vertical_segment::get(int i) const 
{
	//SLOW !!!!! and useless 
	// call _data[i]
    const int n = _parameters.dimX() + _parameters.dimY();
    vec res = _data[i].head(n);
    
    return res ;
}

void vertical_segment::set(int i, const vec& x)
{
   // Check if the input data 'x' has the size of a vertical segment (i.e. dimX+3*dimY),
   // if it only has the size of a single value, then create the segment.
   if(x.size() == _parameters.dimX() + 3*_parameters.dimY()) {
      _data[i] = x;

   } else if(x.size() == _parameters.dimX() + _parameters.dimY()) {
      _data[i] = vs(x);

   } else {
      std::cerr << "<<ERROR>> Passing an incorrect element to vertical_segment::set" << std::endl;
      throw;
   }
}

vec vertical_segment::vs(const vec& x) const {
   vec y(_parameters.dimX() + 3*_parameters.dimY());

   // Copy the head of each vector
   y.head(_parameters.dimX() + _parameters.dimY()) = x.head(_parameters.dimX() + _parameters.dimY());
   
   for(unsigned int i=0; i<_parameters.dimY(); ++i) {
      const double val = x[i + _parameters.dimX()];
      if(_is_absolute) {
         y[i + _parameters.dimX()+1*_parameters.dimY()] = val - _dt;
         y[i + _parameters.dimX()+2*_parameters.dimY()] = val + _dt;
      } else {
         y[i + _parameters.dimX()+1*_parameters.dimY()] = val * (1.0 - _dt);
         y[i + _parameters.dimX()+2*_parameters.dimY()] = val * (1.0 + _dt);
      }
   }

   return y;
}
