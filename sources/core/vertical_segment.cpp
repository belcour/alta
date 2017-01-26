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
#include <cstring>

using namespace alta;
using namespace Eigen;

//#define RELATIVE_ERROR

// Note: For some reason returning an rvalue here doesn't work, so let's hope
// RVO comes into play.
static vec data_min(Ref<MatrixXd> data)
{
    vec min = VectorXd::Constant(data.rows(),
                                 std::numeric_limits<double>::max());

    for(int i = 0; i < data.cols(); i++)
    {
        auto col = data.col(i);
        for (int k = 0; k < data.rows(); ++k)
        {
            min[k] = std::min(min[k], col[k]);
        }
    }

    return min;
}

static vec data_max(Ref<MatrixXd> data)
{
    vec max = VectorXd::Constant(data.rows(),
                                 -std::numeric_limits<double>::max());

    for(int i = 0; i < data.cols(); i++)
    {
        auto col = data.col(i);
        for (int k = 0; k < data.rows(); ++k)
        {
            max[k] = std::max(max[k], col[k]);
        }
    }

    return max;
}

// Return a matrix view of PTR showing only the 'dimX' first rows of that
// matrix.
static Ref<MatrixXd>
x_view(double *ptr, size_t cols, const parameters& params,
       vertical_segment::ci_kind kind)
{
    auto stride = params.dimX() + params.dimY()
        + vertical_segment::confidence_interval_columns(kind, params);

    return Map<MatrixXd, 0, OuterStride<> >(ptr, params.dimX(), cols,
                                            OuterStride<>(stride));
}

vertical_segment::vertical_segment(const parameters& params,
                                   size_t size,
                                   std::shared_ptr<double> input_data,
                                   ci_kind kind)
    : data(params, size,
           data_min(x_view(input_data.get(), size, params, kind)),
           data_max(x_view(input_data.get(), size, params, kind))),
      _data(input_data),
      _ci_kind(kind),
      _is_absolute(true), _dt(0.1)
{
}

// Work around the lack of array support in C++11's 'shared_ptr'.
static void delete_array(double *thing)
{
    delete[] thing;
}

vertical_segment::vertical_segment(const parameters& params, unsigned int rows):
    // Create a data object with a "value-initialized" (i.e., zeroed) array.
    vertical_segment(params, rows,
                     std::shared_ptr<double>
                     (new double[rows * (params.dimX() + 3 * params.dimY())]{},
                      delete_array))
{
}


void vertical_segment::get(int i, vec& x, vec& yl, vec& yu) const
{
    // Make sure we have the lower and upper bounds of Y.
    assert(confidence_interval_kind() == ASYMMETRICAL_CONFIDENCE_INTERVAL);

    auto row = matrix_view().row(i);

    x.resize(_parameters.dimX());
    x = row.segment(0, _parameters.dimX()).transpose();

    auto y = row.segment(_parameters.dimX(), _parameters.dimY()).transpose();
    yl.resize(_parameters.dimY()); yu.resize(_parameters.dimY());

    yl = row.segment(_parameters.dimX() + _parameters.dimY(),
                     _parameters.dimY()).transpose();

    yl += y;

    yu = row.segment(_parameters.dimX() + 2 * _parameters.dimY(),
                     _parameters.dimY()).transpose();
    yu += y;
}

void vertical_segment::get(int i, vec& yl, vec& yu) const
{
    vec x;
    get(i, x, yl, yu);
}

vec vertical_segment::get(int i) const
{
    return data_view().row(i);
}

void vertical_segment::set(int i, const vec& x)
{
   // Check if the input data 'x' has the size of a vertical segment (i.e. dimX+3*dimY),
   // if it only has the size of a single value, then create the segment.
   if(x.size() == column_number()
      || x.size() == _parameters.dimX() + _parameters.dimY())
   {
       auto row = matrix_view().row(i);
       row.head(x.size()) = x;
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
