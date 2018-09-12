/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2015, 2016 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#pragma once

#include <string>
#include <utility>
#include <iostream>
#include <limits>
#include <fstream>
#include <cmath>
#include <cassert>

#include "common.h"
#include "args.h"
#include "params.h"
#include "clustering.h"
#include "ptr.h"


namespace alta {

/*! \brief A data object. Allows to load data from files.
 *  \ingroup core
 */
class data
{
  public: // methods

    data(const parameters &p, int size)
        : _parameters(p), _size(size) {}

    data(const parameters& p, int size, const vec& min, const vec& max)
        : _parameters(p), _size(size), _min(min), _max(max)
    {
        assert(min.size() == p.dimX());
        assert(max.size() == p.dimX());
    }

    /* TODO: Eventually mark the following constructors as deprecated.  */
#ifdef __GNUG__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
    data() {}

    data(unsigned int dim_X, unsigned int dim_Y)
      : _parameters(dim_X, dim_Y)
    {}
#ifdef __GNUG__
#pragma GCC diagnostic pop
#endif
    // Virtual destructor
    virtual ~data() {}

    // Save the data to a file
    virtual void save(const std::string& filename) const;

    // Acces to data
    virtual vec get(int i) const = 0 ;

    //! \brief Provide an evaluation of the data using interpolation. If
    //! the data object does not provide an interpolation mechanism, it
    //! should throw an exception.
    //!
    //! \details
    //! The input vector must have the parametrization of the data, and
    //! match the total dimension: dimX + dimY.
    virtual vec value(const vec& in) const = 0;

    //! \brief Put the sample inside the data at index I.
    virtual void set(int i, const vec& x) = 0;


    // Get data size, e.g. the number of samples to fit
    int size() const { return _size; };

    //! \brief Return true if this object is equal to DATA ±ε.
    virtual bool equals(const data& data,
                        double epsilon =
                        std::pow(1.0, -int(std::numeric_limits<double>::digits10 - 1)));


    const parameters& parametrization() const {
        return _parameters;
    }

    void setParametrization(const parameters& p) {
        _parameters = p;
    }


    /* Maximum values of the data */

    //! \brief Get the minimum value the input can take
    const vec& min() const { return _min; };

    //! \brief Get the maximum value the input can take
    const vec& max() const { return _max; };

  protected: // data

    parameters _parameters;
    int _size;
    vec _min, _max;
} ;

/*! \brief Change the parametrization of data to fit the parametrization of the
 *  function to be fitted.
 *
 *  \ingroup core
 *  \internal
 *  \todo Finish this class
 */
class data_params : public data
{
  public: // structures

    //! \brief when changing from a parametrization to another, you might
    //! lose some dimensions. This list enumerate the different operators
    //! that can be applied on the raw data to be clusterized.
    //! \note by default we use <em>none</em>, but if the input space
    //! dimension is reduced, the program will halt.
    enum clustering
    {
      MEAN,
      MEDIAN,
      NONE
    };

  public: // methods

    //! \brief contructor requires the definition of a base class that
    //! has a parametrization, and a new parametrization.
    data_params(const ptr<data> d, int size, params::input new_param,
                data_params::clustering method = data_params::NONE) :
      data(parameters(params::dimension(new_param),
                      d->parametrization().dimY(),
                      new_param,
                      d->parametrization().output_parametrization()),
           size),
      _clustering_method(method)
    {
      std::cout << "<<INFO>> Reparametrization of the data" << std::endl;
      //TODO
      //clustering<data>(d, _nY, d->parametrization(), new_param, _data);

      std::cout << "<<INFO>> clustering left " << _data.size() << "/" << d->size() << " elements" << std::endl;
      save(std::string("cluster.gnuplot"));
    }

    virtual vec value(const vec&) const
    {
      NOT_IMPLEMENTED();
    }

    // Acces to data
    virtual vec get(int i) const
    {
      return _data[i];
    }

    virtual void set(int i, const vec& x)
    {
      this->set(i, x);
    }

  protected: // data

    data_params::clustering _clustering_method;

    std::vector<vec> _data;
};
}

/* -*- c++ -*- */
