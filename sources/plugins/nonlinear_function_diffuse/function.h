/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2016 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#pragma once

// Include STL
#include <vector>
#include <string>

// Interface
#include <core/function.h>
#include <core/data.h>
#include <core/args.h>
#include <core/common.h>
#include <core/params.h>

using namespace alta;


/*!
 * \class nonlinear_function_diffuse
 * \ingroup functions
 * \brief A diffuse BRDF model.
 *
 * \details
 * This BRDF model has a constant value for any couple of light-view directions.
 * However, it is not directly used during the fitting procedure. The diffuse
 * function's values are set to be the minimum value in the dataset.
 *
 *  \author Laurent Belcour \<laurent.belcour@umontreal.ca\>
 */
class diffuse_function : public nonlinear_function
{

	public: // methods

    // Set the input parametrization to CARTESIAN to reduce the number
    // of transformations in a compound object.
    diffuse_function(const alta::parameters& params);

		// Overload the function operator
		virtual vec operator()(const vec& x) const ;
		virtual vec value(const vec& x) const ;


		//! \brief Boostrap the function by defining the diffuse term
		virtual void bootstrap(const ptr<data> d, const arguments& args);

		//! \brief Load function specific files
        virtual bool load(std::istream& in) ;
		virtual void save_call(std::ostream& out, const arguments& args) const;

		//! \brief Number of parameters to this non-linear function
		virtual int nbParameters() const ;

		//! \brief Get the vector of parameters for the function
		virtual vec parameters() const ;

		//! \brief Update the vector of parameters for the function
		virtual void setParameters(const vec& p) ;

		//! \brief Obtain the derivatives of the function with respect to the
		//! parameters.
		virtual vec parametersJacobian(const vec& x) const ;

	private: // data

		vec _kd;

    diffuse_function() {};
} ;

