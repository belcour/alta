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
#include <core/fitter.h>
#include <core/args.h>
#include <core/common.h>
#include <core/params.h>

using namespace alta;

/*! 
 * \class nonlinear_function_sg
 * \ingroup plugins
 * \ingroup functions
 * \brief A spherical Gaussian lobe class.
 *
 *  \details
 *  A spherical_gaussian lobe is defined as \f$k_s exp(n (v.p - 
 *  1))\f$. Where \f$n\f$ is the roughness control and \f$p\f$
 *  is the lobe direction and \f$v\f$ the evaluation direction.
 *
 *
 *  <h3>Plugin parameters</h3>
 *
 *  <ul>
 *     <li><b>bootstrap</b> function will read the arguments to set 
 *     lobe and evaluation directions.</li>
 *  </ul>
 */
class spherical_gaussian_function : public nonlinear_function
{

	public: // methods

		enum type
		{
			Mirror,
			Half,
			Retro,
			Back,
			Moment
		};

    spherical_gaussian_function(const alta::parameters& params);

		// Overload the function operator
		virtual vec operator()(const vec& x) const ;
		virtual vec value(const vec& x) const ;

		//! \brief Load function specific files
        virtual bool load(std::istream& in) ;

		//! \brief Export function
		virtual void save_call(std::ostream& out, const arguments& args) const;
      virtual void save_body(std::ostream& out, const arguments& args) const;

		//! \brief Boostrap the function by defining the diffuse term
		//!
		//! \details
		virtual void bootstrap(const ptr<data> d, const arguments& args);

		//! \brief Number of parameters to this non-linear function
		virtual int nbParameters() const ;

		//! \brief Get the vector of parameters for the function
		virtual vec parameters() const ;

		//! \brief get the min values for the parameters
		virtual vec getParametersMin() const;

		//! \brief Update the vector of parameters for the function
		virtual void setParameters(const vec& p) ;

		//! \brief Obtain the derivatives of the function with respect to the
		//! parameters. 
		virtual vec parametersJacobian(const vec& x) const ;

	private: // methods

		//! \brief Compute the cosine for inside the lobe function. 
		//! Depending on the lobe type, the dot product can have 
		//! different evaluations.
		double compute_dot(const vec& in) const;

	private: // data

		vec _n, _ks; // Lobes data
		double _a;   // Scaling factor for the moment based SG

		type _type;  // Lobe type

    spherical_gaussian_function();
} ;

