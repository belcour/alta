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
#include <sstream>

// Interface
#include "function.h"
#include "data.h"
#include "fitter.h"
#include "args.h"
#include "common.h"


namespace alta {
   class rational_function_1d;
   class rational_function;
}
// STL stream ouput
std::ostream& operator<< (std::ostream&, const alta::rational_function_1d&);
std::ostream& operator<< (std::ostream&, const alta::rational_function&);


namespace alta {

/*! \brief A one dimensional rational function class. A rational function has
 *  the form \f$r(x) = \dfrac{\sum_{i} a_i p_i(x)}{b_i q_i(x)}$\f.
 */
class rational_function_1d : public function
{
	public: // methods

		rational_function_1d(const parameters& params,
                         unsigned int np = 0, unsigned int nq = 0,
                         bool separable = false);

      rational_function_1d() ALTA_DEPRECATED;
		rational_function_1d(int nX, unsigned int np, unsigned int nq, 
									bool separable = false) ALTA_DEPRECATED;
		virtual ~rational_function_1d() {}


		/* FUNCTION INHERITANCE */

		//! Overload the function operator
		virtual vec value(const vec& x) const ;
		virtual vec operator()(const vec& x) const { return value(x) ; }

		//! IO function to text files
		virtual bool load(std::istream& in);

		//! \brief Save the rational function expansion. It should
		//! not be store in any variable (e.g. "y = rf(x);") as the
		//! nD rational function can append factor to the 1D rational
		//! function.
		virtual void save_body(std::ostream&, const arguments&) const;


      /* RATIONAL FUNCTION SPECIFIC */

		//! Evaluate the numerator \f$p(\mathbf{x})\f$ of the rational
		//! function. This function is provided to allow fast 
		//! implementation. For example one can use the Clenshaw 
		//! algorithm to efficiently evaluate recursively defined
		//! polynomials. 
		virtual vec p(const vec& x) const ;

		//! Evaluate the denominator \f$q(\mathbf{x})\f$ of the rational
		//! function. This function is provided to allow fast 
		//! implementation. For example one can use the Clenshaw 
		//! algorithm to efficiently evaluate recursively defined
		//! polynomials. 
		virtual vec q(const vec& x) const ;

		//! Evaluate the basis function \f$p_i(\mathbf{x})\f$ for the
		//! numerator of the rational function.
		virtual double p(const vec& x, int i) const ;
		//! Evaluate the basis function \f$q_i(\mathbf{x})\f$ for the
		//! denominator of the rational function.
		virtual double q(const vec& x, int j) const ;


		//! Update the coefficient vectors with new values. The new values
		//! are normalized by the first element of the denominator 
		//! coefficients.
		virtual void update(const vec& in_a, 
		                    const vec& in_b) ;
		
		//! Update the 1D rational function with another one. Note that
		//! this function can change the dimensions of the coefficients
		//! vectors.
		virtual void update(const rational_function_1d* r) ;

		//! Resize the polynomial.
		virtual void resize(unsigned int np, unsigned int nq);


		//! Get the i-th coefficient of the numerator.
		virtual double getP(int i) const { return _p_coeffs[i].a; }

		//! Get the i-th coefficient of the denominator.
		virtual double getQ(int i) const { return _q_coeffs[i].a; }

		//! Get the vector of coefficient for the numerator.
		virtual vec getP() const 
		{
			const int np = _p_coeffs.size();
			vec t(np);
			for(int i=0; i<np; ++i) {t[i] = _p_coeffs[i].a; }
			return t; 
		}
		
		//! Get the vector of coefficient for the denominator.
		virtual vec getQ() const 
		{ 
			const int nq = _q_coeffs.size();
			vec t(nq);
			for(int i=0; i<nq; ++i) {t[i] = _q_coeffs[i].a; }
			return t; 
		}


		//! Convert a 1D index into a vector of degree for a
		//! multinomial coeffcient. The resulting vector v should
		//! be used as prod_k x[k]^v[k] for the monomial basis
		std::vector<int> index2degree(int i) const ;

	protected: // functions

		static int estimate_dk(int k, int d);

		static void populate(std::vector<int>& vec, int N, int M, int j);


	protected: // data

		// Structure to store a multi-dimensional coefficient. The
		// coeffcient, a, is associated with the vector of degree
		// for each dimension.
		struct coeff
		{
			coeff() {}
			coeff(double a, std::vector<int> deg) :
				a(a), deg(deg) { }

			double a;
			std::vector<int> deg;
		};

		// Table of coefficients and indices, sorted with respect
		// to the indices.
		std::vector<coeff> _p_coeffs;
		std::vector<coeff> _q_coeffs;

		//! Is the function separable with respect to its input dimensions?
		//! \todo Make possible to have only part of the dimensions
		//! separable.
		bool _separable;
} ;

/*! \brief Rational function define a BRDF model using a ratio of polynomials
 *  \ingroup core
 *
 *  \details
 *  A rational function define a BRDF model using a ratio of polynomials.
 */
class rational_function : public function
{

	public: // methods

		rational_function() ALTA_DEPRECATED;
		rational_function(int np, int nq) ALTA_DEPRECATED;

    rational_function(const parameters& params,
                      int np = 0, int nq = 0);

		virtual ~rational_function() ;

		// Overload the function operator
		virtual vec value(const vec& x) const ;
		virtual vec operator()(const vec& x) const { return value(x) ; }

		// IO function to text files
		virtual bool load(std::istream& in) ;

		// Update the function
		virtual void update(const ptr<rational_function>& r) 
		{
      assert(r->parametrization().dimX() == _parameters.dimX());
      assert(r->parametrization().dimY() == _parameters.dimY());

			for(int k=0; k < _parameters.dimY(); ++k)
			{
				get(k)->update(r->get(k));
			}
		}
		virtual void update(int i, rational_function_1d* r) ;

		//! Get the 1D function associated with color channel i. If no one exist, 
		//! this function allocates a new element. If i > nY, it returns NULL.
		virtual rational_function_1d* get(int i) ;
		virtual rational_function_1d* get(int i) const ALTA_DEPRECATED;

		//! Set the dimension of the output space of the function. This function 
		//! will update the size of the rs vector size.
		virtual void setDimY(int nY)
		{
      parameters new_params(_parameters.dimX(), nY,
                            _parameters.input_parametrization(),
                            _parameters.output_parametrization());
      _parameters = new_params;
      rs.resize(nY);
		}

		//! \brief Set the size of the rational function. Any newly created 1D 
		//! function will have np and nq fixed as defined.
		virtual void setSize(int np, int nq)
		{
			this->np = np;
			this->nq = nq;
			clear();
		}
		
		//! \brief Clear the vector of 1D rational functions.
		virtual void clear()
		{
			rs.clear();
			rs.resize(_parameters.dimY());
		}

		virtual void setMin(const vec& min)
		{
			function::setMin(min);
		}

		virtual void setMax(const vec& max)
		{
			function::setMax(max);
		}

		//! \brief Save the rational function to the rational format (see 
		//! \ref formating).
		virtual void save_call(std::ostream& out, const arguments& args) const ;

	protected: // data

		//! Store the y \in R rational functions. Each channel is a distinct 
		//! polynomial and should be fitted separately.
		std::vector<rational_function_1d*> rs ;

		//! Size of the polynomials
		//! \todo Change it by a more adaptive scheme, with different np, nq per 
		//! color channel?
		int np, nq;
};
}
