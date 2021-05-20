/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2016 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include "rational_function.h"

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

using namespace alta;

rational_function_1d::rational_function_1d(const parameters& params,
                                           unsigned int np, unsigned int nq,
                                           bool separable)
    : function(params), _separable(separable)
{
    resize(np, nq);
}

rational_function_1d::rational_function_1d(int nX, unsigned int np, unsigned int nq, 
                                           bool separable):
    function(parameters(nX, 1, params::UNKNOWN_INPUT, params::UNKNOWN_OUTPUT))
{
	resize(np, nq);
	_separable = separable;
}

bool rational_function_1d::load(std::istream& in)
{
	// Variables
	int _np, _nq;
	std::string token;

	// Shoudl have the #NP [int]
	in >> token;
	if(token != "#NP") { return false; }
	in >> _np;

	// Shoudl have the #NQ [int]
	in >> token >> _nq;
	if(token != "#NQ") { return false; }
	vec a(_np), b(_nq);

	// Parse the p_i coefficients
	// TODO: Check if the indices match the current index?
	for(int j=0; j<_np; ++j)
	{
		for(int k=0; k<_parameters.dimX(); ++k)
		{
			in >> token;
		}
		in >> a[j];
	}

	// Parse the q_i coefficients
	// TODO: Check if the indices match the current index?
	for(int j=0; j<_nq; ++j)
	{
		for(int k=0; k<_parameters.dimX(); ++k)
		{
			in >> token;
		}
		in >> b[j];
	}
#ifdef NONE
	std::cout << "p_" << i << " = " << a << std::endl;
	std::cout << "q_" << i << " = " << b << std::endl;
#endif
	// Update the 1D function 
	//RP: WITHOUT RESIZE BEFORE THE .deg vector is never initialize correctly
	this->resize(_np,_nq);
	this->update(a, b);
	return true;
}

void rational_function_1d::save_body(std::ostream& out, const arguments& args) const
{
    bool const is_matlab = args["export"] == "matlab";
	bool const is_cpp      = args["export"] == "C++" || args["export"] == "c++" || args["export"] == "cpp";
    bool const is_alta   = !args.is_defined("export") || args["export"] == "alta";

	 const unsigned int np = _p_coeffs.size();
	 const unsigned int nq = _q_coeffs.size();

	 if(is_alta)
	 {
		 out << "#NP " << np << std::endl;
		 out << "#NQ " << nq << std::endl;
		 for(unsigned int i=0; i<np; ++i)
		 {
			 for(int j=0; j<_parameters.dimX(); ++j)
			 {
				 out << _p_coeffs[i].deg[j] << "\t" ;
			 }
			 out << _p_coeffs[i].a << std::endl ;
		 }

		 for(unsigned int i=0; i<nq; ++i)
		 {
			 for(int j=0; j<_parameters.dimX(); ++j)
			 {
				 out << _q_coeffs[i].deg[j] << "\t" ;
			 }
			 out << _q_coeffs[i].a << std::endl ;
		 }
	 }
	 else if(is_matlab)
	 {
		 out << "(";
		 for(unsigned int i=0; i<np; ++i)
		 {
			 out << _p_coeffs[i].a;
			 for(int k=0; k<_parameters.dimX(); ++k)
			 {
				 if(k != _parameters.dimX()-1) { out << ".*"; }
				 out << "x(k).^" << _p_coeffs[i].deg[k];
			 }
			 if(i != np-1) { out << " + "; }
		 }
		 out << ") / (";

		 for(unsigned int i=0; i<nq; ++i)
		 {
			 out << _q_coeffs[i].a << "x.^" << i;
			 for(int k=0; k<_parameters.dimX(); ++k)
			 {
				 if(k != _parameters.dimX()-1) { out << ".*"; }
				 out << "x(k).^" << _q_coeffs[i].deg[k];
			 }
			 if(i != nq-1) { out << " + "; }
		 }
		 out << ")";
	 }
	else if(is_cpp)
	 {
		 out << "\tdouble num = ";
		 for(unsigned int i=0; i<np; ++i)
		 {
			 out << _p_coeffs[i].a;
			 for(int k=0; k<_parameters.dimX(); ++k)
			 {
				 out << "*pow(x[" << k << "], " << _p_coeffs[i].deg[k] << ")";
			 }
			 if(i != np-1) { out << " + "; }
		 }
		 out << ";" << std::endl;

		out << "\tdouble den = ";
		 for(unsigned int i=0; i<nq; ++i)
		 {
			 out << _q_coeffs[i].a;
			 for(int k=0; k<_parameters.dimX(); ++k)
			 {
				 out << "*pow(x[" << k << "], " << _q_coeffs[i].deg[k] << ")";
			 }
			 if(i != nq-1) { out << " + "; }
		 }
		 out << ";" << std::endl;
		 out << "\tres[0] = num / den;" << std::endl;
	 }
	 else
	 {
		 NOT_IMPLEMENTED();
	 }
}

void rational_function_1d::update(const vec& in_a,
                                  const vec& in_b)
{
	// Get the size of the input vector
  const unsigned int np = in_a.size();
  const unsigned int nq = in_b.size();

	// Resize the coefficient vector if they do not match
	if(np != _p_coeffs.size()) {
		_p_coeffs.resize(np);
	}
	if(nq != _q_coeffs.size()) {
		_q_coeffs.resize(nq);
	}

#define NORMALIZE

#ifdef NORMALIZE
    const double b0 = (std::abs(in_b[0]) > 1.0E-16) ? in_b[0] : 1.0;
#else
    const double b0 = 1.0;
#endif

  for(unsigned int k=0; k<np; ++k)
	{
		_p_coeffs[k].a = in_a[k] / b0;
	}

  for(unsigned int k=0; k<nq; ++k)
	{
		_q_coeffs[k].a = in_b[k] / b0;
	}
}
		
void rational_function_1d::update(const rational_function_1d* r)
{
	// Get the size of the input vector
	const unsigned int np = r->_p_coeffs.size();
	const unsigned int nq = r->_q_coeffs.size();

	// Resize the rational function
	resize(np, nq);

	for(unsigned int k=0; k<np; ++k) {
		_p_coeffs[k].a = r->_p_coeffs[k].a;
	}

	for(unsigned int k=0; k<nq; ++k) {
		_q_coeffs[k].a = r->_q_coeffs[k].a;
	}

}

void rational_function_1d::resize(unsigned int np, unsigned int nq)
{
	// It is not possible to resize the multi-dimensional vector
	// if the input dimensions size is not defined. This can
	// happen at the creation of the rational function object.
	if(_parameters.dimX() == 0) { return; }

	// Resize the numerator
	if(_p_coeffs.size() != np)
	{
		_p_coeffs.resize(np);
		for(unsigned int k=0; k<np; ++k)
		{
			std::vector<int> deg = index2degree(k);
			_p_coeffs[k].deg = deg;
		}
	}

	// Resize the denominator
	if(_q_coeffs.size() != nq)
	{
		_q_coeffs.resize(nq);
		for(unsigned int k=0; k<nq; ++k)
		{
			std::vector<int> deg = index2degree(k);
			_q_coeffs[k].deg = deg;
		}
	}
}



// Get the p_i and q_j function
vec rational_function_1d::p(const vec& x) const
{
	vec res(1) ;

	const unsigned int np = _p_coeffs.size() ;
	for(unsigned int i=0; i<np; ++i)
	{
		res[0] += _p_coeffs[i].a * this->p(x, i) ;
	}

	return res ;
}
vec rational_function_1d::q(const vec& x) const 
{
	vec res(1) ;

	const unsigned int np = _q_coeffs.size() ;
	for(unsigned int i=0; i<np; ++i)
	{
		res[0] += _q_coeffs[i].a * this->q(x, i) ;
	}

	return res ;
}

// Estimate the number of configuration for an indice vector of dimension
// d with maximum element value being k. This function is heavy to call and
// should be use lightely. I did the maximum I could to unroll the dimension
// loop.
int rational_function_1d::estimate_dk(int k, int d)
{
	// For dim = 1, only one configuration is possible and only one 
	// configuration with zeroth powers remains
	if (d == 1 || k == 0) {
		return 1;

	// If only one power is available, there is dim configurations
	} else if (k == 1) {
		return d;

/*  // This one is not working. I need to redo the maths here.

	// If max power is two, then there is d^2 configuration minus the
	// number of symmetric configurations: d(d+1)/2.
	} else if (k == 2) {
		return (d*(d-1))/2;
*/

	// For dim = 2, configuration are [k-i, i], i \in [0,k]
	} else if (d == 2) {
		return k+1;

	// For dim = 3, we need to list all the possible cases for the first
	// dimension equals to i \in [0,k] and the remaining dimensions
	// sharing the k-i remaining elements.
	} else if (d == 3) {
		const int k2 = k+1;
		return k2*k2 - (k2*k)/2;

	} else {
		int res = 0;
		for(int i=0; i<=k; ++i)	{
			res += estimate_dk(k-i, d-1);
		}
		return res;
	}
}

// Populate a vector of degrees of dimension N using a
// maximum degree of M. The index at the current level
// is j
void rational_function_1d::populate(std::vector<int>& vec, int N, int M, int j)
{
	// For each dimension, estimate the current level
	// based on the number of configurations in the
	// other dimensions
	int current_M = M ;
	int nb_conf = 0;
	for(int d=0; d<N-1; ++d)
	{
		int k;
		for(k=0; k<=current_M; ++k)
		{
			int oracle = estimate_dk(current_M-k, N-(d+1));
			if(nb_conf <= j && j < nb_conf+oracle)
			{
				break;
			}
			nb_conf += oracle;
		}

		vec[N-1 - d] = k ;	
		current_M -= k ;
	}
	vec[0] = current_M;
}

std::vector<int> rational_function_1d::index2degree(int i) const
{
	std::vector<int> deg ; deg.assign(_parameters.dimX(), 0) ;

	if(i == 0)
		return deg ;
	
	// The case of one dimensional signals is trivial
	if(_parameters.dimX() == 1)
	{
		deg[0] = i;
		return deg;
	}
#ifdef NOT_WORKING
	// Case of two or more dimension signal, we differ the treatment of
	// separable signals and non separable signals.
	if(_separable)
	{
		const int index = i-1;
		const int base  = index / _parameters.dimX() + 1;
		for(int k=0; k<_parameters.dimX(); ++k)
		{
			deg[k] = (index % _parameters.dimX() == k) ? base : 0;
		}
	}
	else
#endif
	{
		if(_parameters.dimX() == 2)
		{
			int Nk = 1 ;
			int k  = 1 ;
			while(!(i >= Nk && i < Nk+k+1))
			{
				Nk += k+1 ;
				++k ;
			}

			int r = i-Nk ;
			deg[0] = k-r;
			deg[1] = r;
		}
		else
		{
			int Nk = 1 ;
			int k  = 1 ;
			int dk = estimate_dk(k, _parameters.dimX()) ;
			while(!(i >= Nk && i < Nk+dk))
			{
				Nk += dk ;
				++k ;
				dk = estimate_dk(k, _parameters.dimX()) ;
			}

			// Populate the vector from front to back
			int j = i-Nk ;
			populate(deg, _parameters.dimX(), k, j) ;
		}
	}

#ifdef CORE_DEBUG
	for(int k=0; k<deg.size(); ++k) {
		std::cout << deg[k] << ", ";
	}
	std::cout << std::endl;
#endif
	return deg ;

}

// Get the p_i and q_j function
double rational_function_1d::p(const vec& x, int i) const
{
	double res = 1.0;
	for(int k=0; k<_parameters.dimX(); ++k)
	{
		const double xp = 2.0*((x[k] - _min[k]) / (_max[k]-_min[k]) - 0.5);
		res *= pow(xp, _p_coeffs[i].deg[k]) ;
	}

	return res ;
}
double rational_function_1d::q(const vec& x, int i) const 
{
	double res = 1.0;
	for(int k=0; k<_parameters.dimX(); ++k)
	{
		const double xp = 2.0*((x[k] - _min[k]) / (_max[k]-_min[k]) - 0.5);
		res *= pow(xp, _q_coeffs[i].deg[k]) ;
	}

	return res ;
}

// Overload the function operator
vec rational_function_1d::value(const vec& x) const 
{
	vec res(1) ;

	unsigned int const np = _p_coeffs.size();
	unsigned int const nq = _q_coeffs.size();

	double p = 0.0f ;
	double q = 0.0f ;

	for(unsigned int i=0; i<np; ++i)
	{
		p += _p_coeffs[i].a*this->p(x, i) ;
	}

	for(unsigned int i=0; i<nq; ++i)
	{
		q += _q_coeffs[i].a*this->q(x, i) ;
	}

	res[0] = p/q ;
	return res ;
}


std::ostream& operator<< (std::ostream& out, const rational_function_1d& r) 
{
   const auto p = r.getP();
   const auto q = r.getQ();
   const unsigned int np = p.size();
   const unsigned int nq = q.size();

	std::cout << "p = [" ;
	for(unsigned int i=0; i<np; ++i)
	{
		if(i != 0)
		{
			std::cout << ", " ;
		}
     std::cout << p[i];
	}
	std::cout << "]" << std::endl ;

	std::cout << "q = [" ;
	for(unsigned int i=0; i<nq; ++i)
	{
		if(i != 0)
		{
			std::cout << ", " ;
		}
     std::cout << q[i];
	}
	std::cout << "]" << std::endl ;

	return out ;
}

#ifndef TODO

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

rational_function::rational_function() : np(0), nq(0)
{
}


rational_function::rational_function(int np, int nq) : np(np), nq(nq)
{
}

rational_function_1d::rational_function_1d()
{
}

#pragma GCC diagnostic pop

rational_function::rational_function(const parameters& params,
                                     int np, int nq):
    function(params), np(np), nq(nq)
{
    rs.resize(params.dimY());
}

//! \todo clean memory here
rational_function::~rational_function()
{
	for(unsigned int i=0; i < rs.size(); i++)
	{
		delete rs[i];
	}
}


void rational_function::update(int i, rational_function_1d* r)
{
	rs[i]->update(r);
}


rational_function_1d* rational_function::get(int i)
{
	// Check for consistency in the index of color channel
	if(i < _parameters.dimY())
	{
		if(rs[i] == NULL)
		{
			rs[i] = new rational_function_1d(_parameters.dimX(), np, nq);

			// Test if the input domain is not empty. If one dimension of
			// the input domain is a point, I manually inflate this dimension
			// to avoid numerical issues.
			vec _min = min();
			vec _max = max();
			for(int k=0; k<_parameters.dimX(); ++k)
			{
				if(_min[k] == _max[k]) 
				{
					_min[k] -= 0.1;
					_max[k] += 0.1;
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

rational_function_1d* rational_function::get(int i) const
{
	// Check for consistency in the index of color channel
	if(i < _parameters.dimY())
	{
		return rs[i];
	}
	else
	{
		std::cout << "<<ERROR>> tried to access out of bound 1D RF" << std::endl;
		return NULL;
	}
}

// Overload the function operator

vec rational_function::value(const vec& x) const
{
	vec res(_parameters.dimY()) ;

	for(int k=0; k<_parameters.dimY(); ++k)
	{
		res[k] = rs[k]->value(x)[0] ;
	}
	return res ;
}

// IO function to text files
bool rational_function::load(std::istream& in)
{

	// Parse line until the next comment
	while(in.peek() != '#')
	{
		char line[256];
		in.getline(line, 256);

		// If we cross the end of the file, or the badbit is
		// set, the file cannot be loaded
		if(!in.good())
			return false;
	}

    // Checking for the comment line #FUNC nonlinear_function_lafortune
	std::string token;
	in >> token;
	if(token.compare("#FUNC") != 0) 
	{ 
		std::cerr << "<<ERROR>> parsing the stream. The #FUNC is not the next line defined." << std::endl; 
		return false;
	}

	in >> token;
   if(token.compare("rational_function") != 0) 
	{
		std::cerr << "<<ERROR>> parsing the stream. function name is not the next token." << std::endl; 
		return false;
	}

	// Check for the MIN and MAX vector
	vec min(_parameters.dimX()), max(_parameters.dimX());
	in >> token;
   if(token.compare("#MIN") != 0) 
	{
		std::cerr << "<<ERROR>> the min value for the input space is not defined." << std::endl; 
		return false;
	}
	for(int k=0; k<_parameters.dimX(); ++k) {in >> min[k];}
	setMin(min);

	in >> token;
   if(token.compare("#MAX") != 0) 
	{
		std::cerr << "<<ERROR>> the max value for the input space is not defined." << std::endl; 
		return false;
	}
	for(int k=0; k<_parameters.dimX(); ++k) {in >> max[k]; }
	setMax(max);

	for(int i=0; i<_parameters.dimY(); ++i)
	{
		// Update the i_th color channel
		if(!get(i)->load(in)) { return false; }
	}

	return true;
}


void rational_function::save_call(std::ostream& out, const arguments& args) const
{
	// out.precision(8);
	out << std::scientific;

	bool is_alta = !args.is_defined("export") || args["export"] == "alta";
	if(is_alta) {
		out << "#FUNC rational_function" << std::endl;
		out << "#MIN "; for(int k=0; k<_parameters.dimX(); ++k) { out << _min[k] << " "; } out << std::endl;
		out << "#MAX "; for(int k=0; k<_parameters.dimX(); ++k) { out << _max[k] << " "; } out << std::endl; 
	}

	for(int k=0; k<_parameters.dimY(); ++k)
	{
      const rational_function_1d* rf = get(k);

      // We're calling the 'get(int) const' method, which can return NULL.
      // Skip them (XXX: is this the right thing to do?).
      if (rf != NULL)
      {
          rf->save_body(out, args);
          out << std::endl;
      }
	}
}

std::ostream& operator<< (std::ostream& out, const rational_function& r)
{
  for(int i=0; i < r.parametrization().dimY(); ++i)
	{
		rational_function_1d* rf = r.get(i);
		out << "dimension " << i << ": ";
		if(rf != NULL)
		{
			out << *rf << std::endl;
		}
		else
		{
			out << "[NULL]" << std::endl;
		}
	}

	return out ;
}


#endif
