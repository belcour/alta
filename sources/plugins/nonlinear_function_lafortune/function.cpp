/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2016 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include "function.h"

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

#include <core/common.h>
#include <core/params.h>

using namespace alta;

ALTA_DLL_EXPORT function* provide_function(const alta::parameters& params)
{
    return new lafortune_function(params);
}

lafortune_function::lafortune_function(const alta::parameters& params) :
    nonlinear_function(params)
{
    auto nY = params.dimY();

    // Update the length of the vectors
    if(_isotropic)
        _C.resize(_n*nY*2) ;
    else
        _C.resize(_n*nY*3) ;
    _N.resize(_n*nY) ;
    _kd.resize(nY);

    for(int i=0; i<nY; ++i)
        _kd[i] = 0.0;
}

// Overload the function operator
vec lafortune_function::operator()(const vec& x) const 
{
	return value(x);
}
vec lafortune_function::value(const vec& x) const 
{
    const int nY = _parameters.dimY();
    vec res(nY);

#ifdef ADAPT_TO_PARAM
	vec y(6);
	params::convert(&x[0], _in_param, params::CARTESIAN, &y[0]);
#endif

	double dx, dy, dz;
#ifdef ADAPT_TO_PARAM
	dx = y[0]*y[3];
	dy = y[1]*y[4];
	dz = y[2]*y[5];
#else
	dx = x[0]*x[3];
	dy = x[1]*x[4];
	dz = x[2]*x[5];
#endif

	// For each color channel
    for(int i=0; i<nY; ++i)
	{
		// Start with the diffuse term
		res[i] = _kd[i];

		// For each lobe
		for(int n=0; n<_n; ++n)
		{
			double Cx, Cy, Cz, N;
			getCurrentLobe(n, i, Cx, Cy, Cz, N);

			const double d = Cx*dx + Cy*dy + Cz*dz;
			if(d > 0.0)
				res[i] += pow(d, N);
		}

	}

    return res;
}
        
vec lafortune_function::value(const vec& x, const vec& p)
{
	// Test input parameters for correct size and update the
	// values of the function.
	assert(p.size() == nbParameters());
	setParameters(p);

    const int nY = _parameters.dimY();
    vec res(nY);

#ifdef ADAPT_TO_PARAM
	vec y(6);
	params::convert(&x[0], _in_param, params::CARTESIAN, &y[0]);
#endif

	double dx, dy, dz;
#ifdef ADAPT_TO_PARAM
	dx = y[0]*y[3];
	dy = y[1]*y[4];
	dz = y[2]*y[5];
#else
	dx = x[0]*x[3];
	dy = x[1]*x[4];
	dz = x[2]*x[5];
#endif

	// For each lobe and for each color channel
    for(int i=0; i<nY; ++i)
	{
		// Start with the diffuse
		res[i] = _kd[i];

		// Add the n lobes
		for(int n=0; n<_n; ++n)
		{
			double N, Cx, Cy, Cz;
			getCurrentLobe(n, i, Cx, Cy, Cz, N);

			const double d = Cx*dx + Cy*dy + Cz*dz;
			if(d > 0.0)
				res[i] += pow(d, N);
		}
	}

	return res;
}

// Set the number of lobes of the Lafortune BRDF
void lafortune_function::setNbLobes(int N)
{
    _n = N;
    const int _nY = _parameters.dimY();

    // Update the length of the vectors
	 if(_isotropic)
		 _C.resize(_n*_nY*2) ;
	 else
		 _C.resize(_n*_nY*3) ;
    _N.resize(_n*_nY) ;
}

//! Number of parameters to this non-linear function
int lafortune_function::nbParameters() const 
{
    const int nY = _parameters.dimY();
#ifdef FIT_DIFFUSE
	if(_isotropic)
        return (3*_n+1)*nY;
	else
        return (4*_n+1)*nY;
#else
	if(_isotropic)
        return (3*_n)*nY;
	else
        return (4*_n)*nY;
#endif
}

//! Get the vector of parameters for the function
vec lafortune_function::parameters() const 
{
    const int nY = _parameters.dimY();
	vec res(nbParameters());
	for(int n=0; n<_n; ++n)
        for(int i=0; i<nY; ++i)
		{
			if(_isotropic)
			{
                res[(n*nY + i)*3 + 0] = _C[(n*nY + i)*2 + 0];
                res[(n*nY + i)*3 + 1] = _C[(n*nY + i)*2 + 1];
                res[(n*nY + i)*3 + 2] = _N[n*nY  + i];
			}
			else
			{
                res[(n*nY + i)*4 + 0] = _C[(n*nY + i)*3 + 0];
                res[(n*nY + i)*4 + 1] = _C[(n*nY + i)*3 + 1];
                res[(n*nY + i)*4 + 2] = _C[(n*nY + i)*3 + 2];
                res[(n*nY + i)*4 + 3] = _N[n*nY  + i];
			}
		}

#ifdef FIT_DIFFUSE
     for(int i=0; i<nY; ++i)
	 {
		 if(_isotropic)
		 {
            res[3*_n*nY + i] = _kd[i];
		 }
	 }
#endif
    return res;
}

//! Update the vector of parameters for the function
void lafortune_function::setParameters(const vec& p) 
{
    const int nY = _parameters.dimY();
	// Safety check the number of parameters
	assert(p.size() == nbParameters());

	for(int n=0; n<_n; ++n)
        for(int i=0; i<nY; ++i)
		{
            _C[(n*nY + i)*3 + 0] = p[(n*nY + i)*4 + 0];
            _C[(n*nY + i)*3 + 1] = p[(n*nY + i)*4 + 1];
            _C[(n*nY + i)*3 + 2] = p[(n*nY + i)*4 + 2];
            _N[n*nY  + i]        = p[(n*nY + i)*4 + 3];
		}
#ifdef FIT_DIFFUSE
    for(int i=0; i<nY; ++i)
	{
        _kd[i] = p[4*_n*nY + i];
	}
#endif
}

//! Obtain the derivatives of the function with respect to the 
//! parameters. 
vec lafortune_function::parametersJacobian(const vec& x) const 
{
    const int nY = _parameters.dimY();

#ifdef ADAPT_TO_PARAM
	vec y(6);
	params::convert(&x[0], _in_param, params::CARTESIAN, &y[0]);
#endif

	double dx, dy, dz;
#ifdef ADAPT_TO_PARAM
	dx = y[0]*y[3];
	dy = y[1]*y[4];
	dz = y[2]*y[5];
#else
	dx = x[0]*x[3];
	dy = x[1]*x[4];
	dz = x[2]*x[5];
#endif

    vec jac(nY*nbParameters());
     for(int i=0; i<nY; ++i)
	 {
		 for(int n=0; n<_n; ++n)
             for(int j=0; j<nY; ++j)
			 {
				 // index of the current monochromatic lobe
                 int index = i*nbParameters() + 4*(n*nY + j);
				 
				 double Cx, Cy, Cz, N;
				 getCurrentLobe(n, j, Cx, Cy, Cz, N);
				 
				 double d  = Cx*dx + Cy*dy + Cz*dz;

				 if(i == j && d > 0.0)
				 {
					 // df / dCx
					 jac[index+0] = dx * N * std::pow(d, N-1.0);

					 // df / dCy
					 jac[index+1] = dy * N * std::pow(d, N-1.0);
					 
					 // df / dCz
					 jac[index+2] = dz * N * std::pow(d, N-1.0);

					 // df / dN
					 if(d <= 0.0)
						 jac[index+3] = 0.0;
					 else
						 jac[index+3] = std::log(d) * std::pow(d, N);
				 }
				 else
				 {
					 jac[index+0] = 0.0;
					 jac[index+1] = 0.0;
					 jac[index+2] = 0.0;
					 jac[index+3] = 0.0;
				 }
			 }

#ifdef FIT_DIFFUSE
         for(int j=0; j<nY; ++j)
		 {
			 // index of the current monochromatic lobe
             int index = i*nbParameters() + 4*_n*nY + j;

			 jac[index] = 1.0;
		 }
#endif
	 }

    return jac;
}
		
void lafortune_function::bootstrap(const ptr<data> d, const arguments& args)
{
    const int nY = _parameters.dimY();

    // Check the arguments for the number of lobes
    this->setNbLobes(args.get_int("lobes", 1));

    // Set the diffuse component
	vec x0 = d->get(0);
    for(int i=0; i<d->parametrization().dimY(); ++i)
        _kd[i] = x0[d->parametrization().dimX() + i];

	for(int i=1; i<d->size(); ++i)
	{
		vec xi = d->get(i);
        for(int j=0; j<d->parametrization().dimY(); ++j)
            _kd[j] = std::min(xi[d->parametrization().dimX() + j], _kd[j]);
	}

    // The remaining data will be equal to one
    for(int n=0; n<_n; ++n)
        for(int i=0; i<_parameters.dimY(); ++i)
        {
            double theta = 0.5 * M_PI * n / (double)_n;

            _C[(n*nY + i)*3 + 0] = -sin(theta);
            _C[(n*nY + i)*3 + 1] = -sin(theta);
            _C[(n*nY + i)*3 + 2] = cos(theta);
            _N[n*nY  + i]        = (double)_n;
        }
}

std::ofstream& type_definition(std::ofstream& out, int nY)
{
    if(nY == 1)
        out << "float " ;
    else
        out << "vec" << nY ;

    return out;
}

std::ofstream& type_affectation(std::ofstream& out, const std::string& name, const vec& x, int  nY, int n=0, int s=0, int S=1)
{

    out << name << " = ";

	if(nY != 1)
		out << "vec" << nY << "(";

	for(int i=0; i<nY; ++i)
	{
		if(i != 0) out << ", ";
        out << x[n*nY*S + i*S+s];
	}

	if(nY != 1)
		out << ")";

	out << ";" << std::endl;

	return out;
}


//! Load function specific files
bool lafortune_function::load(std::istream& in)
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
    if(token != "#FUNC")
    {
        std::cerr << "<<ERROR>> parsing the stream. The #FUNC is not the next line defined." << std::endl;
        return false;
    }

    in >> token;
    if(token != "nonlinear_function_lafortune")
    {
        std::cerr << "<<ERROR>> parsing the stream. function name is not the next token." << std::endl;
        return false;
    }

	 // Shoudl have the #NB_LOBES [int]
	 int nb_lobes;
	 in >> token >> nb_lobes;
	 setNbLobes(nb_lobes);

	// Parse the lobe
    const int _nY = _parameters.dimY();
	for(int n=0; n<_n; ++n)
	{
		for(int i=0; i<_nY; ++i)
		{

			in >> token >> _C[(n*_nY + i)*3 + 0];
            in >> token >> _C[(n*_nY + i)*3 + 1];
			in >> token >> _C[(n*_nY + i)*3 + 2];
			in >> token >> _N[i];
		}

	}
	
	std::cout << "<<INFO>> Cd = " << _C << std::endl;
	std::cout << "<<INFO>> N = " << _N << std::endl;
    return true;
}

void lafortune_function::save_call(std::ostream& out, const arguments& args) const
{
    bool is_alta   = !args.is_defined("export") || args["export"] == "alta";

    const int _nY = _parameters.dimY();
    if(is_alta)
    {
        out << "#FUNC nonlinear_function_lafortune" << std::endl ;
        out << "#NB_LOBES " << _n << std::endl ;

        for(int n=0; n<_n; ++n)
        {

            for(int i=0; i<_nY; ++i)
            {
                out << "Cx " << _C[(n*_nY + i)*3 + 0] << std::endl;
                out << "Cx " << _C[(n*_nY + i)*3 + 1] << std::endl;
                out << "Cz " << _C[(n*_nY + i)*3 + 2] << std::endl;
                out << "N  " << _N[n*_nY + i] << std::endl;

            }


            out << std::endl;
        }
    }
    else
    {
        for(int n=0; n<_n; ++n)
        {
            out << "lafortune(L, V, N, X, Y, vec3(";
            for(int i=0; i<_nY; ++i)
            {
                out << _C[(n*_nY + i)*3 + 0];
                if(i < _nY-1) { out << ", "; }
            }

            out << "), vec3(";
            for(int i=0; i<_nY; ++i)
            {
                out << _C[(n*_nY + i)*3 + 1];
                if(i < _nY-1) { out << ", "; }
            }

            out << "), vec3(";
            for(int i=0; i<_nY; ++i)
            {
                out << _C[(n*_nY + i)*3 + 2];
                if(i < _nY-1) { out << ", "; }
            }

            out << "), vec3(";
            for(int i=0; i<_nY; ++i)
            {
                out << _N[n*_nY + i];
                if(i < _nY-1) { out << ", "; }
            }

            // For multiple lobes, add a sum sign
            out << "))";
            if(n < _n-1) { out << " + "; }
        }
    }
}

void lafortune_function::save_body(std::ostream& out, const arguments& args) const
{
    out << "vec3 lafortune(vec3 L, vec3 V, vec3 N, vec3 X, vec3 Y, vec3 Cx, vec3 Cy, vec3 Cz, vec3 Nl)" << std::endl;
    out << "{" << std::endl;
    out << "\tvec3 ext_dot = Cx * dot(L,X)*dot(V,X) + Cy * dot(L,Y)*dot(V,Y) + Cz * dot(L,N)*dot(V,N);" << std::endl;
    out << "\treturn pow(max(ext_dot, vec3(0,0,0)), Nl);" << std::endl;
    out << "}" << std::endl;
    out << std::endl;
}
