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

using namespace alta;

ALTA_DLL_EXPORT function* provide_function(const alta::parameters& params)
{
    return new spherical_gaussian_function(params);
}

spherical_gaussian_function::spherical_gaussian_function(const alta::parameters& params)
    : nonlinear_function(params.set_input(6, params::CARTESIAN)),
      _a(1.0), _type(Mirror)
{
    // Update the length of the vectors
    _n.resize(_parameters.dimY()) ;
    _ks.resize(_parameters.dimY()) ;
}

// Overload the function operator
vec spherical_gaussian_function::operator()(const vec& x) const 
{
	return value(x);
}
vec spherical_gaussian_function::value(const vec& x) const 
{
    vec res(_parameters.dimY());
	double dot = compute_dot(x);

    for(int i=0; i<_parameters.dimY(); ++i)
	{
		res[i] = _ks[i] * std::exp(_n[i] * (dot-1));
	}

	return res;
}

double spherical_gaussian_function::compute_dot(const vec& x) const
{
	double dot;
	if(_type == Half)
	{
		params::convert(&x[0], params::CARTESIAN, params::COS_TH, &dot);
	}
	else if(_type == Back)
	{
		params::convert(&x[0], params::CARTESIAN, params::COS_TK, &dot);
	}
	else if(_type == Moment)
	{
		double spherical[4];
		params::convert(&x[0], params::CARTESIAN, params::SPHERICAL_TL_PL_TV_PV, &spherical[0]);
		spherical[1] *= _a;
		spherical[1] += M_PI;
		
		double cart[6];
		params::convert(&spherical[0], params::SPHERICAL_TL_PL_TV_PV, params::CARTESIAN, &cart[0]);
		dot = cart[0]*cart[3] + cart[1]*cart[4] + cart[2]*cart[5];
	}
	else
	{
		double spherical[4];
		params::convert(&x[0], params::CARTESIAN, params::SPHERICAL_TL_PL_TV_PV, &spherical[0]);
		spherical[1] += M_PI;
		
		double cart[6];
		params::convert(&spherical[0], params::SPHERICAL_TL_PL_TV_PV, params::CARTESIAN, &cart[0]);
		dot = cart[0]*cart[3] + cart[1]*cart[4] + cart[2]*cart[5];
	}

	return dot;
}

//! Number of parameters to this non-linear function
int spherical_gaussian_function::nbParameters() const 
{
	if(_type == Moment)
	{
        return 2*_parameters.dimY()+1;
	}
	else
	{
        return 2*_parameters.dimY();
	}
}

//! Get the vector of parameters for the function
vec spherical_gaussian_function::parameters() const 
{
    const int nY = _parameters.dimY();
	if(_type == Moment)
	{
        vec res(2*nY+1);
        res[2*nY] = _a;
        for(int i=0; i<nY; ++i)
		{
			res[i*2 + 0] = _ks[i];
			res[i*2 + 1] = _n[i];
		}

		return res;
	}
	else
	{
        vec res(2*nY);
        for(int i=0; i<nY; ++i)
		{
			res[i*2 + 0] = _ks[i];
			res[i*2 + 1] = _n[i];
		}
		return res;
	}
}

//! \brief get the min values for the parameters
vec spherical_gaussian_function::getParametersMin() const
{
    const int nY = _parameters.dimY();
	if(_type == Moment)
	{
        vec res(2*nY+1);
        res[2*nY] = 0.0;
        for(int i=0; i<nY; ++i)
		{
			res[i*2 + 0] = 0.0;
			res[i*2 + 1] = 0.0;
		}

		return res;
	}
	else
	{
        vec res(2*nY);
        for(int i=0; i<nY; ++i)
		{
			res[i*2 + 0] = 0.0;
			res[i*2 + 1] = 0.0;
		}

		return res;
	}
}

//! Update the vector of parameters for the function
void spherical_gaussian_function::setParameters(const vec& p) 
{
	if(_type == Moment)
	{
        _a = p[2*_parameters.dimY()];
	}

    for(int i=0; i<_parameters.dimY(); ++i)
	{
		_ks[i] = p[i*2 + 0];
		_n[i]  = p[i*2 + 1];
	}
}

//! Obtain the derivatives of the function with respect to the 
//! parameters. 
vec spherical_gaussian_function::parametersJacobian(const vec& x) const 
{
	double dot = compute_dot(x);

    vec jac(_parameters.dimY()*nbParameters());
     for(int i=0; i<_parameters.dimY(); ++i)
	 {
         for(int j=0; j<_parameters.dimY(); ++j)
		 {
			 if(i == j)
			 {

				 // df / dk_s
				 jac[i*nbParameters() + j*2+0] = exp(_n[j] * (dot - 1));

				 // df / dN
				 jac[i*nbParameters() + j*2+1] = _ks[j] * (dot-1) * exp(_n[j]*(dot - 1));

			 }
			 else
			 {
				 jac[i*nbParameters() + j*2+0] = 0.0;
				 jac[i*nbParameters() + j*2+1] = 0.0;
			 }
		 }

		 if(_type == Moment)
		 {
			 // df / da
			 double spherical[4];
			 params::convert(&x[0], params::CARTESIAN, params::SPHERICAL_TL_PL_TV_PV, &spherical[0]);
			 spherical[1] += 0.5* M_PI;

			 double cart[6];
			 params::convert(&spherical[0], params::SPHERICAL_TL_PL_TV_PV, params::CARTESIAN, &cart[0]);
			 double dot2 = cart[0]*cart[3] + cart[1]*cart[4] + cart[2]*cart[5];
			 jac[i*nbParameters() + nbParameters()-1] = _ks[i] * _n[i] * exp(_n[i] * (dot - 1)) * _a * dot2;

		 }
	 }

    return jac;
}
		
void spherical_gaussian_function::bootstrap(const ptr<data> d, const arguments& args)
{
    for(int i=0; i<_parameters.dimY(); ++i)
	{
		_ks[i] = 1.0;
		_n[i]  = 1.0;
	}

	// Parse the lobe type
	if(args.is_defined("sg-type"))
	{
		std::string stype = args.get_string("sg-type", "half");
		if(stype == "mirror")
		{
			_type = Mirror;
		}
		else if(stype == "half")
		{
			_type = Half;
		}
		else if(stype == "back")
		{
			_type = Back;
		}
		else if(stype == "retro")
		{
			_type = Retro;
		}
		else if(stype == "moment")
		{
			_type = Moment;
			_a = 1.0;
		}
	}
}

//! Load function specific files
bool spherical_gaussian_function::load(std::istream& in)
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
   if(token.compare("nonlinear_function_spherical_gaussian") != 0) 
	{
		std::cerr << "<<ERROR>> parsing the stream. function name is not the next token." << std::endl; 
        return false;
	}

	// Parse the lobe
    for(int i=0; i<_parameters.dimY(); ++i)
	{

		in >> token >> _ks[i];
		in >> token >> _n[i];
	}
	
	if(_type == Moment)
	{
		in >> token >> _a;
	}

	std::cout << "<<DEBUG>> Spherical gaussians plugin, found p= " << parameters() << std::endl;

    return true;
}


void spherical_gaussian_function::save_call(std::ostream& out, const arguments& args) const
{
    bool is_alta   = !args.is_defined("export") || args["export"] == "alta";

    const int nY = _parameters.dimY();
    if(is_alta)
    {
		out << "#FUNC nonlinear_function_spherical_gaussian" << std::endl ;

         for(int i=0; i<nY; ++i)
		 {
			 out << "Ks " << _ks[i] << std::endl;
			 out << "N  " <<  _n[i] << std::endl;
		 }
	
		 if(_type == Moment)
		 {
			 out << "a  " << _a;
		 }

		 out << std::endl;
	 }
	 else
	 {
		 out << "spherical_gaussian(L, V, N, X, Y, vec3(";
         for(int i=0; i<nY; ++i)
		 {
			 out << _ks[i];
             if(i < nY-1) { out << ", "; }
		 }

		 out << "), vec3(";
         for(int i=0; i<nY; ++i)
		 {
			 out << _n[i];
             if(i < nY-1) { out << ", "; }
		 }

		 out << "))";
	 }

}

void spherical_gaussian_function::save_body(std::ostream& out, const arguments& args) const
{
    bool is_shader = args["export"] == "shader" || args["export"] == "explorer";

    if(is_shader)
    {
        out << "vec3 spherical_gaussian(vec3 L, vec3 V, vec3 N, vec3 X, vec3 Y, vec3 ks, vec3 Nl)" << std::endl;
        out << "{" << std::endl;

		  switch(_type)
		  {
			  case Half:
				  out << "\tvec3 H = normalize(L + V);" << std::endl;
				  out << "\tvec3 ext_dot = vec3(dot(H,N));" << std::endl;
				  break;

			  case Moment:
				  out << "\tfloat t = " << _a << "*acos(dot(V, N));" << std::endl;
				  out << "\tfloat p = atan(dot(V, Y), dot(V, X)) + " << M_PI << ";" << std::endl;
				  out << "\tvec3 R = cos(p)*sin(t)*X + sin(p)*sin(t)*Y + cos(t)*N;" << std::endl;
				  out << "\tvec3 ext_dot = vec3(dot(L,R));" << std::endl;
				  break;

			  case Back:
				  out << "\tvec3 Vp = 2.0f*N*(dot(N,V)) - V;" << std::endl;
				  out << "\tvec3 K  = normalize(L + Vp);" << std::endl;
				  out << "\tvec3 ext_dot = vec3(dot(K,N));" << std::endl;
				  break;
			  
			  case Retro:
				  out << "\tvec3 ext_dot = vec3(dot(L,V));" << std::endl;
				  break;

			  default:
				  out << "\tvec3 R = 2*dot(V,N)*N - V;" << std::endl;
				  out << "\tvec3 ext_dot = vec3(dot(L,R));" << std::endl;
				  break;
		  }

        out << "\treturn ks * exp(Nl * (ext_dot - vec3(1)));" << std::endl;
        out << "}" << std::endl;
    }
}
