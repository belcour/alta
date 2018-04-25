/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2016 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include "function.h"

#include <string>
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
	return new schlick_masking(params);
}

schlick_masking::schlick_masking(const alta::parameters& params)
    : nonlinear_function(params.set_input(6, params::CARTESIAN))
{
    w.resize(params.dimY());
}

//! Load function specific files
bool schlick_masking::load(std::istream& in)
{
//    fresnel::load(in);

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

    // Checking for the comment line #FUNC nonlinear_fresnel_schlick_masking
    std::string token;
    in >> token;
    if(token != "#FUNC")
    {
        std::cerr << "<<ERROR>> parsing the stream. The #FUNC is not the next line defined." << std::endl;
        return false;
    }

    in >> token;
    if(token != "nonlinear_masking_schlick")
    {
        std::cerr << "<<ERROR>> parsing the stream. function name is not the next token." << std::endl;
        return false;
    }

    // R [double]
    for(int i=0; i<_parameters.dimY(); ++i)
    {
        in >> token >> w[i];
    }
    return true;
}

void schlick_masking::save_call(std::ostream& out, const arguments& args) const
{
    bool is_alta   = !args.is_defined("export") || args["export"] == "alta";

    if(is_alta)
    {
        out << "#FUNC nonlinear_masking_schlick" << std::endl ;
        for(int i=0; i<_parameters.dimY(); ++i)
        {
            out << "K " << w[i] << std::endl;
        }
        out << std::endl;
    }
    else
    {
        out << "masking_schlick(L, V, N, X, Y, vec3";
        for(int i=0; i<_parameters.dimY(); ++i)
        {
            out << w[i];
            if(i < _parameters.dimY()-1) { out << ", "; }
        }
        out << "))";
    }
}

void schlick_masking::save_body(std::ostream& out, const arguments& args) const
{
    bool is_shader = args["export"] == "shader" || args["export"] == "explorer";

    if(is_shader)
    {
        out << std::endl;
        out << "vec3 masking_schlick(vec3 L, vec3 V, vec3 N, vec3 X, vec3 Y, vec3 K)" << std::endl;
        out << "{" << std::endl;
        out << "\tconst float dotLN = dot(L, N);" << std::endl;
        out << "\tconst float dotVN = dot(L, N);" << std::endl;
		  out << std::endl;
        out << "\tconst vec3 GL = dotLN / (dotLN + K * (dotLN - 1.0));" << std::endl;
        out << "\tconst vec3 GV = dotVN / (dotVN + K * (dotVN - 1.0));" << std::endl;
		  out << std::endl;
        out << "\treturn vec3(GL*GV);" << std::endl;
        out << "}" << std::endl;
        out << std::endl;
    }

}

vec schlick_masking::value(const vec& x) const
{
    vec res(_parameters.dimY());
	const double u = x[5];
	const double v = x[2];

    for(int i=0; i<_parameters.dimY(); ++i)
	{
		const double Gu = u / (u + w[i] * (1.0 - u));
		const double Gv = v / (v + w[i] * (1.0 - v));
		res[i] = Gu*Gv;
	}
	return res;
}

//! \brief Number of parameters to this non-linear function
int schlick_masking::nbParameters() const
{
    return _parameters.dimY();
}

vec schlick_masking::getParametersMin() const
{
    vec m(_parameters.dimY());
    for(int i=0; i<_parameters.dimY(); ++i) { m[i] = 0.0; }
    return m;
}

//! \brief Get the vector of parameters for the function
vec schlick_masking::parameters() const
{
    vec p(_parameters.dimY());
    for(int i=0; i<_parameters.dimY(); ++i) { p[i] = w[i]; }
    return p;
}

//! \brief Update the vector of parameters for the function
void schlick_masking::setParameters(const vec& p)
{
    for(int i=0; i<_parameters.dimY(); ++i) { w[i] = p[i]; }
}

//! \brief Obtain the derivatives of the function with respect to the
//! parameters.
vec schlick_masking::parametersJacobian(const vec& x) const
{
    const int nY = _parameters.dimY();
	vec jac(nY*nY);

	const double u = x[5];
	const double v = x[2];

	for(int i=0; i<nY; ++i)
		for(int j=0; j<nY; ++j)
		{
			if(i == j) 
			{
//				const double denom_u
				const double Gu = u / (u + w[i] * (1.0 - u));
				const double Gv = v / (v + w[i] * (1.0 - v));

				const double dGu = - u*(1.0 - u) / pow(u + w[i]*(1.0-u), 2);
				const double dGv = - v*(1.0 - v) / pow(v + w[i]*(1.0-v), 2);

                jac[j*_parameters.dimY() + i] = Gu*dGv + Gv*dGu;
			}
			else 
			{
                jac[j*_parameters.dimY() + i] = 0.0;
			}
		}

	return jac;
}


void schlick_masking::bootstrap(const ptr<data>, const arguments&)
{
	// Start with a non occluding value for k
    for(int i=0; i<_parameters.dimY(); ++i) { w[i] = 0.0; }
}
