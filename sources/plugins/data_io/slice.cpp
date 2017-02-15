/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2015, 2016 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include <core/data.h>

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "EXR_IO.h"

using namespace alta;


// Allow for a different parametrization depending on the arguments provided.
static const alta::parameters
brdf_slice_parameters(const arguments& args)
{
    auto result = alta::parameters(2, 3,
                                   params::STARK_2D, params::RGB_COLOR);

    if(args.is_defined("param")) {
				params::input param = params::parse_input(args["param"]);

				// The param is a 2D param
				if(params::dimension(param) == 2) {
            std::cout << "<<INFO>> Specified param \"" << args["param"] << "\"" << std::endl;
            result = alta::parameters(result.dimX(), result.dimY(),
                                      param, result.output_parametrization());

            // The oaram is a 3D param
				} else if(params::dimension(param) == 3) {
            std::cout << "<<INFO>> Specified param \"" << args["param"] << "\"" << std::endl;
            result = alta::parameters(3, result.dimY(),
                                      param, result.output_parametrization());

				} else {
            std::cout << "<<ERROR>> Invalid specified param \"" << args["param"] << "\"" << std::endl;
            std::cout << "<<ERROR>> Must have 2D input dimension" << std::endl;
				}
    }

    return result;
}

ALTA_DLL_EXPORT data* load_data(std::istream& input, const arguments& args);

/*! \ingroup datas
 *  \class data_brdf_slice
 *  \brief Data interface for the BRDF slice file format.
 *
 *  \details
 *
 *  This BRDF data format implements the 2D representation of a BRDF presented
 *  by Brent Burley in BRDF Explorer. It only stores a slice of the BRDF for
 *  different 2D parametrizations. This plugin stores the BRDF data into an
 *  EXR file.
 *
 *  It is possible to select the parametrization using the --param NAME
 *  argument when loading the BRDF. The default parametrization is
 *  \ref params::STARK_2D "STARK_2D"
 *
 *  \author Laurent Belcour <laurent.belcour@umontreal.ca>
 *
 */
class BrdfSlice : public data {
  private:
		const int _width, _height, _slice;
		vec _max, _min;
		double _phi;
		double* _data;
		bool _reverse;

  public:

		BrdfSlice(const arguments& args,
              int width, int height, int slice,
              double* content)
        : data(brdf_slice_parameters(args),
               width * height * slice),
          _width(width), _height(height), _slice(slice),
          _data(content)
		{
			// Allocate data
      if (args.is_defined("param") && parametrization().dimX() == 3)
          _phi = (M_PI / 180.0) * args.get_float("phi", 90);
      else
          _phi = 0.5*M_PI;

			// Is the position of the slice componnent (third coordinate)
			// reversed? This ensure that some params can be displayed.
      auto in_param = parametrization().input_parametrization();
			_reverse = in_param == params::ISOTROPIC_TL_TV_PROJ_DPHI ||
          in_param == params::SCHLICK_TL_TK_PROJ_DPHI   ||
          in_param == params::RETRO_TL_TVL_PROJ_DPHI;

			// Update the domain
			_max = max();
			_min = min();
		}

		BrdfSlice(const arguments& args)
        : BrdfSlice(args, 512, 512, 1,
                    new double[3 * 512 * 512])
    {
    }

		~BrdfSlice()
		{
			delete[] _data;
		}

		void save(const std::string& filename) const
		{
			if(!t_EXR_IO<double>::SaveEXR(filename.c_str(), _width, _slice*_height, _data))
			{
				std::cerr << "<<ERROR>> unable to save image file" << std::endl;
			}
		}

		// Acces to data
		vec get(int id) const
		{
      vec res(parametrization().dimX() + parametrization().dimY());
			const int i = id % _width;
			const int k = id / (_width*_height);
			const int j = (id - k*_width*_height) / _width;

			res[0] = (i+0.5) * (_max[0]-_min[0]) / double(_width)  + _min[0];
			res[1] = (j+0.5) * (_max[1]-_min[1]) / double(_height) + _min[1];
			if(parametrization().dimX() == 3) {
				res[2] = _phi;
			}

			// Reverse the first part of the vector
			if(_reverse) {
				res.segment(0, parametrization().dimX()).reverseInPlace();
			}

			res[parametrization().dimX()+0] = _data[3*id + 0];
			res[parametrization().dimX()+1] = _data[3*id + 1];
			res[parametrization().dimX()+2] = _data[3*id + 2];

			return res ;
		}

		void set(int id, const vec& x)
		{
			assert(x.size() == parametrization().dimX() + parametrization().dimY());

			_data[3*id + 0] = x[parametrization().dimX()+0];
			_data[3*id + 1] = x[parametrization().dimX()+1];
			_data[3*id + 2] = x[parametrization().dimX()+2];
		}

		vec value(const vec& x) const
		{
			// Copy vector is required
			vec _x = x;

			// Reverse the first part of the vector
			if(_reverse) {
				_x.segment(0, parametrization().dimX()).reverseInPlace();
			}

			// Safeguard. We can use either asserting or returning zero in case
			// the query values are not within reach.
			if(_x[0] > _max[0] || x[0] < _min[0] || x[1] > _max[1] || x[1] < _min[1] ||
				isnan(_x[0]) || isnan(x[1])) {
					vec res(3);
					return res;
			}
			/*
			assert(_x[0] <= 1.0 && x[0] >= 0.0);
			assert(_x[1] <= 1.0 && x[1] >= 0.0);
			*/
			const int i  = floor((_x[0]-_min[0]) * _width  / (_max[0] - _min[0]));
			const int j  = floor((_x[1]-_min[1]) * _height / (_max[1] - _min[1]));
			const int k  = 1;
			//const int k  = floor(_x[2] * _slice  / (M_PI));
			const int id = (i + j*_width)*k;

			if(i < 0 || i >= _width)  { std::cerr << "<<ERROR>> out of bounds: " << x << std::endl; }
			if(j < 0 || j >= _height) { std::cerr << "<<ERROR>> out of bounds: " << x << std::endl; }

			vec res(3);
			res[0] = _data[3*id + 0];
			res[1] = _data[3*id + 1];
			res[2] = _data[3*id + 2];
			return res;
		}

  private:

		// Get min and max input space values
		vec min() const
		{
			vec res(parametrization().dimX());

			// First fill the third dimension. It can be overwritten by the next
			// part when the parametrization is reversed (projected ones).
			if(parametrization().dimX() == 3) {
				res[2] = 0.0 ;
			}

			// Fill the first two dimension unless the parametrization is a
			// projected one then it will fill the three components.
			if(parametrization().input_parametrization() == params::ISOTROPIC_TL_TV_PROJ_DPHI ||
				parametrization().input_parametrization() == params::ISOTROPIC_TV_TL_DPHI) {
				res[0] = -0.5*M_PI ;
				res[1] = -0.5*M_PI ;
			} else if(parametrization().input_parametrization() == params::ISOTROPIC_TL_TV_PROJ_DPHI ||
					    parametrization().input_parametrization() == params::SCHLICK_TL_TK_PROJ_DPHI   ||
						 parametrization().input_parametrization() == params::RETRO_TL_TVL_PROJ_DPHI) {
				res[0] = -0.5*M_PI ;
				res[1] = -0.5*M_PI ;
				res[2] = -0.5*M_PI ;
			} else {
				res[0] = 0.0 ;
				res[1] = 0.0 ;
			}
			return res ;
		}
		vec max() const
		{
			vec res(parametrization().dimX());

			// First fill the third dimension. It can be overwritten by the next
			// part when the parametrization is reversed (projected ones).
			if(parametrization().dimX() == 3) {
				res[2] = 2.0*M_PI;
			}

			// Fill the first two dimension unless the parametrization is a
			// projected one then it will fill the three components.
			if(parametrization().input_parametrization() == params::RUSIN_TH_TD ||
				parametrization().input_parametrization() == params::RUSIN_TH_TD_PD ||
				parametrization().input_parametrization() == params::ISOTROPIC_TV_TL_DPHI) {
				res[0] = 0.5*M_PI ;
				res[1] = 0.5*M_PI ;
			} else if(parametrization().input_parametrization() == params::ISOTROPIC_TL_TV_PROJ_DPHI ||
					    parametrization().input_parametrization() == params::SCHLICK_TL_TK_PROJ_DPHI   ||
						 parametrization().input_parametrization() == params::RETRO_TL_TVL_PROJ_DPHI) {
				res[0] = 0.5*M_PI ;
				res[1] = 0.5*M_PI ;
				res[2] = 0.5*M_PI ;
			} else {
				res[0] = 1.0 ;
				res[1] = 1.0 ;
			}
			return res ;
		}

    friend data* load_data(std::istream&, const arguments&);
};

ALTA_DLL_EXPORT data* provide_data(size_t size, const parameters& params,
                                   const arguments& args)
{
    return new BrdfSlice(args);
}


ALTA_DLL_EXPORT data* load_data(std::istream& input, const arguments& args)
{
    int width, height, slice = 1;
    double *data;
    t_EXR_IO<double>::LoadEXR(input, width, height, data);

    BrdfSlice* result = new BrdfSlice(args, width, height, slice, data);

    return result;
}
