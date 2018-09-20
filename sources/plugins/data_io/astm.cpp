/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2015, 2016 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include <core/data.h>
#include <core/vertical_segment.h>
#include <core/common.h>
#include <core/args.h>

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <utility>

namespace alta {

/*! \ingroup datas
 *  \ingroup plugins
 *  \class data_astm
 *  \brief Data interface for the ASTM file format from [Cornell University][cornell].
 *  [cornell]: http://www.graphics.cornell.edu/online/measurements/reflectance/
 *
 *  \details
 *  This plugin enables to load data measurments from Cornell University
 *  program of Computer Graphics (see [here][cornell]). The file format currently
 *  only handles integrated luminance [300, 700] and R,G,B color space due to
 *  ALTA's internal format. The X,Y,Z color space files should not be used with
 *  this plugin.
 *
 *  The resulting data object is a \a vertical_segment and is compatible with all
 *  rational fitting plugins. However, the dimension of the input space is usually
 *  not tight (i.e. larger than the dimensionality of the data). Thus it is advised
 *  to convert those data file before hand.
 *
 *  [cornell]: http://www.graphics.cornell.edu/online/measurements/reflectance/
 *
 *  \author Laurent Belcour <laurent.belcour@umontreal.ca>
 *
 */
class ASTM : public vertical_segment
{

public: //methods
    ASTM(const parameters& params, size_t size,
         std::shared_ptr<double> input_data)
        : vertical_segment(params, size, input_data, vertical_segment::NO_CONFIDENCE_INTERVAL)
    { }

    ASTM(const alta::parameters& params, size_t size)
        : vertical_segment(params, size)
    { }
};

// Parse ASTM header.
// An ASTM header is composed of a KEY in capital letters and a list
// of string/numeric values. An ASTM header finishes after the key VARS.
//
static arguments parse_header(std::istream& in) {

    arguments args;
    std::string key, line;

    while(in.good()) {

        // Get the current key
        in >> key;

        // Extract the current line
        std::getline(in, line);

        // End of header
        if(key == "VARS") {
            args.update(key, std::string("[") + line + std::string("]"));
            break;
        } else {
            args.update(key, line);
        }
    }

    return args;
}

const parameters
compute_parameters(const std::vector<std::string>& vars)
{
    // Update the dimX, dimY and input and output dimension based on the
    // values stored in the VARS list.
    //
    // The input parametrization can take the form:
    //    + theta_i, phi_i, theta_s, phi_s
    //    + theta_i, theta_s, phi_s
    //
    // The output parametrization can take the form:
    //    + Integrated ..
    //    + R, G, B
    //
    unsigned int nX = 0, nY = 0;

    for(auto it=vars.begin(); it!=vars.end(); it++) {
        if(*it == "theta_i" || *it == "phi_i" ||
           *it == "theta_s" || *it == "phi_s") {
            nX += 1;
        } else {
            nY += 1;
        }
    }

    auto out_param = [&]() {
        if(vars.back() == "B") {
            return params::RGB_COLOR;
        } else if(vars.back().compare(0, 10, "Integrated") == 0) {
            return params::INV_STERADIAN;
        } else {
            std::cout << "<<ERROR>> Output format not handled in \'data_astm\'" << std::endl;
            return params::UNKNOWN_OUTPUT;
        }
    };

    return parameters(nX, nY,
                      nX == 4 ? params::SPHERICAL_TL_PL_TV_PV : params::UNKNOWN_INPUT,
                      out_param());
}


ALTA_DLL_EXPORT data* provide_data(size_t size, const parameters& params,
                                   const arguments&)
{
    return new ASTM(params, size);
}

ALTA_DLL_EXPORT data* load_data(std::istream& input, const arguments& args)
{
		std::string line;

		// Parse the header, get the output/input dimensions and params
		arguments header = parse_header(input);
    auto vars = header.get_vec<std::string>("VARS");
    auto params = compute_parameters(vars);

    // Check data size from header
    assert(header.is_defined("NUM_POINTS"));
    int size = header.get_int("NUM_POINTS", 0);

    // Size of the data
    size_t n = params.dimX() + params.dimY();
    double* data = new double[size * n];

    for (size_t i = 0; i < size * n; i += n)
    {
        assert(input.good());

	do {
	    std::getline(input, line);
	} while(line.size() == 0 || line.rfind(',') == std::string::npos);

        std::replace(line.begin(), line.end(), ',', ' ');

        // Create a stream from the line data and extract it as
        // a vec of dim parametrization().dimX() + parametrization().dimY().
        std::stringstream stream(line);
        for(int j = 0; j < n; ++j) {
            stream >> data[i + j];
        }
    }

    return new ASTM(params, size, std::shared_ptr<double>(data));
}

}
