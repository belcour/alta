/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2014, 2017 CNRS
   Copyright (C) 2013, 2014, 2015, 2016 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

/*! \package data2data
 *  \ingroup commands
 *  \brief
 *  This command transform a \ref data object into another data object.
 *  \details
 *  This command is useful to change the parametrization of a data file,
 *  or to perform interpolation of sample values (i.e. to fill gaps).
 *
 *  <h3>Usage</h3>
 *  \verbatim
       data2data --input data.file --output data.file --out-data exporter.so --in-data importer.so
    \endverbatim
 *
 *  <h3>Parameters</h3>
 *  <ul>
 *		<li><b>\-\-input <i>filename</i></b> data file to be loaded. The
 *		loading plugin is defined using the option <b>\-\-in-data <i>
 *		filename</i></b>.
 *		</li>
 *		<li><b>\-\-in-data <i>filename</i></b> specify the data plugin
 *		used to load the input file. If this option is not specified, the
 *		loading plugin will be a \ref vertical_segment plugin. \note If
 *		the input plugin is not interpolating, like \ref vertical_segment,
 *		you can only use the reparametrization tool.</li>
 *		<li><b>\-\-output <i>filename</i></b> the resulting data file.
 *		</li>
 *		<li><b>\-\-out-data <i>filename</i></b></li> specify the data plugin
 *		used to export the data. This parameter is optional. If not defined,
 *		the output format will be ALTA's \ref format.
 *		</li>
 *  </ul>
 */
#include <core/args.h>
#include <core/data.h>
#include <core/params.h>
#include <core/function.h>
#include <core/fitter.h>
#include <core/plugins_manager.h>
#include <core/vertical_segment.h>

#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>
#include <cstdlib>
#include <cmath>

using namespace alta;

static parameters compute_parameters(const data& d_in,
                                     const arguments& args)
{
    params::input param = params::parse_input(args["param"]);

    if(param == params::UNKNOWN_INPUT
       && d_in.parametrization().input_parametrization() != params::UNKNOWN_INPUT)
    {
        std::cout << "<<DEBUG>> using the input parametrization of the input file for the output file as well." << std::endl;
        param = d_in.parametrization().input_parametrization();
    }
    else if(param == params::UNKNOWN_INPUT)
    {
        std::cerr << "<<ERROR>> no parametrization defined for input and output files." << std::endl;
        exit(EXIT_FAILURE);
    }

    // Check the output parametrization
    // For now we propagate OUTPUT_PARAM since there isn't any way to do conversion
    params::output const & out_param = d_in.parametrization().output_parametrization();

    return alta::parameters(params::dimension(param),
                            d_in.parametrization().dimY(),
                            param,
                            out_param);
}

int main(int argc, char** argv)
{
	arguments args(argc, argv) ;

	if(args.is_defined("help")) {
		std::cout << "Usage: data2data [options] --input data.file --output data.file" << std::endl ;
		std::cout << "Convert a data object to another data object."<< std::endl ;
		std::cout << std::endl;
		std::cout << "Mandatory arguments:" << std::endl;
		std::cout << "  --input    [filename]" << std::endl;
		std::cout << "  --output   [filename]" << std::endl;
		std::cout << std::endl;
		std::cout << "Optional arguments:" << std::endl;
		std::cout << "  --out-data [filename]  Name of the plugin used to save the outputed data file" << std::endl ;
    std::cout << "                         If none is provided, data2data will export in ALTA" << std::endl ;
		std::cout << "                         by default." << std::endl ;
		std::cout << "  --in-data  [filename]  Name of the plugin used to load the input data file" << std::endl ;
    std::cout << "                         If none is provided, data2data will import in ALTA" << std::endl ;
		std::cout << "                         by default." << std::endl ;
		std::cout << "  --param    [NAME]      Name of the parametrization used to output data when" << std::endl;
		std::cout << "                         No output data plugin is specified. Please see " << std::endl;
		std::cout << "                         --help-params for the list of available " << std::endl ;
		std::cout << "                         parametrizations." << std::endl ;
		std::cout << "  --data-correct-cosine  Divide the value of the data points by the product of" << std::endl;
		std::cout << "                         the light and view vector dot product with the normal." << std::endl ;
    std::cout << "   --all-values          Export all data regardless of their physical validity." << std::endl;
		std::cout << std::endl;
		std::cout << "Helps:" << std::endl;
		std::cout << "  --help                 Display this help." << std::endl;
		std::cout << "  --help-params          Display the available input parametrizations." << std::endl;
    std::cout << "  --verbose              Enable more information output. Helpful to debug." << std::endl;
		return EXIT_SUCCESS;
	}
	if(args.is_defined("help-params")) {
		params::print_input_params();
		return EXIT_SUCCESS;
	}

	if(! args.is_defined("input")) {
		std::cerr << "<<ERROR>> The input filename is not defined" << std::endl ;
		return EXIT_FAILURE;
	}
	if(! args.is_defined("output")) {
		std::cerr << "<<ERROR>> The output filename is not defined" << std::endl ;
		return EXIT_FAILURE;
	}
	/*
	if(! args.is_defined("out-data")) {
		std::cerr << "<<ERROR>> the data exporter is not defined" << std::endl ;
		return 1 ;
	}
	*/


	// Import data
	ptr<data> d_in;

	try
	{
      alta::timer load_timer;
      load_timer.start();

      d_in = plugins_manager::load_data(args["input"], args["in-data"], args) ;

      load_timer.stop();
      std::cout << "<<INFO>> Data Loaded in "<< load_timer << std::endl;
	}
	CATCH_FILE_IO_ERROR(args["input"]);

	if(!d_in)
	{
		std::cout << "<<INFO>> Input data will be treated as ALTA format" << std::endl;
	}


	ptr<data> d_out = plugins_manager::get_data(args["out-data"],
                                              d_in->size(),
                                              compute_parameters(*d_in, args),
                                              args) ;
	if(!d_out)
	{
		std::cout << "<<INFO>> Data will be outputed to ALTA format" << std::endl;
	}

	if(!d_in || !d_out)
	{
		std::cerr << "<<ERROR>> Cannot import or export data" << std::endl ;
		return EXIT_FAILURE;
	}


  std::cout << "<<INFO>> Conversion from " << params::get_name(d_in->parametrization().input_parametrization())
           << " to " << params::get_name(d_out->parametrization().input_parametrization()) << std::endl;

  std::cout << "<<INFO>> Dimensions for  Input Data [X,Y] = " << d_in->parametrization().dimX()
            << ", " << d_in->parametrization().dimY() << std::endl;


	if(dynamic_pointer_cast<vertical_segment>(d_out) || args.is_defined("splat"))
	{
		std::cout << "<<INFO>> Dimensions for Output Data [X,Y] = " << d_out->parametrization().dimX() << ", " << d_out->parametrization().dimY() << std::endl;

    if(args.is_defined("verbose"))
    {
      std::cout << "<<DEBUG>> d_in->parametrization().input_parametrization(), = " << params::get_name(d_in->parametrization().input_parametrization() ) << std::endl;
      std::cout << "<<DEBUG>> d_in->parametrization().output_parametrization(), = " << params::get_name( d_in->parametrization().output_parametrization() ) << std::endl;
      std::cout << "<<DEBUG>> d_out->parametrization().input_parametrization(), = " << params::get_name( d_out->parametrization().input_parametrization() ) << std::endl;
      std::cout << "<<DEBUG>> d_out->parametrization().output_parametrization() = " << params::get_name( d_out->parametrization().output_parametrization() ) << std::endl;
    }
    
    alta::timer save_timer;
    save_timer.start();

    unsigned int nb_invalid_configs = 0;
		vec temp(d_out->parametrization().dimX() + d_out->parametrization().dimY());
		for(int i=0; i<d_in->size(); ++i)
		{

			// Copy the input vector
			vec x = d_in->get(i);
			params::convert(&x[0],
                      d_in->parametrization().input_parametrization(),
                      d_out->parametrization().input_parametrization(),
                      &temp[0]);
      
      // Todo change this at some point because we cannot handle BTDF stuff here
      // Check if this  BRDF parametrization is valid (over the hemisphere)
      if( params::is_below_hemisphere( &temp[0], d_out->parametrization().input_parametrization()  ) )
      {
        nb_invalid_configs++;

        if( args.is_defined("all-values")) // If the user wants to save invalid values as well
        {
          // Converts the output values from vector x to the output values of temp vector
          params::convert(&x[d_in->parametrization().dimX()],
                          d_in->parametrization().output_parametrization(),
                          d_in->parametrization().dimY(),
                          d_out->parametrization().output_parametrization(),
                          d_out->parametrization().dimY(),
                          &temp[d_out->parametrization().dimX()]);
          // SAVING THE INVALID CONFIGURATION
          d_out->set(i, temp);
        }

      }
      else // Apply the conversion
      {
        params::convert(&x[d_in->parametrization().dimX()],
                        d_in->parametrization().output_parametrization(),
                        d_in->parametrization().dimY(),
                        d_out->parametrization().output_parametrization(),
                        d_out->parametrization().dimY(),
                        &temp[d_out->parametrization().dimX()]);
        d_out->set(i, temp);
      }
		}
    save_timer.stop();

    if( nb_invalid_configs > 0 )
    {
      std::cout << "<<INFO>> Number of Invalid Configurations = " << nb_invalid_configs << " over " << d_in->size() << " configurations" << std::endl;
    }

    std::cout << "<<INFO>> Data saved to file in "<< save_timer << std::endl;

	}
	else
	{
		if(d_out->parametrization().output_parametrization() != d_in->parametrization().output_parametrization())
		{
			std::cerr << "<<WARNING>> data types have different output parametrizations." << std::endl;
			std::cerr << "            This is currently not handled properly by ALTA." << std::endl;
		}

		if(d_out->parametrization().dimY() != d_in->parametrization().dimY())
		{
			std::cerr << "<<WARNING>> data types have different output dimensions (" << d_in->parametrization().dimY()
			          << " and " << d_out->parametrization().dimY() << ")." << std::endl;
			std::cerr << "            This is currently not handled properly by ALTA." << std::endl;
		}

    unsigned int stats_incorrect = 0;

		#pragma omp parallel for
		for(int i=0; i<d_out->size(); ++i)
		{
			vec temp(d_in->parametrization().dimX());
			vec cart(6);
			vec y(d_in->parametrization().dimY());

			// Copy the input vector
			vec x = d_out->get(i);
			params::convert(&x[0],
                      d_out->parametrization().input_parametrization(),
                      params::CARTESIAN, &cart[0]);


         // Check if the output configuration is below the hemisphere when
         // converted to cartesian coordinates. Note that this prevent from
         // converting BTDF data.
			if(cart[2] >= 0.0 || cart[5] >= 0.0) {
				params::convert(&cart[0], params::CARTESIAN,
                        d_in->parametrization().input_parametrization(),
                        &temp[0]);
				y = d_in->value(temp);
			} else {
        ++stats_incorrect;
				y.setZero();
			}

         // Convert the value stored in the input data in the value format of
         // the output data file.
			params::convert(&y[0],
                      d_in->parametrization().output_parametrization(),
                      d_in->parametrization().dimY(),
                      d_out->parametrization().output_parametrization(), d_out->parametrization().dimY(),
                      &x[d_out->parametrization().dimX()]);

			d_out->set(i, x);
		}

      if(stats_incorrect > 0) {
         std::cerr << "<<DEBUG>> Number of incorrect configuration: "
                   << stats_incorrect << " / " << d_out->size() << std::endl;
      }
	}

	d_out->save(args["output"]);
	return EXIT_SUCCESS;
}
