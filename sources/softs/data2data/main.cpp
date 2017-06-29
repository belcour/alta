/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2014, 2017 CNRS
   Copyright (C) 2013, 2014, 2015, 2016, 2017 Inria

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
 *  or to perform interpolation of sample values (i.e., to fill gaps).
 *
 *  <h3>Usage</h3>
 *  \verbatim
       data2data --input data.file --output data.file --out-data  output_plugin --in-data input_plugin  OPTIONS
    \endverbatim
 *
 *  <h3>Parameters</h3>
 *  <ul>
 *		<li><b>\-\-input <i>filename</i></b> data file to be loaded. The
 *		loading plugin is defined using the option <b>\-\-in-data <i>
 *		filename</i></b>.
 *		</li>
 *		<li><b>\-\-in-data <i>filename</i></b> specifies the data plugin (cf. \ref datas)
 *		used to load the input file. If this option is not specified, the
 *		loading plugin will be a \ref alta::vertical_segment plugin. \note If
 *		the input plugin is not interpolating, like \ref alta::vertical_segment,
 *		you can only use the reparametrization tool.</li>
 *		<li><b>\-\-output <i>filename</i></b> the resulting data file.
 *		</li>
 *		<li><b>\-\-out-data <i>filename</i></b></li> specifies the data plugin
 *		used to export the data. This parameter is optional. If not defined,
 *		the output format will be ALTA's \ref format.
 *		</li>
 *  </ul>
 *  <h3>Options</h3>
 *  <ul>
 *    <li>
 *      <b>--params <i>NAME</i> </b> specifies the name of the parametrization to be used
 *               no output data plugin is specified.
 *               Please see --help-params for the list of available parametrizations.
 *    </li>
 *
 *    <li>
 *      <b>--all-values </b> exports all data regardless of their physical validity
 *    </li>
 *
 *    <li>
 *      <b>--data-correct-cosine</b> divides the value of the data points by the product of
 *       the light and view vector dot product with the normal
 *    </li>
 *   </ul>
 *
 * <h4>Data Filtering  (currently only supported with ALTA format) </h4>
 *  <ul>
 *    <li>
 *      <b>--min [min_1, min_2, ...] </b> discards all samples which input values are lower than min_1, min_2,...
 *    </li>
 *    <li>
 *      <b>--max  [max_1, max_2, ...] </b> discards all samples which input values are greater than max_1, max_2,...
 *    </li>
 *    <li>
 *      <b>--ymin [ymin_1, ymin_2, ...] </b> discards all samples which output values are lower than min_1, min_2,...
 *    </li>
 *    <li>
 *      <b> --ymax [ymax_1, ymax_2, ...] </b> discards all samples which ouptut values are greater than max_1, max_2,...
 *    </li>
 *
 *  </ul>
 *
 *  <h3>Usage Examples</h3>
 *    <h4> Converting from MERL to plain text ALTA file format </h4>
 *      \verbatim
        data2data   --in-data data_merl --input gold-metallic-paint.binary --output gold-metallic-paint.alta  --param RUSIN_TH_TD_PD  --all-values
        \endverbatim
 *    <h4> Converting from ALTA to ALTA format by keeping only positive values</h4>
 *      \verbatim
       data2data  --input gold-metallic-paint_all_values.alta --output gold-metallic-paint_only_valid.alta  --param RUSIN_TH_TD_PD  --ymin [0.0, 0.0, 0.0]
        \endverbatim
 *    <h4> Filtering on input dimensions. </h4>
 *      \verbatim
        data2data  --input gold-metallic-paint_only_valid.alta  --output gold-metallic-paint_filtered_tH.alta  --param RUSIN_TH_TD_PD  --max [1.5, 0.01, 0.01]
        \endverbatim
 *
 */
#include <core/args.h>
#include <core/data.h>
#include <core/data_storage.h>
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

// Convert D_IN to D_OUT without doing any interpolation.  That is, call the
// 'get' method of D_IN to grab its 'x' and 'y' values, and write those to
// D_OUT.  When ALL_VALUES is true, disable filtering of values that are not
// in the definition domain.
static void convert_vertical_segment(ptr<data> d_in, ptr<data> d_out,
                                     bool verbose = false,
                                     bool all_values = false)
{
    std::cout << "<<INFO>> Dimensions for Output Data [X,Y] = " << d_out->parametrization().dimX() << ", " << d_out->parametrization().dimY() << std::endl;

    if (verbose)
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

        // TODO: Change this at some point because we cannot handle BTDF stuff here
        // Check if this  BRDF parametrization is valid (over the hemisphere)
        if( params::is_below_hemisphere( &temp[0], d_out->parametrization().input_parametrization()  ) )
        {
            nb_invalid_configs++;

            if (all_values) // If the user wants to save invalid values as well
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

// Convert D_IN to D_OUT, possibly interpolating the value in D_IN using the
// 'value' method of D_IN.  Assume that D_OUT has its 'x' values already
// initialized, and provide the interpolated 'y' value that correspond.
static void convert_with_interpolation(ptr<data> d_in, ptr<data> d_out)
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
        if(cart[2] >= 0.0 && cart[5] >= 0.0) {
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

// Convert D_IN to D_OUT according to ARGS.
static void convert(ptr<data> d_in, ptr<data> d_out, const arguments& args)
{
	if(dynamic_pointer_cast<vertical_segment>(d_out) || args.is_defined("splat"))
      convert_vertical_segment(d_in, d_out,
                               args.is_defined("verbose"),
                               args.is_defined("all-values"));
	else
      convert_with_interpolation(d_in, d_out);
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
		std::cout << "                         by default.  The name \"alta-binary\" specifies ALTA's" << std::endl ;
		std::cout << "                         native binary format." << std::endl;
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


  // "alta-binary" denotes an output format but not a class, so we need to
  // special-case it.
  const std::string data_class =
      args["out-data"] == "alta-binary"
      ? "vertical_segment" : args["out-data"];

  ptr<data> d_out = plugins_manager::get_data(data_class,
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

  convert(d_in, d_out, args);

  // Special-case ALTA's binary output format.  TODO: In the future, 'save'
  // should no longer be a method and we'd have an output format lookup
  // function in plugin_manager.
  if (args["out-data"] == "alta-binary")
  {
      try
      {
          std::ofstream out;
          out.exceptions(std::ios_base::failbit);
          out.open(args["output"]);
          save_data_as_binary(out, *d_out);
      }
      CATCH_FILE_IO_ERROR(args["output"]);
  }
  else
      d_out->save(args["output"]);

	return EXIT_SUCCESS;
}
