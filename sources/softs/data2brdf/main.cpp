/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2015, 2016 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

/*! \package data2brdf
 *  \ingroup commands
 *  \brief
 *  This command allows to convert a \ref data object to a \ref function
 *  object.
 *  \details
 *
 *  Note, each plugin you load can have its own list of arguments. They
 *  are not listed here. Please refer to the documentation of the plugin
 *  to know how to fine tune it.
 *
 *  <h3>Parameters</h3>
 *  <ul>
 *		<li><b>\-\-help</b> print the help and exit.
 *		<li><b>\-\-input <i>filename</i></b> data file to be loaded using
 *		the data plugin specified by the <b>\-\-data <i>filename</i></b>
 *		option. If no plugin is specified, the data file will be loaded
 *		using a \ref vertical_segment data object.</li>
 *		<li><b>\-\-output <i>filename</i></b> function file to be exported
 *		in the format specified by <b>\-\-export</b>. If no export argument
 *		is given, the function will be exported in ALTA \ref format.
 *  </ul>
 *
 *  <h3>Plugins</h3>
 *  <ul>
 *		<li><b>\-\-fitter <i>plugin.lib</i></b> \ref fitter library to be used.
 *		Note that the fitter plugin must be compatible with the \a function
 *		to be used. We make the convention that the filename of the plugin
 *		expose its compatibility. For example,
 *		<i>rational_function_legendre.lib</i> expose a rational function,
 *		compatible with a rational fitter and will not be compatible with
 *		nonlinear fitters.
 *		<li><b>\-\-data <i>plugin.lib</i></b> \ref data library to be used.
 *		<li><b>\-\-func <i>plugin.lib</i></b> \ref function library to be
 *		used.
 *  </ul>
 */
#include <core/args.h>
#include <core/data.h>
#include <core/function.h>
#include <core/fitter.h>
#include <core/plugins_manager.h>

#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>
#include <cstdlib>

#ifdef __GLIBC__
#include <fenv.h>
#endif

using namespace alta;

int main(int argc, char** argv)
{
    arguments args(argc, argv) ;

#ifdef __GLIBC__
	// feenableexcept(FE_DIVBYZERO | FE_OVERFLOW | FE_INVALID);
#endif

	 if(args.is_defined("help")) {
		std::cout << "Usage: data2brdf [options] --input data.file --output data.file" << std::endl ;
		std::cout << "Convert a data object to a function object using a fitting procedure."<< std::endl ;
		std::cout << std::endl;
		std::cout << "Mandatory arguments:" << std::endl;
		std::cout << "  --input    [filename]" << std::endl;
		std::cout << "  --output   [filename]" << std::endl;
		std::cout << "  --fitter   [filename]" << std::endl;
		std::cout << std::endl;
		std::cout << "Optional arguments:" << std::endl;
		std::cout << "  --func     [filename]  Name of the function plugin. If not defined, a" << std::endl ;
		std::cout << "                         monomial rational function will be used." << std::endl ;
		std::cout << "  --data     [filename]  Name of the data plugin used to load the input" << std::endl ;
		std::cout << "                         data file. If no plugin is defined, the data file" << std::endl ;
		std::cout << "                         will be load using ALTA format." << std::endl ;
		return 0 ;
	}

    ptr<fitter> fit = plugins_manager::get_fitter(args["fitter"]) ;
    if(!fit)
    {
        std::cerr << "<<ERROR>> unable to load the fitter plugin \"" << args["fitter"] << "\"" << std::endl;
        return 1;
    }

    if(args.is_defined("available_params"))
    {
        params::print_input_params();
        return 0;
    }

    if(! args.is_defined("input")) {
        std::cerr << "<<ERROR>> the input filename is not defined" << std::endl ;
        return 1 ;
    }
    if(! args.is_defined("output")) {
        std::cerr << "<<ERROR>> the output filename is not defined" << std::endl ;
        return 1 ;
    }

    fit->set_parameters(args) ;

    ptr<data> d;
		try
		{
        d = plugins_manager::load_data(args["input"], args["data"], args);
		}
		CATCH_FILE_IO_ERROR(args["input"]);

    ptr<function> f = ptr<function>(plugins_manager::get_function(args, d->parametrization()));

    if(!f || !d)
    {
        std::cerr << "<<ERROR>> no function or data object correctly defined" << std::endl;
        return 1;
    }

    if(d->size() == 0)
    {
        std::cerr << "<<ERROR>> no data loaded, please check you input file" << std::endl;
        return 1;
    }

    // Check the compatibility between the data and the function
    plugins_manager::check_compatibility(d, f, args);


    // Start a timer
    timer time ;
    time.start() ;

    // Fit the data
    bool is_fitted = fit->fit_data(d, f, args) ;

    // Get the fitting duration
    time.stop();

    // Display the result
    if(is_fitted)
    {
        std::cout << "<<INFO>> total time: " << time << std::endl ;

        double L2   = f->L2_distance(d);
        double Linf = f->Linf_distance(d);
        std::cout << "<<INFO>> L2   distance to data = " << L2   << std::endl;
        std::cout << "<<INFO>> Linf distance to data = " << Linf << std::endl;

        // Export the L2 and Linf values to the command line
        std::stringstream L2string, Linfstring;
        L2string << L2; Linfstring << Linf;
        args.update("L2",   L2string.str());
        args.update("Linf", Linfstring.str());

        f->save(args["output"], args) ;

        return 0;
    }
    else
    {
        std::cout << "<<ERROR>> data2brdf: unable to fit the data" << std::endl ;
        return 1;
    }


    return 0 ;
}
