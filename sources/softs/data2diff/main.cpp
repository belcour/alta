/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2013, 2014, 2015 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include <core/args.h>
#include <core/data.h>
#include <core/params.h>
#include <core/function.h>
#include <core/fitter.h>
#include <core/plugins_manager.h>

#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>
#include <cstdlib>
#include <cmath>

using namespace alta;

int main(int argc, char** argv)
{
    arguments args(argc, argv) ;

	if(args.is_defined("help")) {
		std::cout << "<<HELP>> data2diff --input data.file --output gnuplot.file --data loader.so --param RUSIN_TH_PH_TD_PD --partial 0" << std::endl ;
		std::cout << " - input, output, data and param are mandatory parameters" << std::endl ;
		std::cout << std::endl << "Available parametrizations: " << std::endl ;
		params::print_input_params();
		return 0 ;
	}

	if(! args.is_defined("input")) {
		std::cerr << "<<ERROR>> the input filename is not defined" << std::endl ;
		return 1 ;
	}
	if(! args.is_defined("output")) {
		std::cerr << "<<ERROR>> the output filename is not defined" << std::endl ;
		return 1 ;
	}
	if(! args.is_defined("data")) {
		std::cerr << "<<ERROR>> the data manager is not defined" << std::endl ;
		return 1 ;
	}
	if(! args.is_defined("param")) {
		std::cerr << "<<ERROR>> the parametrization is not defined" << std::endl ;
		return 1 ;
	}

	// Import data
	ptr<data> d = NULL ;
    d = plugins_manager::get_data(args["data"]) ;
	try
	{
		d->load(args["input"]);
	}
	CATCH_FILE_IO_ERROR(args["input"]);

	// Create output file
	std::ofstream file(args["output"].c_str(), std::ios_base::trunc);

	// Retrieve dimension of differenciation
	int dim = 0;
	params::input p_in = params::parse_input(args["param"]);

	if (args.is_defined("partial"))
	{
		int p_idx = atoi(args["partial"].c_str());
		if (p_idx>=0 && p_idx<params::dimension(p_in)) 
			dim = p_idx;
	}	

	// Get the output data object
    ptr<data> out_d = plugins_manager::get_data(args["data"]);

	if(d != NULL)
	{
		float dt = 0.01;
		vec L(3), V(3), tempParam(params::dimension(p_in)), tempCart(6);
		for(int i=0; i<d->size(); ++i)
		{
			// Copy the input vector
			vec x = d->get(i);
			vec out_x = x;
		
			// Convert input to required param
			params::convert(&x[0], d->parametrization(), p_in, &tempParam[0]);

			// Convert perturbed input to CARTESIAN
			tempParam[dim] += dt;
			params::convert(&tempParam[0], p_in, params::CARTESIAN, &tempCart[0]);
			L[0] = tempCart[0]; L[1] = tempCart[1]; L[2] = tempCart[2];
			V[0] = tempCart[3]; V[1] = tempCart[4]; V[2] = tempCart[5];
			//vec y1 = d->value(L, V);
			
			// Convert perturbed input to CARTESIAN
			tempParam[dim] -= 2.0*dt;
			params::convert(&tempParam[0], p_in, params::CARTESIAN, &tempCart[0]);
			L[0] = tempCart[0]; L[1] = tempCart[1]; L[2] = tempCart[2];
			V[0] = tempCart[3]; V[1] = tempCart[4]; V[2] = tempCart[5];
			//vec y2 = d->value(L, V);
/*
			// Compute the diff vector
			for(int j=0; j<d->dimY(); ++j)
				out_x[d->dimX() + j] = (y1[j]-y2[j])/(2.0*dt);
*/
			// Store it into the output data object
			out_d->set(out_x);

/*
			// Print the input vector	
			for(int j=0; j<d->dimX(); ++j)
				file << x[j] << "\t";

			// And the diff values
			for(int j=0; j<d->dimY(); ++j)
				file << (y1[j]-y2[j])/(2.0*dt) << "\t";

			file << std::endl;
*/
		}	

		out_d->save(args["output"]);

		file.close();
	}	
	else
	{
		std::cerr << "<<ERROR>> cannot export data" << std::endl ;
	}

	return 0 ;
}
