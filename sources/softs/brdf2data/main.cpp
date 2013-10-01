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

int main(int argc, char** argv)
{
	arguments args(argc, argv) ;

	if(args.is_defined("help")) {
		std::cout << "<<HELP>> brdf2data --input brdf.file --func importer.so --output data.file --data exporter.so" << std::endl ;
		std::cout << " - input, output, func, data are mandatory parameters" << std::endl ;
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
		std::cerr << "<<ERROR>> the data exporter is not defined" << std::endl ;
		return 1 ;
	}
/*
	if(! args.is_defined("func")) {
		std::cerr << "<<ERROR>> the function importer is not defined" << std::endl ;
		return 1 ;
	}
*/
	// Import data
	data* d = NULL ;
    d = plugins_manager::get_data(args["data"]) ;

	function* f = NULL;
    f = plugins_manager::get_function(args["input"]);

	// Modify function or data to provide coherent
	// interfaces
//	plugins_manager::check_compatibility(d, f, args);	

	if(d != NULL && f != NULL)
	{
		vec temp(f->dimX());
		for(int i=0; i<d->size(); ++i)
		{
			// Copy the input vector
			vec x = d->get(i);
			params::convert(&x[0], d->parametrization(), f->parametrization(), &temp[0]);

			vec y = f->value(temp);

			for(int j=0; j<d->dimY(); ++j)
			{
				x[d->dimX() + j] = y[j];
			}

			d->set(x);
		}	

		d->save(args["output"]);
	}	
	else
	{
		std::cerr << "<<ERROR>> cannot import function or export data" << std::endl ;
	}

	return 0 ;
}
