/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2015, 2016 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

/*! \package brdf2data
 *  \ingroup commands
 *  \brief
 *  This command allows to convert a \ref function object to a \ref data 
 *  object. And to save the \ref data object in a file specified by the
 *  \ref data plugin.
 *  \details
 */
#include <core/args.h>
#include <core/data.h>
#include <core/vertical_segment.h>
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
      std::cout << "Usage: brdf2data --input brdf.file --output data.file [--data exporter.so --data-file data.file]" << std::endl ;
      std::cout << "Convert a function object to a data object."<< std::endl ;
      std::cout << std::endl;
      std::cout << "Mandatory arguments:" << std::endl;
      std::cout << "  --input     [filename]" << std::endl;
      std::cout << "  --output    [filename]" << std::endl;
      std::cout << std::endl;
      std::cout << "Optional arguments:" << std::endl;
      std::cout << "  --data      [filename] Name of the data plugin used to save the output" << std::endl ;
      std::cout << "                         data file. If no plugin is defined, the data file" << std::endl ;
      std::cout << "                         will be load using ALTA format but require an" << std::endl;
      std::cout << "                         ALTA file as template." << std::endl ;
      std::cout << "  --data-file [filename] ALTA data file used as a template if no data" << std::endl ;
      std::cout << "                         plugin is specified to export data." << std::endl ;
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
   if(! args.is_defined("data") && ! args.is_defined("data-file")) {
      std::cerr << "<<ERROR>> the data exporter is not defined" << std::endl ;
      return 1 ;
   }
   
   // Get the associated data object and load the file is any
   ptr<data> d;
   if(args.is_defined("data-file")) {
      try
      {
         d = plugins_manager::load_data(args["data-file"], args["data"], args);
      }
      CATCH_FILE_IO_ERROR(args["data-file"]);
   } else {
      d = plugins_manager::get_data(args["data"], 0, parameters(), args);
   }

   // Get the function file
   function* f = NULL;
   f = plugins_manager::load_function(args["input"]);
   if(f == NULL)
   {
      std::cerr << "<<ERROR>> cannot open the function file" << std::endl;
      return 1;
   }

   if(d && f != NULL)
   {
      const bool output_dif = args.is_defined("export-diff");

      vec temp(f->parametrization().dimX());
      for(int i=0; i<d->size(); ++i)
      {
         // Copy the input vector
         vec x = d->get(i);

         // Convert the data to the function's input space.
         if(f->parametrization().input_parametrization() == params::UNKNOWN_INPUT)
         {
            temp = x.head(f->parametrization().dimX());
         }
         else
         {
            params::convert(&x[0],
                            d->parametrization().input_parametrization(),
                            f->parametrization().input_parametrization(),
                            &temp[0]);
         }

         // Eval BRDF
         vec y = f->value(temp);
         x.tail(d->parametrization().dimY()) = (output_dif) ? x.tail(d->parametrization().dimY())-y : y;
         d->set(i, x);
      }  

      // Save data to file
      d->save(args["output"]);
   }  
   else
   {
      std::cerr << "<<ERROR>> cannot import function or export data" << std::endl ;
   }

   return 0 ;
}
