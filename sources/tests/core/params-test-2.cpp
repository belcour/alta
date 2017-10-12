/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2016 CNRS, romain.pacanowski@institutoptique.fr


   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

// ALTA includes
#include <core/params.h>
#include <tests.h>

using namespace alta;
using namespace alta::tests;

#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <random>

typedef std::mt19937_64 GENERATOR_TYPE;

/* Test different configurations for the Half / Cartesian parametrization.  */
int main(int argc, char** argv) 
{

  unsigned int const NB_TESTS = 1000;
  unsigned int nb_pass = 0;
  // Number of failed tests
  unsigned int nb_fail = 0;

  // Testing is Above with Cartesians Configurations
  std::uniform_real_distribution<> distribution(-1.0, 1.0);
  GENERATOR_TYPE                      generator( std::time(NULL) );

  vec input_lv(6);
  for (unsigned int i = 0; i < NB_TESTS; ++i)
  {
     double const l_x =  distribution( generator );
     double const l_y = distribution( generator );
     double const l_z = std::sqrt(std::max( 1.0 - l_x*l_x - l_y*l_y, 0.0) );

     double const v_x = distribution( generator );
     double const v_y = distribution( generator );
     double const v_z = std::sqrt( std::max( 1.0 - v_x*v_x - v_y*v_y, 0.0) );

     input_lv(0) = l_x; input_lv(1) = l_y; input_lv(2) = l_z;
     input_lv(3) = v_x; input_lv(4) = v_y; input_lv(5) = v_z;

     if(  params::is_above_hemisphere( &input_lv[0], params::CARTESIAN ) )
     {
        nb_pass++;
     }
     else
     {
        nb_fail++;
     }
  }


  if(nb_fail > 0) {
      std::cerr << "<<ERROR>> " << nb_fail << " out of " << NB_TESTS
                << " configurations fail the test is_above_hemisphere "
                << std::endl;
  }

  return nb_fail > 0 ? EXIT_FAILURE : EXIT_SUCCESS;
}
