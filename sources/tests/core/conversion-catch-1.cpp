/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2017 CNRS
    
   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  
*/

/*
 *  Initial author : Romain Pacanowski @ institutoptique.fr
 */


//Catch Framework
// CATCH SHOULD PROVIDE MAIN
#define CATCH_CONFIG_MAIN
#include "catch.hpp"


//ALTA
#include <core/args.h>
#include <core/data.h>
#include <core/params.h>
#include <core/common.h>
#include <core/plugins_manager.h>
#include <core/vertical_segment.h>



TEST_CASE( " Conversion from MERL data to ALTA with all values", "[MERL to ALTA. All Values]")
{
  using namespace alta;

  arguments args;


  std::string const DATA_FILE1 = "./build/tests/gold-metallic-paint.binary";

  //Load data from File
  ptr<data> d_in;

  try
  {
      d_in = plugins_manager::load_data(DATA_FILE1.c_str(),"data_merl", args) ;
  }
  catch(...)
  {
    //Could not load file ?
    REQUIRE(false);
  }


  REQUIRE( d_in->size() == 90*90*180 );

  REQUIRE( d_in->parametrization().dimX() == 3); // RUSINKIEWICZ PARAMETRIZATION AS INPUT
  REQUIRE( d_in->parametrization().dimY() == 3); // COLOR RGB AS VALUES

}

TEST_CASE( "ALTA FILTERING TEST by keeping only positive values for Y", "[ALTA Y Filtering]")
{
  using namespace alta;


  arguments args;
  args.update("ymin","[0.0, 0.0, 0.0]");

  std::string const DATA_FILE1 = "./build/tests/gold-metallic-paint.alta";

  std::cout << " ATTEMTPING TO LOAD  = " << DATA_FILE1 << std::endl;


  //Load data from File
  ptr<data> d_in;
 //Try to convert to VerticalSegment. Should work or we have a problem
  ptr<vertical_segment> vs;
 
  try
  {
      d_in = plugins_manager::load_data(DATA_FILE1.c_str(),"", args) ;
      vs = dynamic_pointer_cast< vertical_segment >( d_in );

  }
  catch(...)
  {
    //Could not load file ?
    REQUIRE(false);
  }

  // vec y_min = args.get_vec("ymin",3, -1);
  // std::cout << " y_min = " << y_min << std::endl;
  // std::cout <<  " is vec for ymin " << args.is_vec("ymin");

  // std::cout << " d_in->parametrization().dimX()  = " << d_in->parametrization().dimX() << std::endl;
  // std::cout << " d_in->size() = " << d_in->size() << std::endl;

  
  for (int i = 0; i < d_in->size(); ++i)
  {
    vec x_y = d_in->get(i);

    for(int d=0; d < d_in->parametrization().dimY(); d++ )
    {
      double y = x_y( d_in->parametrization().dimX() + d );
      REQUIRE( y  >= 0.0 );  
    }
  }


  Eigen::Map<RowMatrixXd>  data_as_matrix = vs->matrix_view();
  
  RowMatrixXd const & y_values =  data_as_matrix.block(0, d_in->parametrization().dimX(), d_in->size(), d_in->parametrization().dimY()  );

  REQUIRE( y_values.all() >= 0.0 );


}

TEST_CASE( "ALTA FILTERING TEST by keeping only values for phi_diff and theta_diff near zero", "[ALTA X Filtering]")
{
    using namespace alta;

  double const THETA_H_MAX_VALUE = 1.55;

  arguments args;
  args.update("ymin","[0.0, 0.0, 0.0]");
  args.update("max","[1.55, 0.01, 0.01]");

  //args.create_arguments( "--ymin [0.0, 0.0, 0.0] --max [1.5, 0.01, 0.01] --param RUSIN_TH_TD_PD");
  

  std::string const DATA_FILE1 = "./build/tests/gold-metallic-paint.alta";

  //Load data from File
  ptr<data> d_in;
  ptr<vertical_segment> vs;

  try
  {
      d_in = plugins_manager::load_data(DATA_FILE1.c_str(),"", args) ;
      vs = dynamic_pointer_cast< vertical_segment >( d_in );
  }
  catch(...)
  {
    //Could not load file ?
    REQUIRE(false);
  }

  Eigen::Map<RowMatrixXd>  data_as_matrix = vs->matrix_view();
  
  RowMatrixXd const & x_values =  data_as_matrix.block(0, 0, d_in->size(), d_in->parametrization().dimX()  );

  // std::cout << " d_in->parametrization().dimX()  = " << d_in->parametrization().dimX() << std::endl;
  // std::cout << " d_in->size()  = " << d_in->size() << std::endl;
  // std::cout << " Nb Rows in x_values = " << x_values.rows() << std::endl;
  // std::cout << " Nb Colums in x_avlues " << x_values.cols() << std::endl;

  vec x_max = args.get_vec("xmax",3, -1);

  // std::cout << " xmax = " << x_max << std::endl;
  // std::cout <<  " is vec for xmax " << args.is_vec("xmax");

  // for( int i=0; i < x_values.rows(); i++)
  // {
  //   std::cout << x_values(i,0) << " " << x_values(i,1) << " " << x_values(i,2)  
  //             << std::endl;
  // } 


  Eigen::ArrayXd col_0 = x_values.col( 0 ).array();

  //CHECK ALL first X values are supposed to be less than 1.5
  //Thats theta_h for MERL material
  REQUIRE( (col_0 <= THETA_H_MAX_VALUE).all() );

  //CHECK ALL Second X values are supposed to be less than 0.01 
  // Thats theta_diff for MERL Material
  REQUIRE( (x_values.col( 1 ).array()   <= 0.01).all() ) ;
  
  //CHECK ALL Third X values are supposed to be less than 0.01
  // Thats phi_diff for MERL Material
  REQUIRE( (x_values.col( 2 ).array()   <= 0.01).all() ) ;

}