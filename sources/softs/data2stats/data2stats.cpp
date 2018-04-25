/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2016 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */


// STL includes
#include <cstdlib>

// ALTA includes
#include <core/common.h>
#include <core/args.h>
#include <core/data.h>
#include <core/vertical_segment.h>
#include <core/function.h>
#include <core/fitter.h>
#include <core/plugins_manager.h>
#include <core/metrics.h>

// Eigen includes
#include <Eigen/Core>

using namespace alta;

//TODO: Move to ALTA CORE
//TODO: Documentation
//TODO:  Add the relative error computation
class Error
{
   public:
      static vec meanSquareError(ptr<data> const & data, ptr<function> const & f )
      {
         vec mse = vec::Zero( f->parametrization().dimY() );

         vec dat = vec::Zero( data->parametrization().dimX() + data->parametrization().dimY() );

         //Note that: data_x = dat.head( data->parametrization().dimX() );
         //           data_y = dat.tail( data->parametrization().dimY() );
         vec f_y = vec::Zero( f->parametrization().dimY() );

         for(unsigned int i=0; i < (unsigned int) data->size(); i++)
         {
            dat = data->get(i);

            f_y = f->value(dat.head( data->parametrization().dimX() ) );

            mse += (f_y - dat.tail( data->parametrization().dimY() ) ).array().square().matrix();
         }

         return mse * (1.0 / data->size() );

      }

      static vec rootMeanSquareError(ptr<data> const & data, ptr<function> const & f )
      {
         vec const rmse = meanSquareError(data, f);
         return rmse.cwiseSqrt();
      }

};


// TODO : Write the documentation
// TODO: Move this class to ALTA Core
// The methods provided by this class assume that the parametrization of
// the data and function are the same.
class Norm
{
   public:
      static vec L2( ptr<data> const & data, ptr<function> const & f )
      {
         vec l2 = vec::Zero( f->parametrization().dimY() );

         vec dat = vec::Zero( data->parametrization().dimX() + data->parametrization().dimY() );

         //Note that: data_x = dat.head( data->parametrization().dimX() );
         //           data_y = dat.tail( data->parametrization().dimY() );
         vec f_y = vec::Zero( f->parametrization().dimY() );

         for(unsigned int i=0; i < (unsigned int) data->size(); i++)
         {
            dat = data->get(i);

            f_y = f->value(dat.head( data->parametrization().dimX() ) );

            l2 += ( f_y - dat.tail( data->parametrization().dimY() )).array().square().matrix();

         }
         return l2.cwiseSqrt();
      }

      static vec weightedL2( ptr<data> const & data, ptr<function> const & f, vec const & weights)
      {
         vec l2 = vec::Zero( f->parametrization().dimY() );

         vec dat = vec::Zero( data->parametrization().dimX() + data->parametrization().dimY() );

         for(unsigned int i=0; i < (unsigned int) data->size(); i++)
         {
            dat = data->get(i);
            l2 += ((f->value(dat.head( data->parametrization().dimX() )) - dat.tail( data->parametrization().dimY() )) * weights[i]).array().square().matrix();

         }
         return l2.cwiseSqrt();
      }

      static vec LInf( ptr<data> const & data, ptr<function> const & f )
      {
         vec max = vec::Zero( f->parametrization().dimY() );

         vec dat = vec::Zero( data->parametrization().dimX() + data->parametrization().dimY() );

         for(unsigned int i=0; i < data->size(); i++)
         {
            dat = data->get(i);
            //vec const data_x = dat.head( data->parametrization().dimX() );
            //vec const data_y = dat.tail( data->parametrization().dimY() );

            max =  max.cwiseMax( (f->value(dat.head( data->parametrization().dimX() )) - dat.tail( data->parametrization().dimY() )).cwiseAbs() );

         }
         return max;
      }

      static vec L1( ptr<data> const & data, ptr<function> const & f )
      {
         vec L1 = vec::Zero( f->parametrization().dimY() );

         vec dat = vec::Zero( data->parametrization().dimX() + data->parametrization().dimY() );

         for(unsigned int i=0; i < (unsigned int ) data->size(); i++)
         {
            dat = data->get(i);
            //vec const data_x = dat.head( data->parametrization().dimX() );
            //vec const data_y = dat.tail( data->parametrization().dimY() );

            L1 += (f->value(dat.head( data->parametrization().dimX() )) - dat.tail( data->parametrization().dimY() )).cwiseAbs();

         }
         return L1;
      }

      /*!
       * \brief Compute the Lp norm from a given p (with \f$ p >= 1 \f$)
       * \f$ \big(\sum_{i=1}^{i=n} |x_i|^p\big)^{1/p} \f$
       */
      static vec Lp(ptr<data> const & data, ptr<function> const & f, double p)
      {
         //if p is not superior or equal to 1.0 this formula does not define a norm
         assert( p >= 1.0 );

         // This is another strategy to compute the norm
         // First we put all the data into Eigen::Array
         // Then we perform the computation on the big array directly
         Eigen::ArrayXd  a_Lp = Eigen::ArrayXd::Zero( f->parametrization().dimY() );
         Eigen::ArrayXXd  all_data_y( data->size(), f->parametrization().dimY() );
         all_data_y.setZero(data->size(), f->parametrization().dimY());

         Eigen::ArrayXXd  all_y_fx( data->size(), f->parametrization().dimY() );
         all_y_fx.setZero(data->size(), f->parametrization().dimY());

         //timer chrono;
         //chrono.start();
         for(unsigned int i=0; i < (unsigned int) data->size(); i++)
         {
            Eigen::ArrayXd const dat = data->get(i);
            Eigen::ArrayXd const data_x = dat.head( data->parametrization().dimX() );
            Eigen::ArrayXd const data_y = dat.tail( data->parametrization().dimY() );
            Eigen::ArrayXd const y_fx = f->value(data_x);

            all_data_y.row(i) = data_y;
            all_y_fx.row(i)   = y_fx;

            //a_Lp += pow( abs(y_fx - data_y), p );
         }
         // chrono.stop();
         // std::cout << "<<INFO>> Time to demultiplex data " << chrono << std::endl;
         // chrono.reset();

         Eigen::ArrayXd pre_cal = pow(abs(all_data_y - all_y_fx), p).colwise().sum();


         double const one_over_p = 1.0 / p;

         return pre_cal.pow( one_over_p );

         //return Lp.array().pow( one_over_p );

         //Implici conversion from Eigen::Array to Matrix
         return a_Lp.pow( one_over_p );

      }

};

//Returns true if the conversion was successful
bool convertDataToFunctionParam(ptr<data> const & d,
      ptr<function> const & f,
      bool & conversion_necessary,
      vertical_segment* & converted_data )
{
   conversion_necessary = true;
   if( d->parametrization().input_parametrization() == f->parametrization().input_parametrization() )
   {
      if( d->parametrization().output_parametrization() == f->parametrization().output_parametrization() )
      {
         conversion_necessary = false;
         return true;
      }
      else // Ouput parametrizations are different.
         // Output parametrization of  the Function prevails
      {
         parameters p(f->parametrization().dimX(), f->parametrization().dimY(),
                      f->parametrization().input_parametrization(), f->parametrization().output_parametrization());
         converted_data = new vertical_segment(p, d->size());

         //Note that: data_x = dat.head( data->parametrization().dimX() );
         //           data_y = dat.tail( data->parametrization().dimY() );
         vec dat        = vec::Zero(d->parametrization().dimX() + d->parametrization().dimY());
         vec data_x     = vec::Zero(d->parametrization().dimX());
         vec data_y     = vec::Zero(d->parametrization().dimY());
         vec new_data_y = vec::Zero(f->parametrization().dimY());
         vec new_data   = vec::Zero(d->parametrization().dimX() + f->parametrization().dimY());

         for( unsigned int i=0; i < converted_data->size(); i++)
         {
            dat = d->get(i);
            data_x = dat.head( d->parametrization().dimX() );
            data_y = dat.tail( d->parametrization().dimY() );

            params::convert(&data_y[0],
                            d->parametrization().output_parametrization(),
                            d->parametrization().dimY(),
                            f->parametrization().output_parametrization(),
                            f->parametrization().dimY(),
                            &new_data_y[0]);

            new_data.head(d->parametrization().dimX()) = data_x;
            new_data.tail(f->parametrization().dimY()) = new_data_y;

            converted_data->set(i, new_data );
         }

      }
   } else {
      // Input parametrization are different & Output param are the same
      // Converting to the function input parametrization
      if( d->parametrization().output_parametrization() == f->parametrization().output_parametrization() )
      {
         parameters p(f->parametrization().dimX(), f->parametrization().dimY(),
                      f->parametrization().input_parametrization(), f->parametrization().output_parametrization());

         converted_data = new vertical_segment(p, d->size());

         //Note that: data_x = dat.head( data->parametrization().dimX() );
         //           data_y = dat.tail( data->parametrization().dimY() );
         vec dat        = vec::Zero(d->parametrization().dimX() + d->parametrization().dimY());
         vec data_x     = vec::Zero(d->parametrization().dimX());
         vec data_y     = vec::Zero(d->parametrization().dimY());

         vec new_data   = vec::Zero(f->parametrization().dimX() + f->parametrization().dimY());
         vec new_data_y = vec::Zero(f->parametrization().dimY());
         vec new_data_x = vec::Zero(f->parametrization().dimX());

         for(auto i=0; i<converted_data->size(); i++) {
            dat    = d->get(i);
            data_x = dat.head(d->parametrization().dimX());
            data_y = dat.tail(d->parametrization().dimY());

            params::convert(data_x.data(),
                            d->parametrization().input_parametrization(),
                            f->parametrization().input_parametrization(),
                            new_data_x.data());

            new_data.head(f->parametrization().dimX()) = new_data_x;
            new_data.tail(f->parametrization().dimY()) = data_y;

            converted_data->set(i, new_data);
         }
      }
      else //EVERYTHING IS DIFFERENT ! ALL TO FUNCTION PARAMETRIZATIONS
      {
         //       unsigned int const new_size = f->parametrization().dimX() + f->parametrization().dimY();
         // converted_data = new vertical_segment( f->parametrization().dimX(), f->parametrization().dimY(), data->size() );
         // converted_data->setParametrizations( f->parametrization().input_parametrization(), f->parametrization().output_parametrization() );

         // for( unsigned int i=0; i < converted_data->size(); i++)
         // {
         //   vec const dat = data->get(i);
         //   vec const data_x = dat.head( data->parametrization().dimX() );
         //   vec const data_y = dat.tail( data->parametrization().dimY() );
         // }
         NOT_IMPLEMENTED();
      }
   }

   return true;
}


// TODO : Add a new parametriztion to params.h for the COS_TL_TV
//
void computeCosineFactorsFromData( ptr<data> const & data,
      vec & cosine_theta_light,
      vec & cosine_theta_view )
{
   cosine_theta_view.setZero( data->size() );
   cosine_theta_light.setZero( data->size() );

   if( data->parametrization().input_parametrization() == params::CARTESIAN)
   {
      vec data_x = vec::Zero( data->parametrization().dimX() );

      for( unsigned int i=0; i < data->size(); i++)
      {
         data_x = (data->get(i)).head( data->parametrization().dimX() );
         cosine_theta_light[i] = data_x[2];
         cosine_theta_view[i]  = data_x[5];
      }
   }
   else
   {
      NOT_IMPLEMENTED();
   }
}


void  demultiplexData(ptr<data> const & data,
                      Eigen::MatrixXd & o_data_x,
                      Eigen::MatrixXd & o_data_y) {

   vec dat        = vec::Zero( data->parametrization().dimX() + data->parametrization().dimY() );
   vec data_x     = vec::Zero( data->parametrization().dimX() );
   vec data_y     = vec::Zero( data->parametrization().dimY() );

   for(auto i=0; i<data->size(); i++) {
      dat    = data->get(i);
      data_x = dat.head( data->parametrization().dimX() );
      data_y = dat.tail( data->parametrization().dimY() );

      o_data_x.row(i) = data_x;
      o_data_y.row(i) = data_y;
   }
}

void evaluateDataAtData(const ptr<data>& ref,
                        const ptr<data>& dat,
                        Eigen::MatrixXd& ref_y,
                        Eigen::MatrixXd& dat_y) {

   // TODO handle the case of output format conversion
   if(ref->parametrization().output_parametrization()
      != dat->parametrization().output_parametrization()) {
      NOT_IMPLEMENTED();
   }

   // Temp variables
   vec ref_xy = vec::Zero(ref->parametrization().dimX() + ref->parametrization().dimY());
   vec ref_x  = vec::Zero(ref->parametrization().dimX());
   vec dat_x  = vec::Zero(dat->parametrization().dimX());
   vec cart   = vec::Zero(6);

   // Constants
   const auto nY = ref->parametrization().dimY();
   const auto nX = ref->parametrization().dimX();

   // Evaluate the input data at each position of data_x configuration
   for(auto i=0; i<ref->size(); i++) {

      ref_xy = ref->get(i);
      ref_x  = ref_xy.head(nX);

      params::convert(ref_x.data(),
                      ref->parametrization().input_parametrization(),
                      params::CARTESIAN,
                      cart.data());

      // Check if the output configuration is below the hemisphere when
      // converted to cartesian coordinates. Note that this prevent from
      // converting BTDF data.
      if(cart[2] >= 0.0 || cart[5] >= 0.0) {

         // Convert to the query data `dat` parametrization
         params::convert(cart.data(),
                         params::CARTESIAN,
                         dat->parametrization().input_parametrization(),
                         dat_x.data());

         ref_y.row(i) = ref_xy.tail(nY);
         dat_y.row(i) = dat->value(dat_x);
      } else {
         ref_y.row(i).setZero();
         dat_y.row(i).setZero();
      }
   }
}


void fastNormComputation(Eigen::MatrixXd const & o_data_y,
                         Eigen::MatrixXd const & f_y,
                         Eigen::VectorXd & L1,
                         Eigen::VectorXd & L2,
                         Eigen::VectorXd & L3,
                         Eigen::VectorXd & LInf,
                         Eigen::VectorXd & mse,
                         Eigen::VectorXd & rmse )
{
   assert(o_data_y.rows() == f_y.rows());
   assert(o_data_y.cols() == f_y.cols());
   const auto dMatrix = (o_data_y - f_y);

   // Check that the output dimensions match in size for all the
   // vectors and the difference matrix.
   const auto ncols = dMatrix.cols();
   assert(L1.size()   == ncols);
   assert(L2.size()   == ncols);
   assert(L3.size()   == ncols);
   assert(LInf.size() == ncols);

   // Compute the different distance metrics per output dimension
   for(auto i=0; i<ncols; i++) {
      L1(i)   = dMatrix.col(i).lpNorm<1>();
      L2(i)   = dMatrix.col(i).lpNorm<2>();
      L3(i)   = dMatrix.col(i).lpNorm<3>();
      LInf(i) = dMatrix.col(i).lpNorm<Eigen::Infinity>();
   }

   // Compute the RMSE and MSE
   mse = dMatrix.array().square().colwise().sum() / o_data_y.rows();
   rmse = mse.cwiseSqrt();
}

void fastWeightedErrors( Eigen::ArrayXXd const & o_data_y,
      Eigen::ArrayXXd const & f_y,
      Eigen::VectorXd const & weights,
      Eigen::VectorXd & mse,
      Eigen::VectorXd & rmse )
{

   Eigen::MatrixXd  const distance_matrix = (o_data_y - f_y).matrix();
   Eigen::MatrixXd tmp = Eigen::MatrixXd( distance_matrix.rows(), distance_matrix.cols() );

   for( unsigned int i=0; i < f_y.cols() ; i++)
   {
      tmp.col(i) = distance_matrix.col(i).cwiseProduct( weights );
   }

   mse = tmp.array().square().colwise().sum();

   mse /= o_data_y.rows();
   rmse = mse.cwiseSqrt();

}

/*! \package data2stats
 *  \ingroup commands
 *  \brief
 *  This command takes at input an interpolant data and outputs statistics
 *  regarding the distance to another data file.
 *
 *  \details
 *  This command takes at input an interpolant data and outputs statistics
 *  regarding the distance to another data file.
 */
int main(int argc, char* argv[])
{

   arguments args(argc, argv) ;

   if(args.is_defined("help") || args.is_defined("h"))
   {
      std::cout << "Usage: data2stats [options] --input [filename] --ref [filename]" << std::endl ;
      std::cout << std::endl;
      std::cout << "Mandatory arguments:" << std::endl;
      std::cout << "  --ref     [data_filename] : needs to be an ALTA data file" << std::endl;
      std::cout << "  --input   [data_filename] : needs to be an ALTA data file" << std::endl;
      std::cout << std::endl;

      std::cout << "Optional arguments:" << std::endl;
      std::cout << "  --ref-data [plugin_name] : a valid ALTA data plugin name" << std::endl;
      std::cout << "  --in-data  [plugin_name] : a valid ALTA data plugin name" << std::endl;

      return EXIT_SUCCESS;
   }

   std::cout << "<<INFO>> Total System Memory:"
      << plugins_manager::get_system_memory()  / (1024*1024*1024.0) << " Gigabytes "
      << std::endl;


   if(!args.is_defined("ref")) {
      std::cerr << "<<ERROR>> Reference data filename is undefined" << std::endl;
      return EXIT_FAILURE ;
   }

   if(!args.is_defined("input")) {
      std::cerr << "<<ERROR>> Input data filename is undefined" << std::endl ;
      return EXIT_FAILURE ;
   }

   std::cout << "<<INFO>> Loading data ..." << std::endl;
   ptr<data> input;
   //ptr<data> input = plugins_manager::get_data("vertical_segment");

   try
   {
       input = plugins_manager::load_data(args["input"], args["in-data"]);
   }
   CATCH_FILE_IO_ERROR (args["input"]);

   if(!input) {
      std::cout << "<<ERROR>> Could not load data file \'" << args["input"]
                << "\'"  << std::endl;
      return EXIT_FAILURE;
   }

   // Load a function file representing the BRDF
   ptr<data> ref;

   try
   {
       ref = plugins_manager::load_data(args["ref"], args["ref-data"]);
   }
   CATCH_FILE_IO_ERROR (args["ref"]);

   if(!ref) {
      std::cout << "<<ERROR>> Could not load data file \'" << args["ref"]
                << "\'"  << std::endl;
      return EXIT_FAILURE;
   }
   std::cout << "<<INFO>> Starting to compute statistics ... " << std::endl;


   /* Compute the different metrics using the CORE functionalities */
   errors::metrics result;
   errors::compute(input.get(), ref.get(), nullptr, result);

   const auto L1_norm   = result["L1"];
   const auto L2_norm   = result["L2"];
   const auto L3_norm   = result["L3"];
   const auto LInf_norm = result["LInf"];
   const auto mse       = result["MSE"];
   const auto rmse      = result["RMSE"];

   std::cout << "<<INFO>> L1_norm "   << L1_norm   << std::endl
             << "<<INFO>> L2_norm "   << L2_norm   << std::endl
             << "<<INFO>> L3_norm "   << L3_norm   << std::endl
             << "<<INFO>> Linf_norm " << LInf_norm << std::endl
             << "<<INFO>> Mse  "      << mse       << std::endl
             << "<<INFO>> Rmse "      << rmse      << std::endl;


   /* If output is not void we output the different metrics to a file */
   if(args.is_defined("output")) {
      std::string output_filename = args["output"];
      if( output_filename != "") {
         std::ofstream  fwriter(output_filename.c_str() );

         if(fwriter.good()) {
            fwriter << "#FIT2STATS EXPORT" << std::endl;
            fwriter << "#CMD " << args.get_cmd() << std::endl;
            fwriter << "#DIM " << ref->parametrization().dimX() << " "
                               << ref->parametrization().dimY() << std::endl;
            fwriter << std::endl;

            fwriter << "L1 :"   << L1_norm   << std::endl;
            fwriter << "L2 :"   << L2_norm   << std::endl;
            fwriter << "L3 :"   << L3_norm   << std::endl;
            fwriter << "LINF :" << LInf_norm << std::endl;
            fwriter << "MSE :"  << mse       << std::endl;
            fwriter << "RMSE :" << rmse      << std::endl;
            fwriter << std::endl;

         } else {
            std::cerr << "<<ERROR>> Could not save metrics to "
                      << output_filename << std::endl;
            return EXIT_FAILURE;
         }
      }
   }

   return EXIT_SUCCESS;
}
