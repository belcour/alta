/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2013, 2014, 2016 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */


// Include STL
#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

// ALTA Interface
#include <core/common.h>
#include <core/function.h>
#include <core/data.h>
#include <core/fitter.h>
#include <core/args.h>
#include <core/rational_function.h>
#include <core/vertical_segment.h>

// Eigen includes
#include <Eigen/Dense>
#include <Eigen/SVD>

using namespace alta;

/*! \brief A least square fitter for rational function
 *  \ingroup plugins
 *  \ingroup fitters
 *  \todo :  WRITE MORE
 */
class rational_fitter_leastsquare : public fitter
{
   protected: // data
		
      int _np, _nq ;
		int _max_iter;

   public: // methods

      rational_fitter_leastsquare() : _np(10), _nq(10), _max_iter(1)
      {
      }
      ~rational_fitter_leastsquare() 
      {
      }

      bool fit_data(const ptr<data>& dat, ptr<function>& fit, const arguments &args)
      {
         ptr<rational_function> r = dynamic_pointer_cast<rational_function>(fit) ;
         const ptr<vertical_segment> d = dynamic_pointer_cast<vertical_segment>(dat) ;
         if(!r || !d
         || d->confidence_interval_kind() != vertical_segment::ASYMMETRICAL_CONFIDENCE_INTERVAL)
         {
            std::cerr << "<<ERROR>> not passing the correct class to the fitter" << std::endl ;
            return false ;
         }

      // XXX: FIT and D may have different values of dimX() and dimY(), but
      // this is fine: we convert values as needed in operator().
         r->setMin(d->min());
         r->setMax(d->max());

         _np = args.get_int("np", _np);
         _nq = args.get_int("nq", _nq);
         r->setSize(_np, _nq);

         std::cout << "<<INFO>> Fitting using np = " << _np <<  " & nq = " << _nq  << std::endl ;
         

         timer time ;
         time.start() ;


         if(fit_data(d, _np, _nq, r))
         {
            time.stop();
            std::cout << "<<INFO>> got a fit, it took " << time << std::endl ;
            return true ;
         }

         std::cout << "<<INFO>> fit failed\r"  ;
         std::cout.flush() ;

         return false ;
      }

      void set_parameters(const arguments& args)
      {
         _np = args.get_int("np", 10) ;
         _nq = args.get_int("nq", 10) ;
         _max_iter = args.get_int("max-iter", 1) ;
      }
            
      bool fit_data(const ptr<vertical_segment>& d, int np, int nq, const ptr<rational_function>& r) 
      {
         // For each output dimension (color channel for BRDFs) perform
         // a separate fit on the y-1D rational function.
         for(int j=0; j<d->parametrization().dimY(); ++j)
         {
            rational_function_1d* rs = r->get(j);

            if(!fit_data(d, np, nq, j, rs))
            {
                  return false ;
            }
         }

         return true ;
      }

      // dat is the data object, it contains all the points to fit
      // np and nq are the degree of the RP to fit to the data
      // y is the dimension to fit on the y-data (e.g. R, G or B for RGB signals)
      // the function return a ration BRDF function and a boolean
      bool fit_data(const ptr<vertical_segment>& d, int np, int nq, int ny, rational_function_1d* r)
      {
      using namespace Eigen;
      
      double scale = 1;
      
      MatrixXd D(d->size(), np+nq);
      VectorXd Y(d->size());
      for(int i=0; i<d->size(); ++i) 
      {
         const double y = d->get(i)[d->parametrization().dimX() + ny];
         Y[i] = y;
         vec x = d->get(i);
         VectorXd::Map(&x[0], 1) /= scale;
         // A row of the constraint matrix has this 
         // form: [p_{0}(x_i), .., p_{np}(x_i), q_{0}(x_i), .., q_{nq}(x_i)]
         for(int j=0; j<np+nq; ++j)
         {
            if(j<np)
            {
            const double pi = r->p(x, j) ;
            D(i, j) =  pi ;
            }
            else
            {
            const double qi = r->q(x, j-np) ;
            D(i,j) = qi ;
            }
         }
      }

      VectorXd pq(np+nq);
      
      // initialize with a constant numerator:
      pq.head(np).setZero();
      pq(0) = 1;
      
      // Alternate between fixed numerator and fixed denominator
      // By default we iterate only once (usually successive iteration introduced more errors)
      for(int iter=0; iter<_max_iter; ++iter)
      {
         // Step 1 fit 1/f_i using 1/q only
         {
            // Theoretically best weighting scheme to approcimate the true LS problem, which is: sum_i (1/q(x_i) - f_i)^2
            MatrixXd A = Y.cwiseAbs2().asDiagonal() * D.rightCols(nq);
            VectorXd b = Y.asDiagonal() * (D.leftCols(np) * pq.head(np));
            VectorXd x = A.colPivHouseholderQr().solve(b);
            
            // Sometimes, this scheme works better:
      //       MatrixXd A = Y.asDiagonal() * D.rightCols(nq);
      //       VectorXd b = D.leftCols(np) * pq.head(np);
      //       VectorXd x = A.colPivHouseholderQr().solve(b);
            
            pq.tail(nq) = x;
            
            // Statistics on the problem:
      //       VectorXd sv = A.jacobiSvd().singularValues();
      //       std::cout << "Singular values: " << sv.transpose() << std::endl;
      //       std::cout << "Rcond: " << sv(0)/sv(nq-1) << std::endl;
      //       std::cout << "<<INFO>> LS error " << (A*x-b).norm()/b.norm() << std::endl;
         }
         
         VectorXd res = ((D.leftCols(np) * pq.head(np)).array() / (D.rightCols(nq) * pq.tail(nq)).array() - Y.array());
      /*
      std::cout << "<<INFO>> Real LS norm (1): "
            << res.norm() / Y.norm()
            << " ; inf " << res.lpNorm<Infinity>() / Y.lpNorm<Infinity>()
            << " ; L1  " << res.lpNorm<1>() / Y.lpNorm<1>() << std::endl;
      */
         // Step 2 fit f_i using p/q with q fix
         {
            // This scheme gives more weights to small values: (not good)
      //     VectorXd V = D.rightCols(nq) * pq.tail(nq);
      //     MatrixXd A = Y.asDiagonal().inverse() * (V.asDiagonal().inverse() * D.leftCols(np));
      //     VectorXd b = VectorXd::Ones(Y.size());
      //     VectorXd x = A.colPivHouseholderQr().solve(b);
            
            VectorXd Q = D.rightCols(nq) * pq.tail(nq);
            MatrixXd A = Q.asDiagonal().inverse() * D.leftCols(np);
            VectorXd x = A.colPivHouseholderQr().solve(Y);
            
            pq.head(np) = x;

            // Statistics on the problem:
      //       VectorXd sv = A.jacobiSvd().singularValues();
      //       std::cout << "Singular values: " << sv.transpose() << std::endl;
      //       std::cout << "Rcond: " << sv(0)/sv(np-1) << std::endl;
      //       std::cout << "<<INFO>> LS error " << (A*x-Y).norm()/Y.norm() << std::endl;
         }
      } // iterations
      
      vec p(np), q(nq);
      Eigen::VectorXd::Map(&p[0], np) = pq.head(np);
      Eigen::VectorXd::Map(&q[0], nq) = pq.tail(nq);
      
      double s = 1;
      for(int i=0; i<np; ++i)
      {
         p[i] /= s;
         s *= scale;
      }
      s = 1;
      for(int i=0; i<nq; ++i)
      {
         q[i] /= s;
         s *= scale;
      }
      
      // Evaluate true LS error: sum_i (p(x_i)/q(x_i) - f_i)^2
      VectorXd res = ((D.leftCols(np) * pq.head(np)).array() / (D.rightCols(nq) * pq.tail(nq)).array() - Y.array());
      std::cout << "<<INFO>> L_2 "
         << res.norm() / Y.norm()
         << " ; L_inf " << res.lpNorm<Infinity>() / Y.lpNorm<Infinity>()
         << " ; L_1  " << res.lpNorm<1>() / Y.lpNorm<1>() << std::endl;

      r->update(p, q) ;
      std::cout << "<<INFO>> got solution " << *r << std::endl ;
      return true;

      }
};

ALTA_DLL_EXPORT fitter* provide_fitter()
{
    return new rational_fitter_leastsquare();
}
