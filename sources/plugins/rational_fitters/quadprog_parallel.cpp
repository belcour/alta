/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2014 CNRS
   Copyright (C) 2013, 2014, 2016 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

// Include STL
#include <vector>
#include <string>
#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>
#include <string>
#include <list>
#ifdef _OPENMP
#include <omp.h>
#endif

// Interface
#include <core/function.h>
#include <core/rational_function.h>
#include <core/data.h>
#include <core/vertical_segment.h>
#include <core/fitter.h>
#include <core/args.h>
#include <core/plugins_manager.h>

// Eigen
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <QuadProg++.hh>

#include "quadratic_program.h"

using namespace alta;

/*! \brief A vertical segment fitter for rational functions that search for a solution
 *  for a fixed number of coefficient. This plugin can run in parallel with OpenMP and
 *  is using QuadProg++ quadratic solver.
 *  \ingroup plugins
 *  \ingroup fitters
 *
 *  \details
 *  You can find QuadProg++ <a href="http://quadprog.sourceforge.net/">here</a>.
 *  <br />
 *
 *  <h3>Plugin parameters</h3>
 *
 *  We provide the following command line arguments to manipulate this plugin:
 *  <ul>
 * <li><b>--np</b> <em>[int]</em> controls the maximum number of total coefficients
 *   an interpolation should have. By default, this number is 10.</li>
 *   <li><b>--min-np</b> <em>[int]<em></li> controls the starting value for the
 *   number of coefficients for the rational function. by default, this number
 *   is 10.</li>
 *   <li><b>--np-step</b> <em>[int]</em> stepping for the number of coefficients
 *   of the rational function. By default, this number is 1.</li>
 *   <li><b>--nb-cores</b> <em>[int]</em> number of core allocated to perform
 *   the seach. By default, this is equal to the number of processors.</li>
 *   <li><b>--use_delta</b> use the strategy of Pacanowski et al. [2012] to
 *   modify the constraint vector by the condition number of the constraint
 *   matrix. We did not experience any benefit from using it.</li>
 *  </ul>
 */
class rational_fitter_parallel : public fitter
{
   protected: // data

      //! Number of points used when starting the adaptive interpolation
      //! procedure. By default, this value is 100.
      int nb_starting_points;

   public: // methods
      rational_fitter_parallel() : nb_starting_points(100)
      {
      }
      ~rational_fitter_parallel()
      {
      }

      bool fit_data(const ptr<data>& dat, ptr<function>& fit, const arguments &args)
      {
         ptr<rational_function> r = dynamic_pointer_cast<rational_function>(fit) ;
         if(!r)
         {
            std::cerr << "<<ERROR>> not passing the correct function class to the fitter: must be a rational_function" << std::endl ;
            return false ;
         }

         ptr<vertical_segment> d = dynamic_pointer_cast<vertical_segment>(dat) ;
         if(!d
               || d->confidence_interval_kind() != vertical_segment::ASYMMETRICAL_CONFIDENCE_INTERVAL)
         {
            std::cerr << "<<WARNING>> automatic convertion of the data object to vertical_segment," << std::endl;
            std::cerr << "<<WARNING>> we advise you to perform convertion with a separate command." << std::endl;

            size_t elem_size =
               dat->parametrization().dimX() + 3*dat->parametrization().dimY();
            double* content = new double[dat->size() * elem_size];

            for(int i=0; i<dat->size(); ++i)
            {
               const vec x = dat->get(i);

               for(int k=0; k<x.size(); ++k) {
                  content[i + k] = x[k];
               }
               for(int k=0; k<dat->parametrization().dimY(); ++k) {
                  content[i + k + dat->parametrization().dimX() + dat->parametrization().dimY()] =
                     (1.0 - args.get_float("dt", 0.1)) * x[k + dat->parametrization().dimX()];
               }
               for(int k=0; k<dat->parametrization().dimY(); ++k) {
                  content[i + k + dat->parametrization().dimX() + 2*dat->parametrization().dimY()] =
                     (1.0 + args.get_float("dt", 0.1)) * x[k + dat->parametrization().dimX()];
               }
            }

            ptr<vertical_segment> vs(new vertical_segment(dat->parametrization(),
                     dat->size(),
                     std::shared_ptr<double>(content)));

            d = vs;
         }

         // XXX: FIT and D may have different values of dimX() and dimY(), but
         // this is fine: we convert values as needed in operator().
         r->setMin(d->min());
         r->setMax(d->max());

         const int _min_np = args.get_int("min-np", 10);
         const int _max_np = args.get_int("np", _min_np);
         std::cout << "<<INFO>> N in  [" << _min_np << ", " << _max_np << "]"  << std::endl ;

         const int nb_starting_points = args.get_int("nb-starting-points", 100);
         std::cout << "<<INFO>> number of data point used in start: " << nb_starting_points << std::endl;

         const int  step      = args.get_int("np-step", 1);
         const bool use_delta = args.is_defined("use_delta");

         for(int i=_min_np; i<=_max_np; i+=step)
         {
            std::cout << "<<INFO>> fit using np+nq = " << i << std::endl ;
            std::cout.flush() ;
            timer time ;
            time.start() ;

#ifdef _OPENMP
            const int nb_cores = args.get_int("nb-cores", omp_get_num_procs());
#ifdef DEBUG
            std::cout << "<<DEBUG>> will use " << nb_cores << " threads to compute the quadratic programs" << std::endl ;
#endif

            omp_set_num_threads(nb_cores) ;
#endif

            double min_delta   = std::numeric_limits<double>::max();
            double min_l2_dist = std::numeric_limits<double>::max();
            double mean_delta = 0.0;
            int nb_sol_found  = 0;
            int nb_sol_tested = 0;

#pragma omp parallel for shared(r, args, nb_sol_found, nb_sol_tested, min_delta, mean_delta), schedule(dynamic,1)
            for(int j=1; j<i; ++j)
            {
               // Compute the number of coefficients in the numerator and in the denominator
               // from the current number of coefficients i and the current index in the
               // loop j.
               int temp_np = i - j;
               int temp_nq = j;

               //vec p(temp_np*r->dimY()), q(temp_nq*r->dimY());

               // Allocate a rational function and set it to the correct size, dimensions
               // and parametrizations.
               ptr<rational_function> rk(NULL);
#pragma omp critical (args)
               {
                  rk = dynamic_pointer_cast<rational_function>(
                        ptr<function>(plugins_manager::get_function(args,
                              r->parametrization())));
               }
               if(!rk)
               {
                  std::cerr << "<<ERROR>> unable to obtain a rational function from the plugins manager" << std::endl;
                  throw;
               }

               rk->setMin(r->min()) ;
               rk->setMax(r->max()) ;

               // Set the rational function size
               rk->setSize(temp_np, temp_nq);

               double delta = 1.0;
               double linf_dist, l2_dist;
               bool is_fitted = fit_data(d, temp_np, temp_nq, rk, args, delta, linf_dist, l2_dist);
               if(is_fitted)
               {
#pragma omp critical (r)
                  {
                     ++nb_sol_found ;
                     mean_delta += delta ;

                     std::cout << "<<INFO>> found a solution with np=" << temp_np
                        << ", nq = " << temp_nq << std::endl;
                     std::cout << "<<INFO>> Linf error = " << linf_dist << std::endl;
                     std::cout << "<<INFO>> L2   error = " << l2_dist << std::endl;
                     std::cout << "<<INFO>>      delta = " << delta << std::endl;
                     std::cout << std::endl;

                     // Get the solution with the minimum delta or the minimum L2 distance,
                     // and update the main rational function r.
                     if((use_delta && delta < min_delta) || (!use_delta && l2_dist < min_l2_dist))
                     {
                        min_delta   = delta ;
                        min_l2_dist = l2_dist ;
                        r->setSize(temp_np, temp_nq);
                        r->update(rk);
                     }
                  }
               }

#pragma omp critical (nb_sol_tested)
               {
                  // Update the solution
                  nb_sol_tested++;

                  std::cout << "<<DEBUG>> nb solutions tested: " << nb_sol_tested << " / " << i << "\r";
                  std::cout.flush();
               }
            }

            if(min_delta < std::numeric_limits<double>::max())
            {
               std::cout << "<<INFO>> mean delta = " << mean_delta/nb_sol_found << std::endl;
               std::cout << "<<INFO>>  min delta = " << min_delta << std::endl;
               std::cout << *(r.get()) << std::endl;

               time.stop();
               std::cout << "<<INFO>> got a fit using N = " << i << std::endl ;
               std::cout << "<<INFO>> it took " << time << std::endl ;
               std::cout << "<<INFO>> I got " << nb_sol_found << " solutions to the QP" << std::endl ;
               return true ;
            }
         }

         return false ;
      }

      void set_parameters(const arguments&)
      {
      }


      bool fit_data(const ptr<vertical_segment>& d, int np, int nq,
            const ptr<rational_function>& r, const arguments &args,
            double& delta, double& linf_dist, double& l2_dist)
      {
         // Fit the different output dimension independantly
         for(int j=0; j<d->parametrization().dimY(); ++j)
         {
            vec p(np), q(nq);
            rational_function_1d* rf = r->get(j);
            rf->resize(np, nq);

            if(!fit_data(d, np, nq, j, rf, args, p, q, delta))
            {
               return false ;
            }

            rf->update(p, q);
         }

         linf_dist = r->Linf_distance(dynamic_pointer_cast<data>(d));
         l2_dist   = r->L2_distance(dynamic_pointer_cast<data>(d));

         return true ;
      }

      // dat is the data object, it contains all the points to fit
      // np and nq are the degree of the RP to fit to the data
      // y is the dimension to fit on the y-data (e.g. R, G or B for RGB signals)
      // the function returns a rational BRDF function and a boolean
      bool fit_data(const ptr<vertical_segment>& d, int np, int nq, int ny,
            rational_function_1d* r, const arguments& args,
            vec& p, vec& q, double& delta)
      {
         const int m = d->size(); // 2*m = number of constraints
         const int n = np+nq;     // n = np+nq

         quadratic_program qp(np, nq, args.is_defined("use_delta"));

         // Starting with only a nb_starting_points vertical segments
         std::list<unsigned int> training_set;
         const int di = std::max((m-1) / (nb_starting_points-1), 1);
         for(int i=0; i<m; ++i)
         {
            if(i % di == 0)
            {
               // Create two vector of constraints
               vec c1(n), c2(n);
               get_constraint(i, np, nq, ny, d, r, c1, c2);

               qp.add_constraints(c1);
               qp.add_constraints(c2);

            }
            else
            {
               training_set.push_back(i);
            }
         }
         qp.set_training_set(training_set);

         do
         {
#ifdef _OPENMP
#ifdef DEBUG
            std::cout << "<<DEBUG>> thread " << omp_get_thread_num() << ", number of intervals tested = " << qp.nb_constraints()/2 << std::endl ;
#endif
#endif
            Eigen::VectorXd x(n);
            bool solves_qp = qp.solve_program(x, delta, p, q);
            r->update(p, q);

            if(solves_qp)
            {
               if(qp.test_constraints(ny, r, d))
               {
#ifdef DEBUG
                  std::cout << "<<INFO>> got solution " << *r << std::endl ;
#endif
                  return true;
               }
            }
            else
            {
#ifdef DEBUG
               std::cout << "<<DEBUG>> not enough coefficients" << std::endl;
#endif
               return false;
            }
         } while(qp.nb_constraints() < 2*m);

         return false;
      }

      void get_constraint(int i, int np, int nq, int ny,
            const ptr<vertical_segment>& data,
            const rational_function_1d* func,
            vec& cu, vec& cl)
      {
         const vec xi = data->get(i) ;
         cu.resize(np+nq);
         cl.resize(np+nq);

         // Create two vectors of constraints
         for(int j=0; j<np+nq; ++j)
         {
            // Filling the p part
            if(j<np)
            {
               const double pi = func->p(xi, j) ;
               cu[j] =  pi ;
               cl[j] = -pi ;
            }
            // Filling the q part
            else
            {
               vec yl, yu ;
               data->get(i, yl, yu) ;
               const double qi = func->q(xi, j-np) ;

               cu[j] = -yu[ny] * qi ;
               cl[j] =  yl[ny] * qi ;
            }
         }
      }
};

ALTA_DLL_EXPORT fitter* provide_fitter()
{
   return new rational_fitter_parallel();
}

