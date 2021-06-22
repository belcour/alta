/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#pragma once

#include <Eigen/Dense>
#include <QuadProg++.hh>

#include <core/rational_function.h>
#include <core/vertical_segment.h>

using namespace alta;

class quadratic_program
{
   public:
      //! \brief Constructor need to specify the number of coefficients
      quadratic_program(int np, int nq, bool compute_delta = false) :
        _np(np), _nq(nq), _compute_delta(compute_delta), CI(_np+_nq, 0)
      { }

      //! \brief Remove the already defined constraints
      void clear_constraints()
      {
         CI.resize(_np+_nq, 0);
      }

      //! \brief Add a constraint by specifying the vector
      void add_constraints(const vec& c)
      {
         const int m = CI.rows();
         const int n = CI.cols();

         if(n > 0)
         {
            // Resize matrix CI
            CI.conservativeResize(m, n+1);
            CI.col(n) = c;
         }
         else
         {
            CI = c;
         }
      }

      //! \brief Provide the number of constraints
      int nb_constraints() const
      {
         return CI.cols();
      }

      //! Set the indices of the remaining data
      void set_training_set(const std::list<unsigned int>& ts)
      {
         this->training_set = ts;
      }

      //! \brief Solves the quadratic program and update the p and
      //! q vector if necessary.
      inline bool solve_program(Eigen::VectorXd& x, double& delta, vec& p, vec& q)
      {
         bool solves_qp = solve_program(x, delta) ;

         if(solves_qp)
         {
            double norm = 0.0;
            for(int i=0; i<_np+_nq; ++i)
            {
               const double v = x[i];
               norm += v*v ;
               if(i < _np)
               {
                  p[i] = v ;
               }
               else
               {
                  q[i-_np] = v ;
               }
            }
            return norm > 0.0;
         }
         else
         {
            return false ;
         }
      }

      //! \brief Solves the quadratic program
      inline bool solve_program(Eigen::VectorXd& v, double& delta)
      {
         const int m = CI.rows();
         const int n = CI.cols();

         Eigen::MatrixXd G (m, m) ;
         Eigen::VectorXd g (m) ;
         Eigen::VectorXd ci(n) ;
         Eigen::MatrixXd CE(m, 0) ;
         Eigen::VectorXd ce ;

         if(_compute_delta)
         {
            // Update the ci column with the delta parameter
            // (See Celis et al. 2007 p.12)
            Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::HouseholderQRPreconditioner> svd(CI.transpose());
            const double sigma_m = svd.singularValues()(std::min(m, n)-1) ;
            const double sigma_M = svd.singularValues()(0) ;
            delta = sigma_M / sigma_m ;
         }

         // Select the size of the result vector to
         // be equal to the dimension of p + q
         G.setIdentity();

         // Each constraint (fitting interval or point
         // add another dimension to the constraint
         // matrix
         for(int i=0; i<n; ++i)
         {
            // Norm of the row vector
            double norm = 0.0 ;

            for(int j=0; j<m; ++j)
            {
               norm += CI(j,i)*CI(j,i);
            }

            // Set the c vector, will later be updated using the
            // delta parameter.
            ci[i] = -delta * sqrt(norm) ;
         }

         // Compute the solution
         const double cost = QuadProgPP::solve_quadprog(G, g, CE, ce, CI, ci, v);

         bool solves_qp = !(cost == std::numeric_limits<double>::infinity());
         return solves_qp;
      }

#define PACANOWSKI2012

        //! \brief Test all the constraints of the data.
        //! Add the sample that is farest away from the function.
      bool test_constraints(int ny, const rational_function_1d* r, const ptr<vertical_segment>& data)
      {
#ifdef PACANOWSKI2012
            int nb_failed = 0;
            double max_dev = 0.0; // Maximum absolute distance of the current
                                  // solution to the data.
            std::list<unsigned int>::iterator max_ind;
            vec cu, cl;

            std::list<unsigned int>::iterator it;
            for(it = training_set.begin(); it != training_set.end(); it++)
            {
               vec x, yl, yu;
               data->get(*it, x, yl, yu);

               vec y = r->value(x);
               bool fail_upper = y[0] > yu[ny];
               bool fail_lower = y[0] < yl[ny];
               if(fail_lower || fail_upper)
               {
                  const double dev = std::abs(0.5*(yu[ny]+yl[ny]) - y[0]);

                  nb_failed++;

                  if(max_dev < dev)
                  {
                     get_constraint(x, yl, yu, ny, r, cu, cl);
                     max_dev = dev;
                     max_ind = it;
                  }
               }
            }
#ifdef DEBUG
            std::cout << "<<TRACE>> " << nb_failed << " constraints where not satified." << std::endl;
            std::cout << "<<TRACE>> an interval failed the test with distance = " << max_dev << std::endl;
#endif

            if(nb_failed > 0)
            {
                add_constraints(cu);
                add_constraints(cl);
                training_set.erase(max_ind);
#ifdef DEBUG
                std::cout << "<<DEBUG>> number of remaining training elements: " << training_set.size() << std::endl;
#endif

                return false;
            }
            else
            {
                return true;
            }
#else
            int n = next_unmatching_constraint(0, ny, r, data);
            if(n < data->size())
            {
                vec x, yl, yu;
                data->get(n, x, yl, yu);

                vec cu, cl;
                get_constraint(x, yl, yu, ny, r, cu, cl);

                add_constraints(cu);
                add_constraints(cl);

                return false;
            }
            else
            {
                return true;
            }
#endif
      }

      //! \brief Generate two constraint vectors from a vertical segment and a
      //! ration function type.
      inline void get_constraint(const vec& xi, const vec& yl, const vec& yu, int ny,
                                 const rational_function_1d* func,
                                 vec& cu, vec& cl);

      //! \brief Give the next position in the data that is not satisfied.
      //! This method works only for a single color channel ny !
      static int next_unmatching_constraint(int i, int ny, const rational_function_1d* r,
                                            const vertical_segment* data);

   protected:
      int _np, _nq;
      bool _compute_delta;
      Eigen::MatrixXd CI;

      //! Contains the indices of the vertical segment unused during the
      //! rational interpolation.
      std::list<unsigned int> training_set;
};


inline void quadratic_program::get_constraint(const vec& xi, const vec& yl, const vec& yu,
                                              int ny, const rational_function_1d* func,
                                              vec& cu, vec& cl)
{
   cu.resize(_np+_nq);
   cl.resize(_np+_nq);

   // Create two vector of constraints
   for(int j=0; j<_np+_nq; ++j)
   {
      // Filling the p part
      if(j<_np)
      {
         const double pi = func->p(xi, j) ;
         cu[j] =  pi ;
         cl[j] = -pi ;

      }
      // Filling the q part
      else
      {
         const double qi = func->q(xi, j-_np) ;

         cu[j] = -yu[ny] * qi ;
         cl[j] =  yl[ny] * qi ;
      }
   }
}

int quadratic_program::next_unmatching_constraint(int i, int ny, const rational_function_1d* r,
                                                  const vertical_segment* data)
{
   for(int n=i; n<data->size(); ++n)
   {
      vec x, yl, yu;
      data->get(n, x, yl, yu);

      vec y = r->value(x);
        if(y[0] < yl[ny] || y[0] > yu[ny])
      {
         return n;
      }
   }
   return data->size();
}
