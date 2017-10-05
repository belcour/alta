/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2016 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include <Eigen/Dense>
//#include <bench/BenchTimer.h>
#include "rational_fitter.h"
#include <QuadProg++.hh>

#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <set>

// #ifdef WIN32
// #define isnan(X) ((X != X))
// #endif

#include <core/common.h>

using namespace alta;

ALTA_DLL_EXPORT fitter* provide_fitter()
{
  return new rational_fitter_quadprog();
}

rational_fitter_quadprog::rational_fitter_quadprog() : _boundary(1.0)
{
}
rational_fitter_quadprog::~rational_fitter_quadprog()
{
}

bool rational_fitter_quadprog::fit_data(const ptr<data>& dat, ptr<function>& fit, const arguments &args)
{
  ptr<rational_function> r = dynamic_pointer_cast<rational_function>(fit) ;
  const ptr<vertical_segment>& d = dynamic_pointer_cast<vertical_segment>(dat) ;
  if(!r || !d
     || d->confidence_interval_kind() != vertical_segment::ASYMMETRICAL_CONFIDENCE_INTERVAL)
  {
    std::cerr << "<<ERROR>> not passing the correct class to the fitter" << std::endl ;
    return false ;
  }

  // XXX: FIT and D may have different values of dimX() and dimY(), but
  // this is fine: we convert values as needed in operator().
  r->setMin(d->min()) ;
  r->setMax(d->max()) ;

  std::cout << "<<INFO>> np in  [" << _min_np << ", " << _max_np
            << "] & nq in [" << _min_nq << ", " << _max_nq << "]" << std::endl ;

  int temp_np = _min_np, temp_nq = _min_nq ;
  while(temp_np <= _max_np || temp_nq <= _max_nq)
  {
    timer time ;
    time.start() ;

    r->setSize(temp_np, temp_nq);
    if(fit_data(d, temp_np, temp_nq, r))
    {
            time.stop() ;
            std::cout << "<<INFO>> got a fit using np = " << temp_np << " & nq =  " << temp_nq << "      " << std::endl ;
            std::cout << "<<INFO>> it took " << time << std::endl ;

      return true ;
    }


    std::cout << "<<INFO>> fit using np = " << temp_np << " & nq =  " << temp_nq << " failed" << std::endl  ;
      time.stop() ;
      std::cout << "<<INFO>> it took " << time << std::endl ;
    std::cout.flush() ;

    if(temp_np == _max_np && temp_nq == _max_nq)
    {
      return false;
    }

    if(temp_np < _max_np)
    {
      ++temp_np ;
    }
    if(temp_nq < _max_nq)
    {
      ++temp_nq ;
    }
  }
  return false ;
}

void rational_fitter_quadprog::set_parameters(const arguments& args)
{
  _max_np = args.get_float("np", 10) ;
  _max_nq = args.get_float("nq", 10) ;
  _min_np = args.get_float("min-np", _max_np) ;
  _min_nq = args.get_float("min-nq", _max_nq) ;

  _max_np = std::max<int>(_max_np, _min_np);
  _max_nq = std::max<int>(_max_nq, _min_nq);

  _boundary = args.get_float("boundary-constraint", 1.0f);

  _scheduling_mode = args.get_string("scheduling-mode","SlidingWindows");
  _scheduling_chunk_size = args.get_int("scheduling-chunk-size",-1);
  _scheduling_grow_factor = args.get_float("scheduling-grow-factor",-1);

  _export_qp = args.is_defined("export-qp");
  _delta = args.get_float("delta", 1) ;
  _add_ls_energy = args.is_defined("add-ls-energy");
}


bool rational_fitter_quadprog::fit_data(const ptr<vertical_segment>& d, int np, int nq, const ptr<rational_function>& r)
{
  // For each output dimension (color channel for BRDFs) perform
  // a separate fit on the y-1D rational function.
  for(int j=0; j<d->parametrization().dimY(); ++j)
  {
    rational_function_1d* rs = r->get(j);
    rs->resize(np, nq);

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
bool rational_fitter_quadprog::fit_data(const ptr<vertical_segment>& d, int np, int nq, int ny, rational_function_1d* r)
{
  using namespace Eigen;
  // Size of the problem
  const int N = np+nq ;
  const int M = d->size() ;

  // Matrices of the problem
  MatrixXd G(N, N) ;    G.setZero();
  VectorXd g(N) ;       g.setZero();
  MatrixXd CI(N, 2*M) ; CI.setZero();
  VectorXd ci(2*M) ;    ci.setZero();
  MatrixXd CE(N, 0) ;   CE.setZero();
  VectorXd ce;

  // Select the size of the result vector to
  // be equal to the dimension of p + q

  MatrixXd Cls(N,M);
  VectorXd x(N);
  // initial guess:
  x.setZero();
  double cost = 0;
  G.setIdentity();
  vec yl, yu, xi;
  // Each constraint (fitting interval or point
  // add another dimension to the constraint
  // matrix
  for(int i=0, k0=0, k1=1; i<M; ++i)
  {
    // Norm of the row vector
    double a0_norm = 0.0 ;
    double a1_norm = 0.0 ;

    xi = d->get(i) ;

    int i0 = i;
    int i1 = i+M;

    // shuffle the constraints to help the solver when working of sub-sets of sequential constraints
    i0 = 2*k0+0;
    i1 = 2*k0+1;

    k0+=32;
    if(k0>=M)
    {
      k0 = k1;
      k1++;
    }

    // A row of the constraint matrix has this
    // form: [p_{0}(x_i), .., p_{np}(x_i), -f(x_i) q_{0}(x_i), .., -f(x_i) q_{nq}(x_i)]
    // For the lower constraint and negated for

    // the upper constrain
    for(int j=0; j<N; ++j)
    {
      // Filling the p part
      if(j<np)
      {
        const double pi = r->p(xi, j);

        CI(j,i0) =  pi;
        CI(j,i1) = -pi;

        Cls(j,i) = pi;
      }
      // Filling the q part
      else
      {
        d->get(i, yl, yu);

        const double qi = r->q(xi, j-np);

        CI(j,i0) = -yu[ny] * qi;
        CI(j,i1) =  yl[ny] * qi;

        Cls(j,i) = -qi*(yu[ny]+yl[ny])/2.0;
      }

      // Update the norm of the row
      a0_norm += CI(j,i0)*CI(j,i0);
      a1_norm += CI(j,i1)*CI(j,i1);
    }

    // Set the c vector, will later be updated using the
    // delta parameter.
    ci[i0] = -sqrt(a0_norm);
    ci[i1] = -sqrt(a1_norm);
  }

  // Add least-square linearized quadratic energy
  if(_add_ls_energy)
    G.noalias() += 1e2 * Cls * Cls.transpose();

#ifdef DEBUG
  std::cout << "CI = [" ;
  for(int j=0; j<2*M; ++j)
  {
    for(int i=0; i<N; ++i)
    {
      std::cout << CI(i,j);
      if(i != N-1) std::cout << ", ";
    }
    if(j != 2*M-1)
      std::cout << ";" << std::endl;
    else
      std::cout << "]" << std::endl ;
  }
#endif

  double delta = _delta;
#if 0
  // Update the ci column with the delta parameter
  // (See Celis et al. 2007 p.12)
  Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::HouseholderQRPreconditioner> svd(CI.transpose());
  const double sigma_m = svd.singularValues()(std::min(2*M, N)-1) ;
  const double sigma_M = svd.singularValues()(0) ;

#ifdef DEBUG
  std::cout << "<<DEBUG>> SVD = [ " ;
  for(int i=0; i<std::min(2*M, N); ++i)
  {
    std::cout << svd.singularValues()(i) << ", " ;
  }
  std::cout << " ]" << std::endl ;
#endif

  delta = sigma_m / sigma_M ;

#ifndef DEBUG
  std::cout << "<<DEBUG>> delta factor: " << sigma_m << " / " << sigma_M << " = " << delta << std::endl ;
#endif

  if(isnan(delta) || (std::abs(delta) == std::numeric_limits<double>::infinity()))
  {
    std::cerr << "<<ERROR>> delta factor is NaN of Inf" << std::endl ;
    return false ;
  }
  else if(delta <= 0.0)
  {
    delta = 1.0 ;
  }
  std::cout << "delta = " << delta << "\n";
#endif
  ci *= delta;

  // Compute the solution
  QuadProgPP::Scheduling scheduling;

  if(_scheduling_mode=="SlidingWindows")
  {
    scheduling.initSlidingWindows(N,2*M);
    if(_scheduling_grow_factor>=1)
      scheduling.grow_factor = _scheduling_grow_factor;
  }
  else if(_scheduling_mode=="WorstSetFirst")
  {
    scheduling.initWorstSetFirst(N,2*M);
  }
  if(_scheduling_chunk_size>0)
    scheduling.initial_chunk_size = _scheduling_chunk_size;

  if(_export_qp)
  {
    // export quadratic objective and all constraints
    std::ofstream Gf("G.txt");
    Gf << std::setprecision(16) << G;
    std::ofstream CIf("CI0.txt");
    CIf << std::setprecision(16) << CI;
    std::ofstream cif("ci1.txt");
    for(int k=0; k<ci.size(); ++k)
      cif << std::setprecision(16) << ci[k] << "\n";
  }

  VectorXi active_set;
//   BenchTimer t;
//   t.reset(); t.start();
  QuadProgPP::init_qp(G);
//   t.stop(); std::cout << "init_qp: " << t.value() << "s\n";
//   t.reset(); t.start();
  cost = QuadProgPP::solve_quadprog_with_guess(G, g, CE, ce, CI, ci, x, scheduling, &active_set);
//   t.stop(); std::cout << "solve: " << t.value() << "s\n";
//   std::cout << "active_set.size(): " << active_set.size() << "\n";

  if(_export_qp)
  {
    // export active constraints
    MatrixXd CI_as(N, active_set.size());
    VectorXd ci_as(active_set.size());
    for(int i=0; i<active_set.size();++i)
    {
      CI_as.col(i) = CI.col(active_set[i]);
      ci_as(i) = ci(active_set[i]);
    }

    std::ofstream CIf_as("CI0_as.txt");
    CIf_as << std::setprecision(16) << CI_as;
    std::ofstream cif_as("ci1_as.txt");
    for(int k=0; k<ci_as.size(); ++k)
      cif_as << std::setprecision(16) << ci_as[k] << "\n";

    // export solution
    std::ofstream Xf("X.txt");
    for(int k=0; k<x.size(); ++k)
      Xf << std::setprecision(16) << x[k] << "\n";
  }

  bool solves_qp = !(cost == std::numeric_limits<double>::infinity());
  for(int i=0; i<np+nq; ++i)
  {
    const double v = x[i];
    solves_qp = solves_qp && !alta::isnan(v) && (v != std::numeric_limits<double>::infinity()) ;
  }

  if(solves_qp)
  {
    // Recopy the vector d
    vec p(np), q(nq);
    double norm = 0.0 ;
    for(int i=0; i<N; ++i)
    {
      const double v = x[i];
      norm += v*v ;
      if(i < np)
      {
        p[i] = v ;
      }
      else
      {
        q[i - np] = v ;
      }
    }

    r->update(p, q);
#ifdef DEBUG
    std::cout << "<<INFO>> got solution " << *r << std::endl ;
#endif
    return norm > 0.0;
  }
  else
  {
    return false;
  }
}
