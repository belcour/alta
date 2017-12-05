/*
  LICENSE

  This file is part of QuadProg++: a C++ library implementing
  the algorithm of Goldfarb and Idnani for the solution of a (convex)
  Quadratic Programming problem by means of an active-set dual method.

  Copyright (C) 2007-2009 Luca Di Gaspero <l.digaspero@uniud.it>
  Copyright (C) 2009 Eric Moyer.
  Copyright (C) 2014-2015 Gael Guennebaud <gael.guennebaud@inria.fr>

  QuadProg++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  QuadProg++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with QuadProg++. If not, see <http://www.gnu.org/licenses/>.

*/

// #undef NDEBUG
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <Eigen/Cholesky>
#include <Eigen/Jacobi>
#include "QuadProg++.hh"

#ifdef DEBUG
#define TRACE_SOLVER_(CODE) std::cout << __LINE__ << ": "; CODE; std::cout << std::endl
// #define TRACE_SOLVER(CODE) TRACE_SOLVER_(CODE)
#define TRACE_SOLVER(CODE) do {} while(0)
#else
#define TRACE_SOLVER_(CODE)
// #define TRACE_SOLVER(CODE) TRACE_SOLVER_(CODE)
#define TRACE_SOLVER(CODE)
#endif
// #define QP_SKIP_DEGENERATE

namespace QuadProgPP{

// Utility functions for updating some data needed by the solution method
EIGEN_DONT_INLINE void compute_step_direction(Ref<VectorXd> z, Ref<VectorXd> r, Ref<const VectorXd> np, Ref<const MatrixXd> J, Ref<const MatrixXd> R, int iq, Ref<VectorXd> d);
EIGEN_DONT_INLINE bool add_constraint(Ref<MatrixXd> R, Ref<MatrixXd> J, Ref<VectorXd> d, int& iq, double& rnorm);
EIGEN_DONT_INLINE void delete_constraint(Ref<MatrixXd> R, Ref<MatrixXd> J, VectorXi& A, Ref<VectorXd> u, int n, int p, int& iq, int l);

// Compute sqrt(a^2 + b^2) without underflow or overflow
double hypot(double a, double b);

void init_qp(Ref<MatrixXd> G)
{
  assert(G.rows()==G.cols());

  internal::llt_inplace<double, Lower>::blocked(G);
}

// The Solving function, implementing the Goldfarb-Idnani method
double solve_quadprog(Ref<MatrixXd> G, Ref<const VectorXd> g0,
                      Ref<const MatrixXd> CE, Ref<const VectorXd> ce0,
                      Ref<const MatrixXd> CI, Ref<const VectorXd> ci0,
                      Ref<VectorXd> x)
{
  std::ostringstream msg;
  {
    // Ensure that the dimensions of the matrices and vectors can be
    // safely converted from unsigned int into to int without overflow.
    unsigned mx = std::numeric_limits<int>::max();
    if(G.cols() >= mx   || G.rows() >= mx   ||
       CE.rows() >= mx  || CE.cols() >= mx  ||
       CI.rows() >= mx  || CI.cols() >= mx  ||
       ci0.size() >= mx || ce0.size() >= mx || g0.size() >= mx)
    {
      msg << "The dimensions of one of the input matrices or vectors were "
          << "too large." << std::endl
          << "The maximum allowable size for inputs to solve_quadprog is:"
          << mx << std::endl;
      throw std::logic_error(msg.str());
    }
  }
  int n = G.cols(), p = CE.cols(), m = CI.cols();
  if ((int)G.rows() != n)
  {
    msg << "The matrix G is not a square matrix (" << G.rows() << " x " << G.cols() << ")";
    throw std::logic_error(msg.str());
  }
  if ((int)CE.rows() != n)
  {
    msg << "The matrix CE is incompatible (incorrect number of rows " << CE.rows() << " , expecting " << n << ")";
    throw std::logic_error(msg.str());
  }
  if ((int)ce0.size() != p)
  {
    msg << "The vector ce0 is incompatible (incorrect dimension " << ce0.size() << ", expecting " << p << ")";
    throw std::logic_error(msg.str());
  }
  if ((int)CI.rows() != n)
  {
    msg << "The matrix CI is incompatible (incorrect number of rows " << CI.rows() << " , expecting " << n << ")";
    throw std::logic_error(msg.str());
  }
  if ((int)ci0.size() != m)
  {
    msg << "The vector ci0 is incompatible (incorrect dimension " << ci0.size() << ", expecting " << m << ")";
    throw std::logic_error(msg.str());
  }

  // decompose the matrix G in the form L^T L
  init_qp(G);

  // Find the unconstrained minimizer of the quadratic form 0.5 * x G x + g0 x
  // this is a feasible point in the dual space
  // x = G^-1 * g0
  x = G.triangularView<Lower>().solve(g0);
  x = G.triangularView<Lower>().adjoint().solve(x);
  x = -x;

  return solve_quadprog_with_guess(G, g0, CE, ce0, CI, ci0, x);
}

template <typename ValuesType, typename MatrixType, typename VectorType>
void partial_sort(ValuesType &values, MatrixType &M, VectorType &m, int start, int ncut)
{
  typedef typename ValuesType::RealScalar RealScalar;
  typedef typename ValuesType::Index Index;
  using std::swap;
  Index mid;
  Index n = values.size()-start;
  Index first, last ;

  ncut--; // to fit the zero-based indices
  first = start;
  last = n-1;
  if (ncut < first || ncut > last ) return;

  do {
    mid = first;
    RealScalar abskey = values(mid);
    for (Index j = first + 1; j <= last; j++) {
      if ( values(j) < abskey) {
        ++mid;
        swap(values(mid), values(j));
        swap(m(mid), m(j));
        M.col(mid).swap(M.col(j));
      }
    }
    /* Interchange for the pivot element */
    swap(values(mid), values(first));
    swap(m(mid), m(first));
    M.col(mid).swap(M.col(first));

    if (mid > ncut) last = mid - 1;
    else if (mid < ncut ) first = mid + 1;
  } while (mid != ncut );
}

VectorXd g_u;

double solve_quadprog_with_guess(Ref<const MatrixXd> L, Ref<const VectorXd> g0,
                                 Ref<const MatrixXd> CE, Ref<const VectorXd> ce0,
                                 Ref<const MatrixXd> CI_, Ref<const VectorXd> ci0_,
                                 Ref<VectorXd> x,
                                 Scheduling scheduling,
                                 VectorXi *active_set)
{
  MatrixXd CI(CI_);
  VectorXd ci0(ci0_);
  int n = L.cols(), p = CE.cols(), m = CI.cols();

  int ip; // this is the index of the constraint to be added to the active set
  MatrixXd R(n, n), J(n, n);
  VectorXd s(m + p), z(n), r(m + p), d(n), np(n), u(m + p);
  double f_value, psi, c1, c2, ss, R_norm;
  double inf;
  if (std::numeric_limits<double>::has_infinity)
      inf = std::numeric_limits<double>::infinity();
  else
      inf = 1.0E300;
  double t, t1, t2; // t is the step length, which is the minimum of the partial step length t1
                    // and the full step length t2
  VectorXi A(m + p), iai(m + p);
  int q, iq;

#ifdef QP_SKIP_DEGENERATE
  std::vector<bool> iaexcl(m + p, true);
  VectorXd x_old(n), u_old(m + p);
  VectorXi A_old(m + p);
#endif

  // p is the number of equality constraints
  // m is the number of inequality constraints
  q = 0;  // size of the active set A (containing the indices of the active constraints)

  // Preprocessing phase
  c1 = L.diagonal().squaredNorm(); // == trace(LL^T) == trace(G)

  // initialize the matrix R
  d.setZero();
  R.setZero();
  R_norm = 1.0; /* this variable will hold the norm of the matrix R */

  // compute the inverse of the factorized matrix G^-1, this is the initial value for H
  c2 = 0.0; // trace of J
  for (int i = 0; i < n; i++)
  {
    d[i] = 1.0;
    J.col(i) = L.transpose().triangularView<Upper>().solve(d);
    c2 += J(i,i);
    d[i] = 0.0;
  }
  // c1 * c2 is an estimate for cond(G)

  // and compute the current solution value
  f_value = 0.5 * g0.dot(x);

  TRACE_SOLVER( std::cout << "Unconstrained solution: " << f_value << std::endl );
  TRACE_SOLVER( std::cout << "x: " << x.transpose() << std::endl );

  // Add equality constraints to the working set A
  iq = 0;
  for (int i = 0; i < p; i++)
  {
    np = CE.col(i);
    compute_step_direction(z, r, np, J, R, iq, d);

    // compute full step length t2:
    // i.e., the minimum step in primal space s.t. the constraint becomes feasible
    t2 = 0.0;
    if(z.squaredNorm() > std::numeric_limits<double>::epsilon()) // i.e. z != 0
      t2 = (-(np.dot(x)) - ce0[i]) / z.dot(np);

    x += t2 * z;

    // set u = u+
    u[iq] = t2;
    u.head(iq) -= t2 * r.head(iq);

    // compute the new solution value
    f_value += 0.5 * (t2 * t2) * z.dot(np);
    A[i] = -i - 1;

    if (!add_constraint(R, J, d, iq, R_norm))
    {
      // Equality constraints are linearly dependent
      throw std::runtime_error("Constraints are linearly dependent");
      return f_value;
    }
  }

  // set iai = K \ A
  for (int i = 0; i < m; i++) iai[i] = i;

  int chunk_start=0, chunk_size = scheduling.initial_chunk_size;
  if(scheduling.mode==Scheduling::WorstFirst)
    chunk_size = m;

//   if(active_set)
//   {
//     // restore active set (assume initial guess satisfies given active set)
//     for(int i=0; i<active_set->size(); ++i)
//     {
//       int j = active_set->coeff(i);
//       iai[j] = -1;
//       A[i] = j;
//       //u[i] = r[i] = d[i] = z[i] = 0;
//       add_constraint(R, J, d, iq, R_norm);
//     }
//     if(active_set->size()>0)
//       u.head(iq+1) = g_u.head(iq+1);
//     //u[iq] = r[iq] = d[iq] = z[iq] = 0;
//   }

  bool done_with_active_set = active_set==0 || active_set->size()==0;
  bool done_with_current_subset = false;
  long flops0 = 0, flops1 = 0;
  int iter = 0;
  int nb_chunks = (m+chunk_size-1)/chunk_size;
  int nb_satisfied_chunks = 0;
  while(iter<=10*m) // make sure we do not run into an infinite loop
  {
    TRACE_SOLVER( std::cout << "x: " << x.transpose() << std::endl );
    // step 1: choose a violated constraint

    // FIXME: should be removed:
    for (int i = p; i < iq; i++) {
      if(iai[A[i]]!=-1) { std::cerr << "A <-> iai missmatch\n"; abort();}
      //iai[A[i]] = -1;
    }

    // compute s[x] = ci^T * x + ci0 for all elements of K \ A
    ss = 0.0;
    psi = 0.0; // this value will contain the sum of all infeasibilities
    ip = 0;    // ip will be the index of the chosen violated constraint

    if(done_with_active_set)
    {
      s.segment(chunk_start,chunk_size)            = ci0.segment(chunk_start,chunk_size);
      s.segment(chunk_start,chunk_size).noalias() += CI.middleCols(chunk_start,chunk_size).transpose() * x;

      flops0 += chunk_size * n;

      for (int i = chunk_start; i < chunk_start+chunk_size; i++)
      {
  #ifdef QP_SKIP_DEGENERATE
        iaexcl[i] = true;
  #endif
        if(iai[i]!=-1)
          psi += std::min(0.0, s[i]);
      }
    }
    else
    {
      for(int i=0; i<active_set->size(); ++i)
      {
        int ii = active_set->coeff(i);
        s(ii) = ci0(ii) + CI.col(ii).transpose() * x;
        if(iai[ii]!=-1)
          psi += std::min(0.0, s[ii]);
      }
    }

    if (std::abs(psi) <= std::max(chunk_size,iq) * std::numeric_limits<double>::epsilon() * c1 * c2 * 100.0)
    {
      if(!done_with_active_set)
      {
        done_with_active_set = true;
        continue;
      }
      // numerically there are not infeasibilities anymore for the current subset
      if(scheduling.mode==Scheduling::SlidingWindows && chunk_size<m)
      {
        // If we haven't considered all inequalities, then enlarge the current windows, and restart over
        TRACE_SOLVER( std::cout << "increase windows size after " << iter << " iterations to " << chunk_size*2 );

        if(scheduling.grow_factor>1)
        {
          // in this case, let's enlarge the current window size, and re-start sliding from the head,
          // until the window covers the entire set of constriants
          chunk_start = 0;
          chunk_size = std::min(m,int(chunk_size*scheduling.grow_factor));
          continue;
        }
        else
        {
          // otherwise, the window does not grow,
          // so we increase the number of sequential chunks which are satisfied,
          // and if we are not done, let's move to the next chunk.
          ++nb_satisfied_chunks;
          if(nb_satisfied_chunks<nb_chunks)
          {
            Scheduling::shift_window(chunk_start, chunk_size, m);
            continue;
          }
        }
      }
      if(scheduling.mode==Scheduling::WorstSetFirst && !done_with_current_subset)
      {
        // Pack the active set to the front
        int j = 0;
        for(int i=p; i<iq; ++i)
        {
          int Ai = A[i];
          if(Ai>=iq-p)
          {
            while(iai[j]==-1) j++;
            std::swap(s(Ai), s(j));
            std::swap(ci0(Ai), ci0(j));
            CI.col(Ai).swap(CI.col(j));
            A[i] = j;
            j++;
          }
        }
        for (int i = 0; i < m; i++) iai[i] = i;
        for(int i=p; i<iq; ++i)
          iai[i] = -1;

        // Update all remaining inequalities
        s.tail(m-iq)            = ci0.tail(m-iq);
        s.tail(m-iq).noalias() += CI.rightCols(m-iq).transpose() * x;
        flops1 += (m-iq) * n;

        // Pack worst inequalities right after the active set
        chunk_size = std::min(m,std::max(chunk_size,iq*2));
        partial_sort(s,CI,ci0,iq,chunk_size-iq);

        done_with_current_subset = true;
        continue;
      }

      q = iq;
      break;
    }
    ++iter;
    done_with_current_subset = false; // for WorstSetFirst scheduling
    nb_satisfied_chunks = 0; // for sliding window scheduling

#ifdef QP_SKIP_DEGENERATE
    // save old values for u, A, and x
    u_old.head(iq) = u.head(iq);
    A_old.head(iq) = A.head(iq);
    x_old = x;
#endif

//  l2:
    // Step 2: check for feasibility and determine a new S-pair

    // Find largest violated inequality -> (ip,ss)
    if(done_with_active_set)
    {
      for (int i = chunk_start; i < chunk_start+chunk_size; i++)
      {
  #ifdef QP_SKIP_DEGENERATE
        if (s[i] < ss && iai[i] != -1 && iaexcl[i])
  #else
        if (s[i] < ss && iai[i] != -1)
  #endif
        {
          ss = s[i];
          ip = i;
        }
      }
    }
    else
    {
      for(int ii=0; ii<active_set->size(); ++ii)
      {
        int i = active_set->coeff(ii);
  #ifdef QP_SKIP_DEGENERATE
        if (s[i] < ss && iai[i] != -1 && iaexcl[i])
  #else
        if (s[i] < ss && iai[i] != -1)
  #endif
        {
          ss = s[i];
          ip = i;
        }
      }
    }
    if (ss >= 0.0)
    {
      std::cerr << "no violated constraint found psi=" << psi << "  m=" << m << " / " << chunk_size << "\n";
      if(scheduling.mode==Scheduling::WorstFirst)
        break;
      abort();
    }

    np = CI.col(ip);    // set np = n[ip]
    u[iq] = 0.0;        // set u = [u 0]^T
    A[iq] = ip;         // add ip to the active set A

    TRACE_SOLVER( std::cout << "Trying with constraint " << ip << std::endl );
    TRACE_SOLVER( std::cout << "np: " << np.transpose() << std::endl );

    // Project onto constraint ip
    int count = 0;
    while(true)
    {
      count++;
      // Step 2a: determine step direction

      // Compute z = H n+ (step direction in the primal space), and if q>0, r = N* n+ (the negative of the step direction in the dual space)
      compute_step_direction(z, r, np, J, R, iq, d);

      TRACE_SOLVER( std::cout << "Found a new step direction z" << std::endl );

      // Step 2b: compute step length
      int l = 0;
      // Compute t1: partial step length (maximum step in dual space without violating dual feasibility)
      t1 = inf;
      // find the index l s.t. it reaches the minimum of u+[x] / r
      for (int k = p; k < iq; k++)
      {
        double u_over_r;
        if((r[k] > 0.0) && ( (u_over_r = u[k] / r[k]) < t1))
        {
          t1 = u_over_r;
          l = A[k];
        }
      }

      // Compute t2: full step length (minimum step in primal space such that the constraint ip becomes feasible)
      if (z.squaredNorm() > std::numeric_limits<double>::epsilon()) // i.e. z != 0
        t2 = -s[ip] / z.dot(np);
      else
        t2 = inf;

      // the step length is chosen as the minimum of t1 and t2
      t = std::min(t1, t2);
      TRACE_SOLVER( std::cout << "Step sizes: " << t << " (t1 = " << t1 << ", t2 = " << t2 << ") " );

      // Step 2c: determine new S-pair and take step:

      if (t >= inf)       // case (i): no step in primal or dual space
      {
        TRACE_SOLVER_( std::cout << "  no step in primal or dual space (" << iq << ")");
        // QPP is infeasible
        // FIXME: unbounded to raise
        q = iq;
        return inf;
      }
      else if (t2 >= inf) // case (ii): step in dual space using t1
      {
        // set u = u +  t * [-r 1] and drop constraint l from the active set A
        u.head(iq) -= t * r.head(iq);
        u[iq] += t;
        iai[l] = l;
        delete_constraint(R, J, A, u, n, p, iq, l);
        TRACE_SOLVER_( std::cout << " step in dual space only " );
      }
      else                // case (iii): step in both primal and dual space
      {
        x += t * z;                                       // update solution
        f_value += t * z.dot(np) * (0.5 * t + u[iq]);     // update the solution value
        u.head(iq) -= t * r.head(iq);                     // u = u + t * [-r 1]
        u[iq] += t;

        TRACE_SOLVER( std::cout << " step in both spaces, f_value=" << f_value );

        if (std::abs(t - t2) >= std::numeric_limits<double>::epsilon())
        {
          // Partial step has been taken
          if(m<1000){TRACE_SOLVER( std::cout << "Partial step of length " << t << " (" << std::abs(t-t2) << ")" );}

          // drop constraint l
          iai[l] = l;
          delete_constraint(R, J, A, u, n, p, iq, l);

          // update value of ip-th inequality (needed to compute full step length t2)
          s[ip] = ci0[ip] + CI.col(ip).dot(x);
        }
        else
        {
          // Full step has been taken
          TRACE_SOLVER( std::cout << "Full step of length " << t << "  (" << std::abs(t-t2) << ")" );

          // add constraint ip to the active set
      #ifdef QP_SKIP_DEGENERATE
          if (!add_constraint(R, J, d, iq, R_norm))
          {
            TRACE_SOLVER_( std::cout << "Failed to add constraint in active set" );
            iaexcl[ip] = false;
            delete_constraint(R, J, A, u, n, p, iq, ip);

            for (int i = 0; i < m; i++) iai[i] = i;
            for (int i = p; i < iq; i++)
            {
              A[i] = A_old[i];
              u[i] = u_old[i];
              iai[A[i]] = -1;
            }
            x = x_old;
            goto l2; // go to step 2
          }
      #else
          if(!add_constraint(R, J, d, iq, R_norm))
          {
            TRACE_SOLVER_( std::cout << "Numerical issue when adding constraint in active set, abort" );
            q = iq;
            return inf;
          }
      #endif

          iai[ip] = -1;

          if(scheduling.mode==Scheduling::SlidingWindows && done_with_active_set)
            Scheduling::shift_window(chunk_start, chunk_size, m);

          // we are done with current constraint
          break;
        }
      }
    } // end projection loop
  } // end main iteration loop

  if(iter>10*m)
  {
    std::cout << "too many iterations...\n";
    return std::numeric_limits<double>::max();
  }

  g_u = u;

  if(active_set)
    *active_set = A.head(iq);

  TRACE_SOLVER_(std::cout << "#iterations = " << iter << " ; inequality update costs: " << 1e-9*flops0 << " + " << 1e-9*flops1 << " = " << 1e-9*(flops0+flops1) << " GFlops");

  return f_value;
}

void compute_step_direction(Ref<VectorXd> z, Ref<VectorXd> r, Ref<const VectorXd> np, Ref<const MatrixXd> J, Ref<const MatrixXd> R, int iq, Ref<VectorXd> d)
{
  // Compute z = H np and r = N* np using equations 4.7 and 4.8
  int n = z.size();
  d.noalias() = J.transpose() * np;
  z.noalias() = J.middleCols(iq,n-iq) * d.segment(iq,n-iq);
  r.head(iq)  = R.topLeftCorner(iq,iq).triangularView<Upper>().solve(d.head(iq));
}

bool add_constraint(Ref<MatrixXd> R, Ref<MatrixXd> J, Ref<VectorXd> d, int& iq, double& R_norm)
{
  int n = d.size();
  TRACE_SOLVER( std::cout << "Add constraint " << iq << '/' );

  // we have to find the Givens rotation which will reduce the element d[j] to zero.
  // if it is already zero we don't have to do anything, except of decreasing j
  for (int j = n - 1; j >= iq + 1; j--)
  {
    // The Givens rotation is done with the matrix (cc cs, cs -cc).
    // If cc is one, then element (j) of d is zero compared with element
    // (j - 1). Hence we don't have to do anything.
    // If cc is zero, then we just have to switch column (j) and column (j - 1)
    // of J. Since we only switch columns in J, we have to be careful how we
    // update d depending on the sign of gs.
    // Otherwise we have to apply the Givens rotation to these columns.
    // The i - 1 element of d has to be updated to h.
    double cc = d[j - 1];
    double ss = d[j];
    double h = hypot(cc, ss);
    if (h < std::numeric_limits<double>::epsilon()) // h == 0
      continue;
    d[j] = 0.0;
    if (cc < 0.0) h = -h;
    ss = ss / h;
    cc = cc / h;
    d[j-1] = h;

    Ref<VectorXd> J0 = J.col(j-1);
    Ref<VectorXd> J1 = J.col(j);
    internal::apply_rotation_in_the_plane(J0, J1, JacobiRotation<double>(cc,ss));
  }
  // update the number of constraints added
  iq++;

  // Update R:
  R.col(iq-1).head(iq) = d.head(iq);

  TRACE_SOLVER( std::cout << iq << std::endl );
  TRACE_SOLVER( std::cout << "R:\n" << R.topLeftCorner(iq,iq) << std::endl );
  TRACE_SOLVER( std::cout << "J:\n" << J << std::endl  );
  TRACE_SOLVER( std::cout << "d: " << d.head(iq).transpose() << std::endl );

  if (std::abs(d[iq - 1]) <= std::numeric_limits<double>::epsilon() * R_norm)
  {
    // problem degenerate
    return false;
  }
  R_norm = std::max<double>(R_norm, std::abs(d[iq - 1]));
  return true;
}

void delete_constraint(Ref<MatrixXd> R, Ref<MatrixXd> J, VectorXi& A, Ref<VectorXd> u, int n, int p, int& iq, int l)
{
  TRACE_SOLVER( std::cout << "Delete constraint " << l << ' ' << iq );
  int qq = -1; // just to prevent warnings from smart compilers

  // Find the index qq for active constraint l to be removed
  // FIXME should be doable in O(1)
  for (int i = p; i < iq; i++)
    if (A[i] == l)
    {
      qq = i;
      break;
    }
  assert(qq>=0);

  // remove the constraint from the active set and the duals
  for (int i = qq; i < iq - 1; i++)
  {
    A[i] = A[i + 1];
    u[i] = u[i + 1];
    R.col(i) = R.col(i+1); // FIXME R is upper triangular so we can probably reduce this copy
  }

  A[iq - 1] = A[iq];
  u[iq - 1] = u[iq];
  A[iq] = 0;
  u[iq] = 0.0;
  R.col(iq-1).head(iq).setZero();
  // constraint has been fully removed
  iq--;
  TRACE_SOLVER( std::cout << '/' << iq << std::endl );

  if (iq == 0)
    return;

  for (int j = qq; j < iq; j++)
  {
    double cc = R(j,j);
    double ss = R(j + 1,j);
    double h = hypot(cc, ss);
    if (h < std::numeric_limits<double>::epsilon()) // h == 0
      continue;

    if (cc < 0.0) h = -h;
    cc = cc / h;
    ss = ss / h;
    R(j + 1,j) = 0.0;
    R(j,j) = h;

    Ref<VectorXd,0,InnerStride<> > R0 = R.row(j).segment(j+1,iq-j-1);
    Ref<VectorXd,0,InnerStride<> > R1 = R.row(j+1).segment(j+1,iq-j-1);;
    internal::apply_rotation_in_the_plane(R0, R1, JacobiRotation<double>(cc,ss));

    Ref<VectorXd> J0 = J.col(j);
    Ref<VectorXd> J1 = J.col(j+1);
    internal::apply_rotation_in_the_plane(J0, J1, JacobiRotation<double>(cc,ss));
  }
}

inline double hypot(double a, double b)
{
  double a1, b1, t;
  a1 = std::abs(a);
  b1 = std::abs(b);
  if (a1 > b1)
  {
    t = (b1 / a1);
    return a1 * std::sqrt(1.0 + t * t);
  }
  else if (b1 > a1)
  {
    t = (a1 / b1);
    return b1 * std::sqrt(1.0 + t * t);
  }
  return a1 * std::sqrt(2.0);
}

}
