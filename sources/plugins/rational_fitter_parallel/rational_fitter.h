#pragma once

// Include STL
#include <vector>
#include <string>

// Interface
#include <core/function.h>
#include <core/rational_function.h>
#include <core/data.h>
#include <core/vertical_segment.h>
#include <core/fitter.h>
#include <core/args.h>

/*! \brief A vertical segment fitter for rational functions that can run in parallel
 *  using QuadProg++ quadratic solver.
 *  \ingroup plugins
 *
 *  \details
 *  You can find QuadProg++ here: http://quadprog.sourceforge.net/
 *  <br />
 *
 *  <h3>Plugin parameters</h3>
 *
 *  We provide the following command line arguments to manipulate this plugin:
 *  <ul>
 *		<li><b>--np</b> <em>[int]</em> controls the maximum number of total coefficients
 *		an interpolation should have. By default, this number is 10.</li>
 *		<li><b>--min-np</b> <em>[int]<em></li> controls the starting value for the 
 *		number of coefficients for the rational function. by default, this number
 *		is 10.</li>
 *		<li><b>--np-step</b> <em>[int]</em> stepping for the number of coefficients
 *		of the rational function. By default, this number is 1.</li>
 *		<li><b>--nb-cores</b> <em>[int]</em> number of core allocated to perform
 *		the seach. By default, this is equal to the number of processors.</li>
 *  </ul>
 */
class rational_fitter_parallel : public fitter
{
	public: // methods

		rational_fitter_parallel() ;
		virtual ~rational_fitter_parallel() ;

		// Fitting a data object
		//
		virtual bool fit_data(const data* d, function* fit, const arguments& args) ;

		// Provide user parameters to the fitter
		//
		virtual void set_parameters(const arguments& args) ;

	protected: // methods

		// Fitting a data object using np elements in the numerator and nq 
		// elements in the denominator
		virtual bool fit_data(const vertical_segment* d, int np, int nq, 
		                      rational_function* fit, const arguments &args, 
									 vec& p, vec& q, double& delta, 
									 double& linf_dist,double& l2_dist) ;
		virtual bool fit_data(const vertical_segment* dat, int np, int nq, 
		                      int ny, rational_function_1d* fit, vec& p, vec& q, 
									 double& delta) ;

		//! \brief Create a constraint vector given its index i in the data
		//! object and the rational function object to fit. This function
		//! returns two rows of the constraints matrix A, cu and cl, 
		//! corresponding to the lower constraint and the upper constraint
		//! of the vertical segment.
		virtual void get_constraint(int i, int np, int nq, int ny, 
		                            const vertical_segment* data, 
											 const rational_function_1d* func, 
											 vec& cu, vec& cl);

	protected: // data

		//! Number of points used when starting the adaptive interpolation
		//! procedure. By default, this value is 100.
		int nb_starting_points;
} ;

