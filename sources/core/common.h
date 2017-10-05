/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2014 CNRS
   Copyright (C) 2013, 2014, 2015, 2016, 2017 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#pragma once

#include <vector>
#include <iostream>
#include <cassert>
#include <cstring>

// Math Constants exist in Windows but they need to be included
// by defining the constant _USE_MATH_DEFINES
// Cf https://msdn.microsoft.com/en-us/library/4hwaceh6%28v=vs.110%29.aspx
#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#else
#include <cmath>
#endif

#include <cstring>
#include <algorithm>

#include <string>
#include <map>

#include <Eigen/Core>

typedef Eigen::VectorXd vec;
typedef Eigen::Ref<vec> vecref;
typedef Eigen::Ref<const vec> const_vecref;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    RowMatrixXd;

// Convenience functions.
static inline double norm(const vec& v)
{
    return v.norm();
}

static inline vec normalize(const vec& v)
{
    vec result = v;
    result.normalize();
    return result;
}

static inline double dot(const vec& a, const vec &b)
{
    return a.dot(b);
}

//! \brief If A and B have equal sizes, return a vector of the same size
//! whose elements result from pairwise multiplications of the elements of A
//! and B.  Otherwise, if A or B is a one-element vector, return the scalar
//! product of the other vector with that element.  It is an error to supply
//! vectors of different sizes.
extern vec product(const vec& a, const vec& b);

std::ostream& operator<<(std::ostream& out, const vec& v);
template<typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& v) {
    out << "[";
    for(int i=0; i<v.size(); ++i)
    {
        out << v[i];
        if(i != v.size()-1) { out << ", "; }
    }
    out << "]";
    return out;
}

namespace alta {

//! \brief locate the first index of value v in vector vec. Complexity in
//! O(n) is the worst case.
template<typename T> int is_in(std::vector<T> ve, T v)
{
    int res = -1;
    for(unsigned int i=0; i<ve.size(); ++i)
    {
        if(ve[i] == v)
            return i;
    }

    return res;
}

template<typename T> T clamp(T x, T a, T b)
{
	return std::max<T>(std::min<T>(x, b), a);
}

template<typename T>
static bool close_to(const T a, const T b, const T epsilon = 1E-7)
{
    return std::abs(a - b) < epsilon;
}


/* Mark a type, class, method, function, or variable as deprecated.  */
#ifdef __GNUC__
# define ALTA_DEPRECATED  __attribute__((__deprecated__))
#else
# define ALTA_DEPRECATED
#endif

#ifdef _MSC_VER
#define NOT_IMPLEMENTED() \
std::cerr << "<<ERROR>> not implemented " << __FUNCDNAME__ << " in file " << __FILE__ \
          << ":" << __LINE__ << std::endl; \
throw
#else
#define NOT_IMPLEMENTED() \
std::cerr << "<<ERROR>> not implemented " << __PRETTY_FUNCTION__ << " in file " << __FILE__ \
          << ":" << __LINE__ << std::endl; \
throw
#endif

// Mathematical definition not provided on the Window plateform
#ifdef _MSC_VER
template<typename T> T erf(T x)
{
    // constants
    const double a1 =  0.254829592;
    const double a2 = -0.284496736;
    const double a3 =  1.421413741;
    const double a4 = -1.453152027;
    const double a5 =  1.061405429;
    const double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0) {
        sign = -1;
	 }
    x = fabs(x);

    // A&S formula 7.1.26
    const double t = 1.0/(1.0 + p*x);
    const double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return sign*y;
}

template<typename T> T erfc( T x )
{
  return 1 - erf(x);
}
#endif // END OF erf and erc not provided by Windows plateform 


#if __cplusplus >= 201103L || (__cplusplus < 200000 && __cplusplus > 199711L)
template<typename T> bool isnan(T x)
{
    return std::isnan(x);
}
#else
template<typename T> bool isnan(T x)
{
    return (x != x );
}
#endif


#ifdef _WIN32
#define ALTA_DLL_EXPORT extern "C" __declspec(dllexport)
#else
#define ALTA_DLL_EXPORT extern "C"
#endif
}

namespace alta {
   class timer;
}
std::ostream& operator<<(std::ostream& out, const alta::timer& t);
namespace alta {

/*! \brief The timer class: a cross plateform timing solution
 */
class timer
{
    public:

        //! \brief Constructor
        timer();

        //! \brief Start the timer
        void start();

        //! \brief Stop the timer
        void stop();

        //! \brief Reset the timer
        void reset();

        //! \brief Compute the time taken since the start, accounting for
        //! all pauses.
        int elapsed() const;

        //! \brief Print the time to the standard output
        void print(std::ostream& out) const;

        //! \brief ostream compliant operator
        friend std::ostream& ::operator<<(std::ostream& out, const timer& t);

    private:

        //! \brief get the current time in seconds. This is a system dependant
        //! function.
        //! \internal
        unsigned int current_time() const;

        unsigned int _start, _stop;
        unsigned int _elapsed;
};


// I/O error handling in user interfaces.

#include <iostream>
#include <cstdlib>

#ifdef _WIN32

// We cannot rely on 'errno' on Windows.

# define ALTA_FILE_IO_ERROR_STRING(e)						\
    "(unspecified I/O error)"

#else

# define ALTA_FILE_IO_ERROR_STRING(e)						\
    (strerror(errno))

#endif

# define CATCH_FILE_IO_ERROR(file)												\
		catch (std::ios_base::failure& e) {										\
				std::cerr << "<<ERROR>> failed to load '"					\
									<< (file) << "'"												\
									<< ": " << ALTA_FILE_IO_ERROR_STRING(e)	\
									<< std::endl;														\
				exit(EXIT_FAILURE);																\
		}
}
