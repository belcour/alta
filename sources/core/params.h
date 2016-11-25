/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2014 CNRS
   Copyright (C) 2013, 2014, 2015, 2016 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#pragma once

#include <string>
#include <map>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cassert>
#include <iostream>

#include <Eigen/Geometry>

#include "common.h"

namespace alta {

/*! \class params
 *  \ingroup core
 *  \brief a static class allowing to change from one parametrization
 *  to another.
 *
 *  Any function object or data object should have an associated
 *  parametrization.
 *
 *  We use the following convention to defined the tangent, normal and
 *  bi-normal of the surface:
 *   * The normal is the upper vector (0, 0, 1)
 *   * The tangent direction is along x direction (1, 0, 0)
 *   * The bi-normal is along the y direction (0, 1, 0)
 */
class params
{
    public: // data

     //! \brief list of all supported parametrization for the input space.
     //! An unsupported parametrization will go under the name ! *unknown*. We
     //!use the following notations:
     //!   * The View vector is \f$V\f$
     //!   * The Light vector is \f$L\f$
     //!   * The Normal vector is \f$N\f$
     //!   * The Reflected vector is \f$R = 2\mbox{dot}({V}, {N}){N} - {V}\f$
     //!   * The Half vector is \f$H = \frac{V+L}{||V+L||}\f$
     //!   * The Back vector is \f$K = \frac{V-L}{||V-L||}\f$
     //!   * The elevation angle of vector \f$V\f$ is \f$\theta_V\f$
     //!   * The azimuth angle of vector \f$V\f$ is \f$\phi_V\f$
       //!
       //! [stark]: http://dl.acm.org/citation.cfm?id=1042387
       //! [neumann]: https://www.cg.tuwien.ac.at/research/publications/1996/Neumann-1996-NCB/
       //!
     enum input
     {
       RUSIN_TH_PH_TD_PD,     /*!< Half-angle parametrization as described by Rusinkiewicz 
                                  [[1998]][rusin] described in spherical coordinates. Coordinates
                                  are: \f$ \left[ \theta_H, \phi_H, \theta_D, \phi_D \right] \f$.
                                  [rusin]: http://link.springer.com/chapter/10.1007/978-3-7091-6453-2_2
                               */
       RUSIN_TH_PH_TD,
     RUSIN_TH_TD_PD,      /*!< Isotropic Half-angle parametrization. This
                      parametrization is defined by the tuple \f$
                      \vec{x} = [\theta_H, \theta_D, \phi_D] \f$*/
       RUSIN_TH_TD,           /*!< Half-angle parametrization with no azimutal information */
       RUSIN_VH_VD,           /*!< Half-angle parametrization in vector format. Coordinates are:
                                   \f$ [H_x, H_y, H_z, D_x, D_y, D_z] \f$.*/
       RUSIN_VH,              /*!< Half-angle parametrization with no difference direction in 
                             vector format. Coordinates are: [\f$\vec{h}_x, \vec{h}_y, 
                            \vec{h}_z\f$]. */
       COS_TH_TD,             /*!< Cosine of the RUSIN_TH_TD parametrization: Coordinates are in 
                                   \f$ [\cos_{\theta_H},\cos_{\theta_D}] \f$. */
       COS_TH,

       SCHLICK_TK_PK,         /*!< Schlick's back vector parametrization */
       SCHLICK_VK,            /*!< Schlick's back vector */
       SCHLICK_TL_TK_PROJ_DPHI,/*!< 3D Parametrization where the phi component is projected and
                              the parametrization is centered around the back direction.
                     \f$[\theta_L, x, y] = [\theta_L, \theta_K \cos(\phi_K), \theta_K \sin(\phi_K)]\f$*/
       COS_TK,                /*!< Schlick's back vector dot product with the normal */


       RETRO_TL_TVL_PROJ_DPHI,/*!< 2D Parametrization where the phi component is projected and
                               the parametrization is centered around the retro direction
                      \f$[x, y] = [\theta_{VL} \cos(\Delta\phi), \theta_{VL} 
                      \sin(\Delta\phi)]\f$.*/

       STEREOGRAPHIC,         /*!< Stereographic projection of the Light and View vectors */


       SPHERICAL_TL_PL_TV_PV, /*!< Light and View vectors represented in spherical coordinates */
       COS_TLV,               /*!< Dot product between the Light and View vector */
       COS_TLR,               /*!< Dot product between the Light and Reflected vector */
       ISOTROPIC_TV_TL,       /*!< Light and View vectors represented in spherical coordinates, */
       ISOTROPIC_TV_TL_DPHI,  /*!< Light and View vectors represented in spherical coordinates,
                                   with the difference of azimutal coordinates in the last component  */
       ISOTROPIC_TV_PROJ_DPHI,/*!< 2D Parametrization where the phi component is projected.
                               Coordinates are: \f$\left[\theta_v \cos(\Delta\phi), \theta_v 
                      \sin(\Delta\phi)\right]\f$.*/
       ISOTROPIC_TL_TV_PROJ_DPHI,/*!< 3D Parametrization where the phi component is projected.
                                  Coordinates are: \f$\left[\theta_l, \theta_v \cos(\Delta\phi), 
                        \theta_v \sin(\Delta\phi)\right]\f$.*/
       ISOTROPIC_TD_PD,       /*!< Difference between two directions such as R and H */

        STARK_2D,             /*!< Modified Stark et al. [[2005]][stark] 2D parametrization by 
                                   Barla et al. [[2015]][barla]. This parametrization is defined 
                                   by the couple \f$ \vec{x} = ||\tilde{H}_p||, ||\tilde{B}|| \f$, 
                                   where \f$ \tilde{H} = \frac{1}{2}(L+V) \f$
                        and \f$ \tilde{B} = \frac{1}{2}(L-V) \f$. \f$ \tilde{H}_p \f$
                        is the projected coordinates of \f$ \tilde{H} \f$ on the
                        tangent plane.
                                      [stark]: http://dl.acm.org/citation.cfm?id=1042387 
                                      [barla]: http://todo/
                               */

        STARK_3D,             /*!< Modified Stark et al. [[2005]][stark] 2D parametrization by 
                                   Barla et al. [[2015]][barla]. This parametrization is defined 
                                   by the tuple \f$ \vec{x} = \left[||\tilde{H}_p||, ||\tilde{B}||,
                                   \phi_B-\phi_H \right] \f$. 
                                      [stark]: http://dl.acm.org/citation.cfm?id=1042387 
                                      [barla]: http://todo/
                               */
        NEUMANN_2D,           /*!< Neumann and Neumann [[1996]][neumann] parametrization. This 
                                   parametrization is defined by the couple \f$ \vec{x} = \left[
                                   ||\tilde{H}_p||, ||\tilde{B}_p||\right] \f$, where \f$ 
                                   \tilde{H} = \frac{1}{2}(L+V) \f$
                        and \f$ \tilde{B} = \frac{1}{2}(L-V) \f$. \f$ \tilde{H}_p \f$
                        is the projected coordinates of \f$ \tilde{H} \f$ on the
                        tangent plane.
                                      [neumann]: https://www.cg.tuwien.ac.at/research/publications/1996/Neumann-1996-NCB/
                               */
        NEUMANN_3D,           /*!< Neumann and Neumann [[1996]][neumann] 3D parametrization. This 
                                   parametrization is defined by the tuple \f$ \vec{x} = 
                                   [||\tilde{H}_p||, ||\tilde{B}_p||, \phi_B-\phi_H] \f$.
                                      [neumann]: https://www.cg.tuwien.ac.at/research/publications/1996/Neumann-1996-NCB/
                               */
       CARTESIAN,             /*!< View and Light vectors represented in cartesian coordinates.
                               We always pack the view vector first: \f$[V.x, V.y, 
                      V.z, L.x, L.y, L.z] \f$*/

       UNKNOWN_INPUT          /*!< Default behaviour. Only use this is you do not fit BRDF data. */
     };

     //! \brief list of all supported parametrization for the output space.
     //! An unsupported parametrization will go under the name
     //! <em>unknown</em>.
     enum output
     {
       INV_STERADIAN,                /*!< Output values in inverse steradian (sr-1). 
                               This is the standard definition for a BRDF. */
       INV_STERADIAN_COSINE_FACTOR,  /*!< Output values in inverse steradian (sr-1)
                                          weighted by the cosine factor of the output
                              direction. */
       ENERGY,
       RGB_COLOR,
       XYZ_COLOR,
       UNKNOWN_OUTPUT
     };

    public: // methods

        //! \brief parse a string to provide a parametrization type.
        static params::input parse_input(const std::string& txt);

        //! \brief parse a string to provide a parametrization type.
        static params::output parse_output(const std::string& txt);

      //! \brief look for the string associated with a parametrization
      //! type.
      static const std::string& get_name(const params::input param);

      //! \brief look for the string associated with a parametrization
      //! type.
      //! \todo Finish this implementation. It requires another static
      //! object.
      static const std::string& get_name(const params::output);

        //! \brief static function for input type convertion. This
        //! function allocate the resulting vector.
        static double* convert(const double* invec, params::input intype,
                               params::input outtype)
        {
            int dim = dimension(outtype); // Get the size of the output vector

            if(dim > 0)
            {
                double* outvec = new double[dim];
                double  temvec[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // Temp CARTESIAN vectors
                to_cartesian(invec, intype, temvec);
                from_cartesian(temvec, outtype, outvec);

                return outvec;
            }
            else
            {
                return NULL;
            }
        }

        //! \brief static function for input type convertion. The outvec
        //! resulting vector should be allocated with the correct
        //! output size.
        //! TODO:  RP this function is weird the outtype and intype ARE NEVER USED
        static void convert(const double* invec, params::output intype, int indim,
                            params::output outtype, int outdim, double* outvec)
        {
          // The convertion is done using the cartesian parametrization as
          // an intermediate one. If the two parametrizations are equals
          // there is no need to perform the conversion.
          if(outdim == indim)
          {
            for(int i=0; i<outdim; ++i) { outvec[i] = invec[i]; }
          }
          // If the dimension of the output data is bigger than the 
          // dimensions of the input domain, and the input domain is of
          // dimension one, spread the data over all dimensions.
          else if(indim == 1)
          {
            for(int i=0; i<outdim; ++i) { outvec[i] = invec[0]; }
          }
          // If the output dimension is one, compute the average of the
          // input vector values.
          else if(outdim == 1)
          {
            for(int i=0; i<indim; ++i) 
            { 
                outvec[0] += invec[i]; 
            } 
            outvec[0] /= static_cast<double>(indim);
          }
          else
          {
            NOT_IMPLEMENTED();
          }
        }

        //! \brief static function for input type convertion. The outvec
        //! resulting vector should be allocated with the correct
        //! output size.
        static void convert(const double* invec, params::input intype,
                            params::input outtype, double* outvec)
        {
        // The convertion is done using the cartesian parametrization as
        // an intermediate one. If the two parametrizations are equals
        // there is no need to perform the conversion.
        if(intype == outtype)
        {
          int dim = dimension(outtype);
          for(int i=0; i<dim; ++i) { outvec[i] = invec[i]; }
        }
        // If the input parametrization is the CARTESIAN param, then 
        // there is no need to transform the input data.
        else if(intype == params::CARTESIAN)
        {
          from_cartesian(invec, outtype, outvec);
        }
        // If the output parametrization is the CARTESIAN param, then
        // there is no need to convert back to another param.
        else if(outtype == params::CARTESIAN)
        {
          to_cartesian(invec, intype, outvec);
        }
        else
        {
          // temporary CARTESIAN vector
          double  temvec[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

          to_cartesian(invec, intype, temvec);
          from_cartesian(temvec, outtype, outvec);
        }
        }

        //! \brief convert a input vector in a given parametrization to an
        //! output vector in a cartesian parametrization, that is two 3d
        //! vectors concatenated.
        static void to_cartesian(const double* invec, params::input intype,
                                 double* outvec);

        //! \brief convert a input CARTESIAN vector, that is two 3d vectors
        //! concatenated  to an output vector in a given parametrization.
        static void from_cartesian(const double* invec, params::input outtype,
                                   double* outvec);

        //! \brief provide a dimension associated with a parametrization
        static int  dimension(params::input t);

        //! \brief provide a dimension associated with a parametrization
        static int  dimension(params::output t)
        {
            switch(t)
            {
                // 1D Parametrizations
                case params::INV_STERADIAN:
                case params::ENERGY:
                    return 1;
                    break;

                // 3D Parametrization
                case params::RGB_COLOR:
                case params::XYZ_COLOR:
                    return 3;
                    break;

                default:
                    assert(false);
                    return -1;
                    break;
            }
        }

      //! \brief Is the value stored weighted by a cosine factor
      static bool is_cosine_weighted(params::output t)
      {
        switch(t)
        {
          case params::INV_STERADIAN_COSINE_FACTOR:
            return true;
            break;

          case params::INV_STERADIAN:
          case params::ENERGY:
          case params::RGB_COLOR:
          case params::XYZ_COLOR:
          default:
            return false;
            break;
        }
      }

      //! \brief from the 4D definition of a half vector parametrization,
      //! export the cartesian coordinates.
      template<typename Real>
      static void half_to_cartesian(Real theta_h, Real phi_h,
                                    Real theta_d, Real phi_d, Real* out)
      {
          // Calculate the half vector
          Real half[3];
          half[0] = sin(theta_h) * cos(phi_h);
          half[1] = sin(theta_h) * sin(phi_h);
          half[2] = cos(theta_h);

          // Compute the light vector using the rotation formula.
          // TODO: Eventually make OUT an Eigen vector.
          out[0] = sin(theta_d)*cos(phi_d);
          out[1] = sin(theta_d)*sin(phi_d);
          out[2] = cos(theta_d);

          // Rotate the diff vector to get the output vector
          Eigen::Matrix<Real, 3, 1> out2(out[0], out[1], out[2]);
          rotate_binormal(out2, theta_h);
          rotate_normal(out2, phi_h);
          out[0] = out2[0];
          out[1] = out2[1];
          out[2] = out2[2];

          // Compute the out vector from the in vector and the half
          // vector.
          const Real dot = out[0]*half[0] + out[1]*half[1] + out[2]*half[2];
          out[3] = -out[0] + 2.0*dot * half[0];
          out[4] = -out[1] + 2.0*dot * half[1];
          out[5] = -out[2] + 2.0*dot * half[2];

#ifdef DEBUG
          assert(out[2] >= 0.0 && out[5] >= 0.0);
#endif
      }

      //! \brief from the 4D definition of a classical vector parametrization,
      //! export the cartesian coordinates.
      static void classical_to_cartesian(double theta_l, double phi_l, 
                                         double theta_v, double phi_v, double* out)
      {
        out[0] = cos(phi_l)*sin(theta_l);
        out[1] = sin(phi_l)*sin(theta_l);
        out[2] = cos(theta_l);
        out[3] = cos(phi_v)*sin(theta_v);
        out[4] = sin(phi_v)*sin(theta_v);
        out[5] = cos(theta_v);
      }

      //! \brief rotate a cartesian vector with respect to the normal of
      //! theta degrees.
      template<class Vector, typename Real>
      static void rotate_normal(Vector& vec, Real theta)
      {
          auto rotz = Eigen::AngleAxis<Real>(theta, Eigen::Matrix<Real, 3, 1>(0, 0, 1));
          vec = rotz * vec;
      }

      //! \brief rotate a cartesian vector with respect to the bi-normal of
      //! theta degrees.
      template<class Vector, typename Real>
      static void rotate_binormal(Vector& vec, Real theta)
      {
          auto roty = Eigen::AngleAxis<Real>(theta, Eigen::Matrix<Real, 3, 1>(0, 1, 0));
          vec = roty * vec;
      }

      //! \brief Rotate VEC by THETA around the normal, and by PHI around the
      //! binormal.
      template<class Vector, typename Real>
      static void rotate(Vector& vec, Real theta, Real phi)
      {
          // Naive implementation.
          rotate_normal(vec, theta);
          rotate_binormal(vec, phi);
      }

      static void print_input_params();

      //! \brief Check whether or not a (light,view) configuration is above the hemisphere
      //! \return true if both vectors (light and view) are above and false otherwise
      //! \remark the function returns true also for grazing, tangential configuration 
      //! (i.e., where light dot normal == 0 or view dot normal == 0)
      static bool is_above_hemisphere( double* invec,         params::input intype );

      static bool inline is_below_hemisphere( double* invec,         params::input intype )
      {
        return ! is_above_hemisphere(invec, intype);
      }

};

/*! \brief A parameters object. Allow to define function object (either data
 *  or functions that are defined over an input space and output space. This
 *  Object allowas to change the parametrization of the input or output space.
 */
class parameters
{
  public:
    parameters()
      : _in_param(params::UNKNOWN_INPUT),
        _out_param(params::UNKNOWN_OUTPUT),
        _nX( 0 ),
        _nY( 0 )
    {
    }

    parameters(unsigned int dim_X, unsigned int dim_Y,
               params::input in_param, params::output out_param)
    : _in_param(in_param),
      _out_param(out_param),
      _nX( dim_X ),
      _nY( dim_Y )
    {}

    parameters( unsigned int dim_X, unsigned int dim_Y)
        ALTA_DEPRECATED
    : _in_param(params::UNKNOWN_INPUT),
      _out_param(params::UNKNOWN_OUTPUT),
      _nX( dim_X ),
      _nY( dim_Y )
    {}

    parameters(params::input in_param, params::output out_param)
        ALTA_DEPRECATED
    : _in_param( in_param ),
      _out_param( out_param ),
      _nX( params::dimension(_in_param) ),
      _nY( params::dimension(_out_param) )
    {
    }

    //! \brief provide the input parametrization of the object.
    params::input input_parametrization() const
    {
      return _in_param;
    }

    //! \brief provide the outout parametrization of the object.
    params::output output_parametrization() const
    {
      return _out_param;
    }

    /* DIMENSION OF THE INPUT AND OUTPUT DOMAIN */

    //! Provide the dimension of the input space of the function
    int dimX() const { return _nX ; }
    //! Provide the dimension of the output space of the function
    int dimY() const { return _nY ; }

    const parameters set_input(int nX, params::input p) const
    {
        return alta::parameters(nX, _nY, p, _out_param);
    }

    const parameters set_output(int nY, params::output p) const
    {
        return alta::parameters(_nX, nY, _in_param, p);
    }

  protected:
    // Input and output parametrization
    params::input  _in_param ;
    params::output _out_param ;

    // Dimension of the function & domain of definition.
    int _nX, _nY ;
};
}
