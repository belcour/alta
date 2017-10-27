/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2016 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include <core/data.h>
#include <core/data_storage.h>
#include <core/vertical_segment.h>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace alta;

static vec get_min(const params::input& param) {

   vec res(params::dimension(param));
   switch(param)
   {
      // 1D Parametrizations
      case params::COS_TH:
         res[0] = 0.0;
         break;
      case params::COS_TK:
         res[0] = 0.0;
         break;
      case params::COS_TLV:
         res[0] = -1.0;
         break;
      case params::COS_TLR:
         res[0] = -1.0;
         break;

         // 2D Parametrizations
      case params::COS_TH_TD:
         res[0] = 0.0;
         res[1] = 0.0;
         break;
      case params::RUSIN_TH_TD:
         res[0] = 0.0;
         res[1] = 0.0;
         break;
      case params::ISOTROPIC_TV_PROJ_DPHI:
         res[0] = -0.5*M_PI;
         res[1] = -0.5*M_PI;
         break;
      case params::STARK_2D:
         res[0] = 0.0;
         res[1] = 0.0;
         break;

      // 3D Params
      case params::RUSIN_TH_TD_PD:
         res[0] = 0.0;
         res[1] = 0.0;
         res[2] = 0.0;
         break;
      case params::STARK_3D:
         res[0] = 0.0;
         res[1] = 0.0;
         res[2] = 0.0;
         break;

      default:
         return res;
         break;
   }
   return res;
}

static vec get_max(const params::input& param) {

   vec res(params::dimension(param));
   switch(param)
   {
      // 1D Parametrizations
      case params::COS_TH:
         res[0] = 1.0;
         break;
      case params::COS_TK:
         res[0] = 1.0;
         break;
      case params::COS_TLV:
         res[0] = 1.0;
         break;
      case params::COS_TLR:
         res[0] = 1.0;
         break;

         // 2D Parametrizations
      case params::COS_TH_TD:
         res[0] = 1.0;
         res[1] = 1.0;
         break;
      case params::RUSIN_TH_TD:
         res[0] = 0.5*M_PI;
         res[1] = 0.5*M_PI;
         break;
      case params::ISOTROPIC_TV_PROJ_DPHI:
         res[0] = 0.5*M_PI;
         res[1] = 0.5*M_PI;
         break;
      case params::STARK_2D:
         res[0] = 1.0;
         res[1] = 1.0;
         break;

      // 3D Params
      case params::RUSIN_TH_TD_PD:
         res[0] = 0.5*M_PI;
         res[1] = 0.5*M_PI;
         res[2] = 2.0*M_PI;
         break;
      case params::STARK_3D:
         res[0] = 1.0;
         res[1] = 1.0;
         res[2] = 2.0*M_PI;
         break;

      default:
         return res;
         break;
   }
   return res;
}


/*! \ingroup datas
 *  \class data_grid
 *  \brief Data object storing BRDF values on a grid. Can perform linear
 *  interpolation between grid cells

 *  \details
 *  It is possible to select the parametrization using the --param NAME
 *  argument when loading the BRDF. The default parametrization is
 *  \ref params::STARK_2D "STARK_2D"
 *
 *  \author Laurent Belcour <laurent.belcour@umontreal.ca>
 *
 */
class BrdfGrid : public vertical_segment {
   public:

      std::vector<int> _size;

      BrdfGrid(size_t size, const parameters& params, const arguments& args)
        : vertical_segment(params, size)
      {
         initialize(params, args);
      }

      // FIXME: Move this to constructor.
      void initialize(const parameters& params, const arguments& args, bool preallocate = true) {

         // Allow to load a different parametrization depending on the
         // parameters provided.
         const params::input  in_param  = params.input_parametrization();
         const params::output out_param = params.output_parametrization();
         const int nX = params.dimX(), nY = params.dimY();

         // Get the grid size, by default use an 10^dimX grid
         int N = 1;
         if(args.is_defined("grid-size") && args.is_vec("grid-size")) {
            _size = args.get_vec<int>("grid-size");
            assert(_size.size() == nX);

            for(int i=0; i<nX; ++i) {
               N *= _size[i];
            }
         } else {
            for(int i=0; i<nX; ++i) {
               _size.push_back(10);
               N *= 10;
            }
         }
         _min = get_min(in_param);
         _max = get_max(in_param);

         // Preallocate the grid data
         if(preallocate) {
            // Set the abscissa data
            vec x(parametrization().dimX() + parametrization().dimY());
            for(int n=0; n<N; ++n) {
               get_abscissa(n, x);
               set(n, x);
            }
         }

         parameters p(nX, nY, in_param, out_param);
         setParametrization(p);
      }

      inline int get_index(const vec& x) const {

         const auto k = get_indices(x);
         return get_index(k);
      }

      inline int get_index(const std::vector<int>& k) const {

         int N = 0;
         for(int i=0; i<parametrization().dimX(); ++i) {
            N = N*_size[i] + k[i];
         }
         return N;
      }

      inline std::vector<int> get_indices(const vec& x) const {

         std::vector<int> k(parametrization().dimX());
         for(int i = 0; i < parametrization().dimX(); ++i) {
            k[i] = floor((_size[i]-1)*(x[i]-_min[i]) / (_max[i]-_min[i]));
         }
         return k;
      }

      inline std::vector<int> get_indices(int N) const {

         // Variables
         std::vector<int> k(parametrization().dimX());
         int K = N;

         // Generate the k-th coordinate
         k[parametrization().dimX()-1] = K % _size[parametrization().dimX()-1];
         for(int i=parametrization().dimX()-2; i>=0; --i) {
            K    = (K-k[i+1]) / _size[i+1];
            k[i] = K % _size[i];
         }

         return k;
      }

      // Obtain the X dimension of vector x from its index in the grid.
      inline void get_abscissa(int N, vec& x) const {

         const auto k = get_indices(N);

         // Evaluate the abscissa from the vector indices and
         // the min and max. The grid exactly map [0..1]
         for(int i = 0; i < parametrization().dimX(); ++i) {
            x[i] = _min[i] + (_max[i]-_min[i]) * k[i] / (_size[i]-1);
         }
      }

      void save(const std::string& filename) const
      {
         std::ofstream file;

         file.exceptions(std::ios_base::failbit);
         file.open(filename.c_str(), std::ios_base::trunc);
         file.exceptions(std::ios_base::goodbit);

         save_header(file);

         file << std::setprecision(std::numeric_limits<double>::digits10);
         for(int i=0; i < size(); ++i)
         {
            vec x = get(i);
            for(int j=0; j < parametrization().dimX() + parametrization().dimY(); ++j)
            {
               file << x[j] << "\t";
            }
            file << std::endl;
         }

         file.close();
      }

      void save_header(std::ostream& out) const {
         out << "#ALTA HEADER BEGIN" << std::endl;
         out << "#DIM " << parametrization().dimX() << " "
             << parametrization().dimY() << std::endl;
         out << "#PARAM_IN  "
             << params::get_name(parametrization().input_parametrization())
             << std::endl;
         out << "#PARAM_OUT "
             << params::get_name(parametrization().output_parametrization())
             << std::endl;
         out << "#grid-size " << _size << std::endl;
         out << "#ALTA HEADER END" << std::endl;
      }

      inline vec operator[](int i) const
      {
         return get(i) ;
      }

      void set(int id, const vec& x)
      {
         vertical_segment::set(id, x);
      }

      virtual vec value(const vec& x) const
      {
         // Get the base index
         const auto k = get_indices(x);
         const int id = get_index(k);

         // Get the base configuration
         const vec x0 = get(id);

         // Compute the weights
         vec alphas(parametrization().dimX());
         for(int i = 0; i < parametrization().dimX(); ++i) {
            const auto xdiff = (_size[i]-1) * (x[i]-x0[i]) / (_max[i]-_min[i]);
            alphas[i] = clamp(xdiff, 0.0, 1.0);
            assert(alphas[i] >= 0.0 && alphas[i] <= 1.0);
         }

         // Result vector
         vec y = vec::Zero(parametrization().dimY());
         const unsigned int D = pow(2, parametrization().dimX());
         for(unsigned int d=0; d<D; ++d) {

            double alpha = 1.0; // Global alpha
            int    cid_s = 0;   // Id shift

            for(int i=0; i<parametrization().dimX(); ++i) {
               bool bitset = ((1 << i) & d);
               auto calpha =  (bitset) ? alphas[i] : 1.0-alphas[i];

               alpha *= calpha;
               cid_s = cid_s*_size[i] + ((bitset) ? 1 : 0);
            }

            if(id+cid_s > size())
               cid_s = 0;

            y += alpha * get(id+cid_s).segment(parametrization().dimX(),
                                               parametrization().dimY());
         }

         return y;
      }
};

ALTA_DLL_EXPORT data* provide_data(size_t size,
                                   const parameters& params,
                                   const arguments& args)
{
    return new BrdfGrid(size, params, args);
}

ALTA_DLL_EXPORT data* load_data(std::istream& input,
                                const arguments& args)
{
    arguments header = arguments::parse_header(input);

#ifdef __GNUC__
# warning FIXME we are not returning a "BrdfGrid"
#endif

    std::cerr << "<<Unfinished>> the data_grid loader should provide a grid object." << std::endl;
    std::cerr << "               see " << __FILE__ << " line " << __LINE__ << std::endl;
    abort();

    return load_data_from_text(input, header);
}
