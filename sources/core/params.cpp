/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2014 CNRS
   Copyright (C) 2013, 2014, 2015, 2016 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include "params.h"
#include "common.h"
#include <cassert>

using namespace alta;

struct param_info
{
	param_info(std::string n, int d, std::string i) :
        name(n), dimension(d), info(i) { }

	std::string name;
	int dimension;
	std::string info;
};

std::map<params::input, const param_info> create_map()
{
	std::map<params::input, const param_info> _map;
	/* 1D Params */
	_map.insert(std::make_pair<params::input, const param_info>(params::COS_TH, param_info("COS_TH", 1, "Cosine of the Half angle")));
	_map.insert(std::make_pair<params::input, const param_info>(params::COS_TK, param_info("COS_TK", 1, "Cosine of the Back angle")));
	_map.insert(std::make_pair<params::input, const param_info>(params::COS_TLV, param_info("COS_TLV", 1, "Cosine of the Light and View directions")));
	_map.insert(std::make_pair<params::input, const param_info>(params::COS_TLR, param_info("COS_TLR", 1, "Cosine of the Light and Reflected directions")));

	/* 2D Params */
	_map.insert(std::make_pair<params::input, const param_info>(params::RUSIN_TH_TD, param_info("RUSIN_TH_TD", 2, "Radialy symmetric Half angle parametrization")));
	_map.insert(std::make_pair<params::input, const param_info>(params::COS_TH_TD, param_info("COS_TH_TD", 2, "Cosines of the elevation angles of the Half angle parametrization")));
	_map.insert(std::make_pair<params::input, const param_info>(params::ISOTROPIC_TV_PROJ_DPHI, param_info("ISOTROPIC_TV_PROJ_DPHI", 2, "Isoptropic projected phi parametrization without a light direction.")));
	_map.insert(std::make_pair<params::input, const param_info>(params::STARK_2D, param_info("STARK_2D", 2, "Stark parametrization H, B but without third component.")));
	_map.insert(std::make_pair<params::input, const param_info>(params::NEUMANN_2D, param_info("NEUMANN_2D", 2, "Neumann parametrization H, B but without third component.")));

	/* 3D Params */
	_map.insert(std::make_pair<params::input, const param_info>(params::RUSIN_TH_TD_PD, param_info("RUSIN_TH_TD_PD", 3, "Isotropic Half angle parametrization")));
	_map.insert(std::make_pair<params::input, const param_info>(params::ISOTROPIC_TV_TL_DPHI, param_info("ISOTROPIC_TV_TL_DPHI", 3, "Isotropic Light/View angle parametrization")));
	_map.insert(std::make_pair<params::input, const param_info>(params::RUSIN_VH, param_info("RUSIN_VH", 3, "Vector representation of the Half angle only")));
	_map.insert(std::make_pair<params::input, const param_info>(params::SCHLICK_VK, param_info("SCHLICK_VH", 3, "Vector representation of the Back angle only")));
	_map.insert(std::make_pair<params::input, const param_info>(params::ISOTROPIC_TL_TV_PROJ_DPHI, param_info("ISOTROPIC_TL_TV_PROJ_DPHI", 3, "Isoptropic projected phi parametrization.")));
	_map.insert(std::make_pair<params::input, const param_info>(params::SCHLICK_TL_TK_PROJ_DPHI, param_info("SCHLICK_TL_TK_PROJ_DPHI", 3, "Isoptropic projected phi parametrization centered around the back vector.")));
	_map.insert(std::make_pair<params::input, const param_info>(params::RETRO_TL_TVL_PROJ_DPHI, param_info("RETRO_TL_TVL_PROJ_DPHI", 3, "Isoptropic retro projected phi parametrization.")));
	_map.insert(std::make_pair<params::input, const param_info>(params::STARK_3D, param_info("STARK_3D", 3, "Stark parametrization H, B.")));
	_map.insert(std::make_pair<params::input, const param_info>(params::NEUMANN_3D, param_info("NEUMANN_3D", 3, "Neumann parametrization H, B.")));

	/* 4D Params */
	_map.insert(std::make_pair<params::input, const param_info>(params::RUSIN_TH_PH_TD_PD, param_info("RUSIN_TH_PH_TD_PD", 4, "Complete Half angle parametrization")));
	_map.insert(std::make_pair<params::input, const param_info>(params::SPHERICAL_TL_PL_TV_PV, param_info("SPHERICAL_TL_PL_TV_PV", 4, "Complete classical parametrization")));
	_map.insert(std::make_pair<params::input, const param_info>(params::STEREOGRAPHIC, param_info("STEREOGRAPHIC", 4, "Light/View vector in stereographic projection"))),

	/* 6D Param */
	_map.insert(std::make_pair<params::input, const param_info>(params::CARTESIAN, param_info("CARTESIAN", 6, "Complete vector parametrization")));

	return _map;
}
static const std::map<params::input, const param_info> input_map = create_map();


static std::map<params::output, std::string> create_output_map()
{
		std::map<params::output, std::string> result;

#define STRINGIFY_(x) #x
#define STRINGIFY(x)  STRINGIFY_(x)
#define DEFINE_MAPPING(name)										\
		result[params::name] = STRINGIFY(name);

		DEFINE_MAPPING(INV_STERADIAN);
		DEFINE_MAPPING(INV_STERADIAN_COSINE_FACTOR);
		DEFINE_MAPPING(ENERGY);
		DEFINE_MAPPING(RGB_COLOR);
		DEFINE_MAPPING(XYZ_COLOR);

#undef DEFINE_MAPPING
#undef STRINGIFY
#undef STRINGIFY_

		return result;
}

static const std::map<params::output, std::string> output_map =
		create_output_map();

void params::to_cartesian(const double* invec, params::input intype,
		double* outvec)
{
	#ifdef DEBUG_PARAM
	std::cout << " ENTERING to_cartesian with inttype = " << get_name(intype) << std::endl;
	#endif

	switch(intype)
	{
		// 1D Parametrizations
		case params::COS_TH:
			half_to_cartesian(acos(invec[0]), 0.0, 0.0, 0.0, outvec);
			break;
		case params::COS_TK:
			outvec[0] = invec[0]*invec[0]-1.0;
			outvec[1] = 0.0;
			outvec[2] = invec[0];
			outvec[3] = 1.0-invec[0]*invec[0];
			outvec[4] = 0.0;
			outvec[5] = invec[0];
			break;
		case params::COS_TLV:
			outvec[0] = sqrt(1.0 - invec[0]*invec[0]);
			outvec[1] = 0.0;
			outvec[2] = invec[0];
			outvec[3] = 0.0;
			outvec[4] = 0.0;
			outvec[5] = 1.0;
			break;
		case params::COS_TLR:
			outvec[0] = - sqrt(1.0 - invec[0]*invec[0]);
			outvec[1] = 0.0;
			outvec[2] = invec[0];
			outvec[3] = 0.0;
			outvec[4] = 0.0;
			outvec[5] = 1.0;
			break;

			// 2D Parametrizations
		case params::COS_TH_TD:
			half_to_cartesian(acos(invec[0]), 0.0, acos(invec[1]), 0.0, outvec);
			break;
		case params::RUSIN_TH_TD:
         half_to_cartesian(invec[0], 0.0, invec[1], 0.5*M_PI, outvec);
			break;
		case ISOTROPIC_TV_PROJ_DPHI:
		{
			const double theta = sqrt(invec[0]*invec[0] + invec[1]*invec[1]);
			if(theta > 0.0)
			{
				outvec[0] = (invec[0]/theta)*sin(theta);
				outvec[1] = (invec[1]/theta)*sin(theta);
			}
			else
			{
				outvec[0] = 0.0;
				outvec[1] = 0.0;
			}
			outvec[2] = cos(theta);
			outvec[3] = 0.0;
			outvec[4] = 0.0;
			outvec[5] = 1.0;
		}
			break;
		// invec[0] = ||Hp|| Norm of the projected normalized Half vector
		// H=(L+V)/2 invec[1] = ||B||  Norm of the unormalized Back vector
		// B=(L-V)/2
		case STARK_2D:
		{
			const double Hx = invec[0];
			const double Hy = 0;
         const double sum = Hx*Hx + invec[1]*invec[1];
			const double Hz = (sum <= 1.0) ? sqrt(1.0 - std::min(sum, 1.0)) : -1.0;

			// Ensuring that <H,B> = 0.
			const double Bx = 0.0;
			const double By = invec[1];
			const double Bz = 0.0;

         //assert(Hz >= 0 && Hz <= 1);
			outvec[0] = Hx-Bx;
			outvec[1] = Hy-By;
			outvec[2] = Hz-Bz;
			outvec[3] = Hx+Bx;
			outvec[4] = Hy+By;
			outvec[5] = Hz+Bz;
		}
		break;
		case NEUMANN_2D:
		{
			outvec[0] =   invec[0];
			outvec[1] =   invec[1];
			outvec[2] =   sqrt(1.0 - outvec[0]*outvec[0] - outvec[1]*outvec[1]);
			outvec[3] =   invec[0];
			outvec[4] = - invec[1];
			outvec[5] =   sqrt(1.0 - outvec[3]*outvec[3] - outvec[4]*outvec[4]);
		}
		break;


			// 3D Parametrization
		case params::RUSIN_TH_PH_TD:
			half_to_cartesian(invec[0], invec[1], invec[2], 0.0, outvec);
			break;
		case params::RUSIN_TH_TD_PD:
			half_to_cartesian(invec[0], 0.0, invec[1], invec[2], outvec);
#ifdef DEBUG
			std::cout << outvec[2] << std::endl;
#endif
			break;
		case params::ISOTROPIC_TV_TL_DPHI:
			classical_to_cartesian(invec[0], 0.0, invec[1], invec[2], outvec);
			break;
		case params::RUSIN_VH:
            half_to_cartesian(acos(invec[2]), atan2(invec[1], invec[0]), 0.0, 0.0, outvec);
			break;
        // \todo I should handle the phi_k in the conversion to CARTESIAN
      case params::SCHLICK_VK:
            outvec[0] = invec[2]*invec[2]-1.0;
            outvec[1] = 0.0;
            outvec[2] = invec[2];
            outvec[3] = 1.0-invec[2]*invec[2];
            outvec[4] = 0.0;
            outvec[5] = invec[2];
            break;
		case ISOTROPIC_TL_TV_PROJ_DPHI:
		{
			const double theta = sqrt(invec[1]*invec[1] + invec[2]*invec[2]);
			// theta can equals zero. In that case, you should not
			// divide by its value, and the relative orientation does
			// not matter (normal direction).
			if(theta > 0.0)
			{
				outvec[0] = (invec[1]/theta)*sin(theta);
				outvec[1] = (invec[2]/theta)*sin(theta);
			}
			else
			{
				outvec[0] = 0.0;
				outvec[1] = 0.0;
			}
			outvec[2] = cos(theta);
			outvec[3] = sin(invec[0]);
			outvec[4] = 0.0;
			outvec[5] = cos(invec[0]);
		}
			break;
		case SCHLICK_TL_TK_PROJ_DPHI:
		{
			// Set the light direction
			outvec[3] = sin(invec[0]);
			outvec[4] = 0.0;
			outvec[5] = cos(invec[0]);

			// The view direction is the symmetric of the reflected direction
			// with respect to the back direction:
			//     v = 2 |r.k| k - r
			//     r = 2 |l.n| n - l = {-lx, -ly, lz }
			const double theta = sqrt(invec[1]*invec[1] + invec[2]*invec[2]);
			const double Kx = (theta > 0.0) ? (invec[1]/theta)*sin(theta) : 0.0;
			const double Ky = (theta > 0.0) ? (invec[2]/theta)*sin(theta) : 0.0;
			const double Kz = cos(theta);
			const double dotKR = outvec[5]*Kz - outvec[3]*Kx;
			outvec[0] = 2.0*dotKR*Kx + outvec[3];
			outvec[1] = 2.0*dotKR*Ky;
			outvec[2] = 2.0*dotKR*Kz - outvec[5];
		}
			break;
		case RETRO_TL_TVL_PROJ_DPHI:
		{
			const double theta_vl = sqrt(invec[1]*invec[1] + invec[2]*invec[2]);
			const double cos_vl   = cos(theta_vl);
			const double sin_vl   = sin(theta_vl);
			const double cos_l    = cos(invec[0]);
			const double sin_l    = sin(invec[0]);
			const double cos_phi  = (theta_vl > 0.0) ? invec[1] / theta_vl : 0.0;
			const double sin_phi  = (theta_vl > 0.0) ? invec[2] / theta_vl : 0.0;

			// Compute the cosine of the outgoing vector using the
			// spherical law of cosines
			const double cos_v = cos_l*cos_vl + cos_phi*sin_l*sin_vl;
			const double sin_v = sqrt(1.0 - std::min(cos_v*cos_v, 1.0));

			if(sin_v != 0.0)
			{
				const double sin_dphi = (sin_vl*sin_phi) / sin_v;

				outvec[0] = sin_v * sqrt(1.0 - std::min(sin_dphi*sin_dphi, 1.0));
				outvec[1] = sin_v * sin_dphi;
			}
			else
			{
				outvec[0] = 0.0;
				outvec[1] = 0.0;
			}
			outvec[2] = cos_v;
			outvec[3] = sin(invec[0]);
			outvec[4] = 0.0;
			outvec[5] = cos(invec[0]);
		}
			break;
		case STARK_3D:
		{
			// Constructing the Half vector
			const double Hx = invec[0];
			const double Hy = 0;
			const double Hz = sqrt(1.0 - invec[0]*invec[0] - invec[1]*invec[1]);

			// Constructing the Back vector using azimuth and elevation angles
			const double cosPhi = cos(invec[2]);
			const double sinPhi = sin(invec[2]);

#define STARK_BUILDING_THETA

#ifdef STARK_BUILDING_THETA
			const double Theta  = atan2(-Hz, invec[0]*cosPhi);
			const double cosThe = cos(Theta);
			const double sinThe = sin(Theta);
#else
			// Pascal Barla's way to do it
			const double cosThe = -cosPhi*invec[0] / sqrt(Hz*Hz + cosPhi*cosPhi*invec[0]*invec[0]);
			const double sinThe = Hz / (sqrt(Hz*Hz + cosPhi*cosPhi * invec[0]*invec[0]));
#endif
			const double Bx = invec[1]*sinThe*cosPhi;
			const double By = invec[1]*sinThe*sinPhi;
			const double Bz = invec[1]*cosThe;

			outvec[0] = Hx-Bx;
			outvec[1] = Hy-By;
			outvec[2] = Hz-Bz;
			outvec[3] = Hx+Bx;
			outvec[4] = Hy+By;
			outvec[5] = Hz+Bz;
		}
			break;
		case NEUMANN_3D:
		{
			const double cosPhi = cos(invec[2]);
			const double sinPhi = sin(invec[2]);

			// Build the projected Half and Back vectors
			outvec[0] =   invec[0] + cosPhi*invec[1];
			outvec[1] =   sinPhi*invec[1];
			outvec[3] =   invec[0] - cosPhi*invec[1];
			outvec[4] = - sinPhi*invec[1];

			// Safeguard, if the vectors are not under unit length return an
			// invalid configuration
			if(outvec[0]*outvec[0]+outvec[1]*outvec[1] > 1.0 ||
				outvec[3]*outvec[3]+outvec[4]*outvec[4] > 1.0) {
				outvec[0] =  0.0;
				outvec[1] =  0.0;
				outvec[2] = -1.0;
				outvec[3] =  0.0;
				outvec[4] =  0.0;
				outvec[5] = -1.0;
				break;
			}

			// Project the vectorx on the hemisphere.
			outvec[2] =   sqrt(1.0 - outvec[0]*outvec[0] - outvec[1]*outvec[1]);
			outvec[5] =   sqrt(1.0 - outvec[3]*outvec[3] - outvec[4]*outvec[4]);
		}
			break;

			// 4D Parametrization
		case params::RUSIN_TH_PH_TD_PD:
			half_to_cartesian(invec[0], invec[1], invec[2], invec[3], outvec);
			break;

		case params::SPHERICAL_TL_PL_TV_PV:
			outvec[0] = cos(invec[3])*sin(invec[2]);
			outvec[1] = sin(invec[3])*sin(invec[2]);
			outvec[2] = cos(invec[2]);
			outvec[3] = cos(invec[1])*sin(invec[0]);
			outvec[4] = sin(invec[1])*sin(invec[0]);
			outvec[5] = cos(invec[0]);
			break;

		case params::STEREOGRAPHIC:
      {
			// Project the 2D direction on the surface invec[0,1] to the View
         // vector outvec[0,1,2]
         const double normL = invec[0]*invec[0] + invec[1]*invec[1];
         const double alphL = (1.0 - normL) / (1.0 + normL);
         outvec[0] = invec[0] + alphL;
         outvec[1] = invec[1] + alphL;
         outvec[2] = alphL;

         // Project the 2D direction on the surface invec[2,3] to the Light
         // vector outvec[3,4,5]
         const double normV = invec[2]*invec[2] + invec[3]*invec[3];
         const double alphV = (1.0 - normV) / (1.0 + normV);
         outvec[3] = invec[2] + alphV;
         outvec[4] = invec[3] + alphV;
         outvec[5] = alphV;
		}
         break;

		// 6D Parametrization
      case params::CARTESIAN:
			memcpy(outvec, invec, 6*sizeof(double));
			break;

		default:
			std::cerr << "<<ERROR>> Transformation not implemented, " << get_name(intype) << " " << __FILE__ << ":" << __LINE__ << std::endl;
			throw;
			break;
	}

}

// Return in HALF the half vector for the two vectors in INVEC.
template<class Vector, typename Real>
static void half_vector(Real const *invec, Vector &half)
{
    half[0] = invec[0] + invec[3];
    half[1] = invec[1] + invec[4];
    half[2] = invec[2] + invec[5];
    Real sqnorm = half[0]*half[0] + half[1]*half[1] + half[2]*half[2];

    if (sqnorm <= 0.) {
        half[0] = 0;
        half[1] = 0;
        half[2] = 1;
    }
    else {
        half /= sqrt(sqnorm);
    }

    assert(half[2] <= 1.);
    assert(half[2] >= 0.);
}

// Return true if X is essentially zero.  We have to resort to this hack
// before using functions that have a high derivative around zero, such as
// 'atan2'.
template<typename Real>
static bool is_zero(Real x)
{
    return close_to(x, 0., 1e-14);
}

// Return true if VEC is alongside z⃗.
template<class Vector>
static bool alongside_z(const Vector& vec)
{
    return (is_zero(vec[0]) && is_zero(vec[1]))
        || close_to(vec[2], 1., 1e-14);
}


void params::from_cartesian(const double* invec, params::input outtype,
		double* outvec)
{
	#ifdef DEBUG_PARAM
	std::cout << " ENTERING from_cartesian with outtype = " << get_name(outtype) << std::endl;
	std::cout << " invec = " << invec[0] <<  " " << invec[1] << " " << invec[2] << " " << invec[3]
						<< "  " << invec[4] << " " << invec[5] << std::endl;
	#endif

	// Compute the half vector.
  Eigen::Vector3d half;
  half_vector(invec, half);

	// Difference vector.
  Eigen::Vector3d diff;

	switch(outtype)
	{
		// 1D Parametrizations
		case params::COS_TH:
			outvec[0] = half[2];
			break;
		case params::COS_TK:
		{
			const double Kx = invec[0]-invec[3];
			const double Ky = invec[1]-invec[4];
			const double Kz = invec[2]+invec[5];
			outvec[0] = (invec[2] + invec[5]) / sqrt(Kx*Kx + Ky*Ky + Kz*Kz);
		}
			break;
		case params::COS_TLV:
			outvec[0] = invec[0]*invec[3] + invec[1]*invec[4] + invec[2]*invec[5];
			break;
		case params::COS_TLR:
			outvec[0] = invec[0]*invec[3] - (invec[1]*invec[4] + invec[2]*invec[5]);
			break;

			// 2D Parametrizations
		case params::COS_TH_TD:
			outvec[0] = half[2];
			outvec[1] = half[0]*invec[0] + half[1]*invec[1] + half[2]*invec[2];
			break;
		case params::RUSIN_TH_TD:
			outvec[0] = acos(half[2]);
			outvec[1] = acos(clamp(half[0]*invec[0] + half[1]*invec[1] + half[2]*invec[2], -1.0, 1.0));
			break;
		case ISOTROPIC_TV_PROJ_DPHI:
		{
			const double theta = acos(invec[2]);
			const double dphi  = atan2(invec[1], invec[0]) - atan2(invec[4], invec[3]);
			outvec[0] = theta * cos(dphi);
			outvec[1] = theta * sin(dphi);
		}
			break;
		// outvec[0] = ||Hp|| Norm of the projected unormalized Half vector
		//             H = (V+L)/2
		// outvec[1] = ||B|| Norm of the unormalized Back vector B = (L-V)/2
		case STARK_2D:
		{
			double Hx = 0.5*(invec[0]+invec[3]);
			double Hy = 0.5*(invec[1]+invec[4]);
			outvec[0] = sqrt(Hx*Hx + Hy*Hy);

			double Bx = 0.5*(invec[3]-invec[0]);
			double By = 0.5*(invec[4]-invec[1]);
			double Bz = 0.5*(invec[5]-invec[2]);
			outvec[1] = sqrt(Bx*Bx + By*By + Bz*Bz);
		}
			break;
		case NEUMANN_2D:
		{
				double Hx = 0.5*(invec[0]+invec[3]);
				double Hy = 0.5*(invec[1]+invec[4]);
				outvec[0] = sqrt(Hx*Hx + Hy*Hy);

				double Bx = 0.5*(invec[3]-invec[0]);
				double By = 0.5*(invec[4]-invec[1]);
				outvec[1] = sqrt(Bx*Bx + By*By);
		}
			break;


			// 3D Parametrization
		case params::RUSIN_TH_PH_TD:
			outvec[0] = acos(half[2]);
			outvec[1] = atan2(half[1], half[0]);
			outvec[2] = acos(half[0]*invec[0] + half[1]*invec[1] + half[2]*invec[2]);
			break;
		case params::RUSIN_TH_TD_PD:
			outvec[0] = acos(half[2]);

			// Compute the diff vector
			diff[0] = invec[0];
			diff[1] = invec[1];
			diff[2] = invec[2];

      rotate(diff, -atan2(half[1], half[0]), -outvec[0]);

      // Note: the following approach is numerically more stable (up to
      // 1e-12) but slower.
      // {
      //     vec view(3);
      //     view[0] = invec[0];
      //     view[1] = invec[1];
      //     view[2] = invec[2];
      //     outvec[1] = acos(clamp(half.dot(view),0.,1.));
      // }
      outvec[1] = acos(diff[2]);

      // By convention, when DIFF is alongside z⃗, return φd = 0.
      if (alongside_z(diff))
          outvec[2] = 0;
      else
          outvec[2] = atan2(diff[1], diff[0]);

			break;
		case params::ISOTROPIC_TV_TL_DPHI:
			outvec[0] = acos(invec[2]);
			outvec[1] = acos(invec[5]);
			outvec[2] = atan2(invec[4], invec[3]) - atan2(invec[1], invec[0]);
			break;
		case params::RUSIN_VH:
			outvec[0] = half[0];
			outvec[1] = half[1];
			outvec[2] = half[2];
			break;
      case params::SCHLICK_VK:
      {
			const double Kx = invec[3]-invec[0];
         const double Ky = invec[4]-invec[1];
         const double Kz = invec[5]+invec[2];

         const double norm =  sqrt(Kx*Kx + Ky*Ky + Kz*Kz);
			if(norm > 1.0E-10)
			{
				outvec[0] = Kx / norm;
		      outvec[1] = Ky / norm;
			   outvec[2] = Kz / norm;
			}
			else
			{
				outvec[0] = 0.0;
		      outvec[1] = 0.0;
			   outvec[2] = 1.0;
			}
      }
			break;
		case ISOTROPIC_TL_TV_PROJ_DPHI:
		{
			const double theta_l = acos(invec[5]);
			const double theta_v = acos(invec[2]);
			const double dphi    = atan2(invec[4], invec[3]) - atan2(invec[1], invec[0]);
			outvec[0] = theta_l;
			outvec[1] = theta_v * cos(dphi);
			outvec[2] = theta_v * sin(dphi);
		}
			break;
		case RETRO_TL_TVL_PROJ_DPHI:
		{
			const double cos_v  = invec[2];
			const double cos_l  = invec[5];
			const double cos_vl = invec[0]*invec[3] + invec[1]*invec[4] + invec[2]*invec[5];
			const double sin_l  = sqrt(1.0 - cos_l * cos_l);
			const double sin_v  = sqrt(1.0 - cos_v * cos_v);
			const double sin_vl = sqrt(1.0 - cos_vl* cos_vl);

			// Sine of the \Delta_{\phi} angle
			const double sin_dp = sin(atan2(invec[4], invec[3]) - atan2(invec[1], invec[0]));

			if(sin_vl != 0.0 && sin_l != 0.0)
			{
				// Cosine and sine of the phi_{L,V} angle: the spherical
				// angle between the arclength {L,N} and {V,N}
				const double cos_phi = (cos_v - cos_l*cos_vl) / (sin_l * sin_vl);
				const double sin_phi = (sin_dp * sin_v) / sin_vl;

				const double theta = acos(cos_vl);

				outvec[0] = acos(cos_l);
				outvec[1] = theta * cos_phi;
				outvec[2] = theta * sin_phi;
			}
			else
			{
				outvec[0] = acos(cos_l);
				outvec[1] = acos(cos_v) - acos(cos_l);
				outvec[2] = 0.0;
			}
		}
			break;
		case SCHLICK_TL_TK_PROJ_DPHI:
		{
			const double vkx     = invec[0]-invec[3];
			const double vky     = invec[1]-invec[4];
			const double vkz     = invec[2]+invec[5];
			const double norm    = sqrt(vkx*vkx + vky*vky + vkz*vkz);
			const double theta_k = acos(vkz / norm);
			const double dphi    = atan2(invec[4], invec[3]) - atan2(vky/norm, vkx/norm);
			outvec[0] = acos(invec[5]);
			outvec[1] = theta_k * cos(dphi);
			outvec[2] = theta_k * sin(dphi);
		}
			break;
		case STARK_3D:
		{
			double Hx = 0.5*(invec[0]+invec[3]);
			double Hy = 0.5*(invec[1]+invec[4]);
			outvec[0] = sqrt(Hx*Hx + Hy*Hy);

			double Bx = 0.5*(invec[3]-invec[0]);
			double By = 0.5*(invec[4]-invec[1]);
			double Bz = 0.5*(invec[5]-invec[2]);
			outvec[1] = sqrt(Bx*Bx + By*By + Bz*Bz);
			outvec[2] = atan2(By, Bx) -  atan2(Hy, Hx);
         if(outvec[2] < 0.0) {
            outvec[2] = 2.0*M_PI + outvec[2];
         }
		}
			break;
		case NEUMANN_3D:
		{
				double Hx = 0.5*(invec[0]+invec[3]);
				double Hy = 0.5*(invec[1]+invec[4]);
				outvec[0] = sqrt(Hx*Hx + Hy*Hy);

				double Bx = 0.5*(invec[3]-invec[0]);
				double By = 0.5*(invec[4]-invec[1]);
				outvec[1] = sqrt(Bx*Bx + By*By);
				outvec[2] = atan2(By, Bx) - atan2(Hy, Hx);
		}
			break;

			// 4D Parametrization
		case params::RUSIN_TH_PH_TD_PD:
			outvec[0] = acos(half[2]);
			outvec[1] = atan2(half[1], half[0]);

			// Compute the diff vector
			diff[0] = invec[0];
			diff[1] = invec[1];
			diff[2] = invec[2];
      rotate(diff, -outvec[1], -outvec[0]);

      // By convention, when DIFF is alongside z⃗, return φd = θd = 0.
      if (alongside_z(diff))
      {
          outvec[2] = 0.;
          outvec[3] = 0.;
      }
      else
      {
          outvec[2] = acos(diff[2]);
          outvec[3] = atan2(diff[1], diff[0]);
      }
			break;

		case params::SPHERICAL_TL_PL_TV_PV:
			outvec[0] = acos(invec[5]);
			outvec[1] = atan2(invec[4], invec[3]);
			outvec[2] = acos(invec[2]);
			outvec[3] = atan2(invec[1], invec[0]);
#ifdef DEBUG
			std::cout << invec[2] << " - acos -> " << outvec[0] << std::endl;
#endif
			break;

      case params::STEREOGRAPHIC:
      {
			// Project the View vector invec[0,1,2] on a 2D direction on the
         // surface outvec[0,1]
         const double dotVN = invec[2];
         outvec[0] = invec[0] / (1.0+dotVN);
         outvec[1] = invec[1] / (1.0+dotVN);

         // Project the Light vector invec[0,1,2] on a 2D direction on the
         // surface outvec[2,3]
         const double dotLN = invec[5];
         outvec[2] = invec[3] / (1.0+dotLN);
         outvec[3] = invec[4] / (1.0+dotLN);

		}
         break;

			// 6D Parametrization
		case params::CARTESIAN:
			memcpy(outvec, invec, 6*sizeof(double));
			break;

		default:
			std::cerr << "<<ERROR>> Transformation not implemented, n°" << outtype << ", " << __FILE__ << ":" << __LINE__ << std::endl;
			assert(false);
			break;
	}
}
params::input params::parse_input(const std::string& txt)
{

	for(std::map<params::input, const param_info>::const_iterator it=input_map.begin(); it != input_map.end(); ++it)
	{
		if(txt.compare(it->second.name) == 0)
		{
			std::cout << "<<INFO>> parsed input parametrization " << it->second.name << " from name \"" << txt << "\"" << std::endl;
			return it->first;
		}
	}

	std::cout << "<<INFO>> the input parametrization is UNKNOWN_INPUT" << std::endl;
	return params::UNKNOWN_INPUT;
}

params::output params::parse_output(const std::string& txt)
{
	if(txt == std::string("ENERGY"))
	{
		return params::ENERGY;
	}
	else if(txt == std::string("INV_STERADIAN"))
	{
		return params::INV_STERADIAN;
	}
	else if(txt == std::string("INV_STERADIAN_COSINE_FACTOR"))
	{
		return params::INV_STERADIAN_COSINE_FACTOR;
	}
	else if(txt == std::string("RGB_COLOR"))
	{
		return params::RGB_COLOR;
	}
	else if(txt == std::string("XYZ_COLOR"))
	{
		return params::XYZ_COLOR;
	}
	else
	{
		return params::UNKNOWN_OUTPUT;
	}
}

const std::string& params::get_name(const params::input param)
{
	std::map<params::input, const param_info>::const_iterator it = input_map.find(param);
	if(it != input_map.end())
	{
		return it->second.name;
	}

#ifdef DEBUG
	std::cerr << "<<WARNING>> Unknown parametrization, n°" << param << ", "<< __FILE__ << ":" << __LINE__ << std::endl;
#endif

	static const std::string unknown = "UNKNOWN_INPUT";
	return unknown;
}

const std::string& params::get_name(const params::output param)
{
		std::map<params::output, std::string>::const_iterator it = output_map.find(param);
		if (it != output_map.end())
		{
				return it->second;
		}

		static const std::string unknown = "UNKNOWN_OUTPUT";
		return unknown;
}

int  params::dimension(params::input t)
{
	std::map<params::input, const param_info>::const_iterator it = input_map.find(t);
	if(it != input_map.end())
	{
		return it->second.dimension;
	}
	else
	{
		return -1;
	}
}

void params::print_input_params()
{
	std::cout << "List of available input parametrizations in the ALTA library:" << std::endl;
	for(std::map<params::input, const param_info>::const_iterator it=input_map.begin(); it != input_map.end(); ++it)
	{
		std::cout << "  - " << it->second.name;
		for(int i=it->second.name.size(); i<26; ++i) { std::cout << " "; }
		std::cout << it->second.info << std::endl;
	}
}


bool 
params::is_above_hemisphere( double* const invec, params::input in_param_type )
{
	vec input_in_cartesian(6);

	params::convert(invec, in_param_type, params::CARTESIAN, &input_in_cartesian[0]);

	return input_in_cartesian[2] >= 0.0  && input_in_cartesian[5] >= 0.0;

}






