/***************************************************************************

  Copyright 2005 Mitsubishi Electric Research Laboratories All Rights Reserved.

  Permission to use, copy and modify this software and its documentation
  without fee for educational, research and non-profit purposes, is hereby
  granted, provided that the above copyright notice and the following three
  paragraphs appear in all copies.

  To request permission to incorporate this software into commercial products
  contact: Vice President of Marketing and Business Development; Mitsubishi
  Electric Research Laboratories (MERL), 201 Broadway, Cambridge, MA 02139 or
  <license@merl.com>.

  IN NO EVENT SHALL MERL BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL,
  INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF
  THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF MERL HAS BEEN ADVISED
  OF THE POSSIBILITY OF SUCH DAMAGES.

  MERL SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
  PURPOSE.  THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND MERL
  HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS OR
  MODIFICATIONS.

*****************************************************************************/

#include <core/data.h>
#include <core/common.h>
#include <core/args.h>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>

#define BRDF_SAMPLING_RES_THETA_H  90
#define BRDF_SAMPLING_RES_THETA_D  90
#define BRDF_SAMPLING_RES_PHI_D    360

#define RED_SCALE   (1.00/1500.0)
#define GREEN_SCALE (1.15/1500.0)
#define BLUE_SCALE  (1.66/1500.0)

using namespace alta;


/*! \ingroup datas
 *  \class data_merl
 *  \brief Data interface for the [MERL][merl] file format.
 *   [merl]: http://people.csail.mit.edu/wojciech/BRDFDatabase/brdfs/
 *
 *  \details
 *  This plugin enables to load data measurments from the Mitsubishi Electric
 *  Research Laboratories. [This dataset][merl] contains 100 reflectance samples
 *  acquired from painted spheres using an imaging system.
 *
 *  The data is stored as RGB tuples in the \ref params::RUSIN_TH_TD_PD
 *  "RUSIN_TH_TD_PD" isotropic parametrization. A non-linear transformation of
 *  the first component is used to improve the storage of specular materials.
 *  The resolution is 90 samples for \f$ \theta_H \f$ and \f$ \theta_D \f$ and
 *  360 samples for \f$ \phi_D \f$
 *
 *  Note that this plugin is not compatible with the anisotropic measurments
 *  of [Ngan et al. [2005]][ngan].
 *
 *  \author Laurent Belcour <laurent.belcour@umontreal.ca>
 *  \author Original code from Mitsubishi Electric Research Laboratories
 *
 *   [merl]: http://people.csail.mit.edu/wojciech/BRDFDatabase/brdfs/
 *   [ngan]: http://people.csail.mit.edu/addy/research/brdf/
 */

ALTA_DLL_EXPORT data* load_data(std::istream& input, const arguments& args);

#define MERL_SIZE                                           \
    (BRDF_SAMPLING_RES_THETA_H * BRDF_SAMPLING_RES_THETA_D  \
     * BRDF_SAMPLING_RES_PHI_D / 2)

class MERL : public data
{
private: // data
	double *brdf ;
	const int _nSlice;

    MERL()
        : MERL(parameters(3, 3, params::RUSIN_TH_TD_PD, params::RGB_COLOR))
    { }

public: // methods

    MERL(const parameters& params) :
      data(parameters(3, 3, params::RUSIN_TH_TD_PD, params::RGB_COLOR), MERL_SIZE), _nSlice(MERL_SIZE) {
		brdf = new double[3*_nSlice];
		std::fill(brdf, brdf + 3*_nSlice, 0.0);

    _min.resize(3);
    _min[0] = 0.0;
    _min[1] = 0.0;
    _min[2] = 0.0;

    _max.resize(3);
    _max[0] = 0.5*M_PI;
    _max[1] = 0.5*M_PI;
    _max[2] = 2.0*M_PI;
    }

    ~MERL() {
    	delete[] brdf;
    }


	void save(const std::string& filename) const
	{
		FILE *f = fopen(filename.c_str(), "wb");

		int dims[3];
		dims[0] = BRDF_SAMPLING_RES_THETA_H;
		dims[1] = BRDF_SAMPLING_RES_THETA_D;
		dims[2] = BRDF_SAMPLING_RES_PHI_D/2;

		const int n = dims[0]*dims[1]*dims[2];

		fwrite(dims, sizeof(int), 3, f);
		fwrite(brdf, sizeof(double), 3*n, f);

		fclose(f);
	}

	// Acces to data
	vec get(int i) const
	{
		int phid_ind = i % (BRDF_SAMPLING_RES_PHI_D / 2);
		int thed_ind = (i / (BRDF_SAMPLING_RES_PHI_D / 2)) % BRDF_SAMPLING_RES_THETA_D ;
		int theh_ind = (i / ((BRDF_SAMPLING_RES_PHI_D / 2) * BRDF_SAMPLING_RES_THETA_D))
			            % BRDF_SAMPLING_RES_THETA_H ;


		vec res(6) ;
		res[2] = phi_diff_from_index(phid_ind);
		res[1] = theta_diff_from_index(thed_ind);
		res[0] = theta_half_from_index(theh_ind);
	#ifdef DEBUG
		std::cout << "get -> " << i << " (" << theh_ind << ", " << thed_ind << ", " << phid_ind << ")" << std::endl;
		std::cout << "       " << res[0] << ", " << res[1]  << ", " << res[2] << std::endl;
	#endif
		res[3] = brdf[i] * RED_SCALE;
		res[4] = brdf[i + BRDF_SAMPLING_RES_THETA_H*BRDF_SAMPLING_RES_THETA_D*BRDF_SAMPLING_RES_PHI_D/2] * GREEN_SCALE;
		res[5] = brdf[i + BRDF_SAMPLING_RES_THETA_H*BRDF_SAMPLING_RES_THETA_D*BRDF_SAMPLING_RES_PHI_D] * BLUE_SCALE;
		return res ;
	}

	void set(int i, const vec& x) {
    	assert(x.size() == (parametrization().dimX() + parametrization().dimY()));
		int iR = i;
		int iG = iR + _nSlice;
		int iB = iG + _nSlice;
		brdf[iR] = x[parametrization().dimX()+0] / RED_SCALE;
		brdf[iG] = x[parametrization().dimX()+1] / GREEN_SCALE;
		brdf[iB] = x[parametrization().dimX()+2] / BLUE_SCALE;
	}

	vec value(const vec& in) const {
	    double r, g, b;

	    lookup_brdf_val(brdf, in[0], in[1], in[2], r, g, b) ;

	    vec res(3);
	    res[0] = r ;
	    res[1] = g ;
	    res[2] = b ;

	    if( res[0] < 0.0 || res[1] < 0.0 || res[2] < 0.0 )
	    {

#ifdef DEBUG
	    	std::cout << __FILE__ << " " << __LINE__ << " in[0] = " << in[0]
	    						<< " in[1] = " << in[1] << " in[2] = " << in[2] << std::endl;
	    	std::cout <<  "res = " << res << std::endl;
#endif
	    	res[0] = 0.0;
	    	res[1] = 0.0;
	    	res[2] = 0.0;

	    	//assert(0);
	    }
	    return res;
	}


private: //methods

	// cross product of two vectors
	void cross_product (double* v1, double* v2, double* out) const
	{
		out[0] = v1[1]*v2[2] - v1[2]*v2[1];
		out[1] = v1[2]*v2[0] - v1[0]*v2[2];
		out[2] = v1[0]*v2[1] - v1[1]*v2[0];
	}

	// normalize vector
	void normalize(double* v) const
	{
		// normalize
		double len = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
		v[0] = v[0] / len;
		v[1] = v[1] / len;
		v[2] = v[2] / len;
	}

	// rotate vector along one axis
	void rotate_vector(double* vector, double* axis, double angle, double* out) const
	{
		double temp;
		double cross[3];
		double cos_ang = cos(angle);
		double sin_ang = sin(angle);

		out[0] = vector[0] * cos_ang;
		out[1] = vector[1] * cos_ang;
		out[2] = vector[2] * cos_ang;

		temp = axis[0]*vector[0]+axis[1]*vector[1]+axis[2]*vector[2];
		temp = temp*(1.0-cos_ang);

		out[0] += axis[0] * temp;
		out[1] += axis[1] * temp;
		out[2] += axis[2] * temp;

		cross_product (axis,vector,cross);

		out[0] += cross[0] * sin_ang;
		out[1] += cross[1] * sin_ang;
		out[2] += cross[2] * sin_ang;
	}


	// convert standard coordinates to half vector/difference vector coordinates
	void std_coords_to_half_diff_coords(double theta_in, double fi_in, double theta_out, double fi_out,
									double& theta_half,double& fi_half,double& theta_diff,double& fi_diff ) const
	{

		// compute in vector
		double in_vec_z = cos(theta_in);
		double proj_in_vec = sin(theta_in);
		double in_vec_x = proj_in_vec*cos(fi_in);
		double in_vec_y = proj_in_vec*sin(fi_in);
		double in[3]= {in_vec_x,in_vec_y,in_vec_z};
		normalize(in);


		// compute out vector
		double out_vec_z = cos(theta_out);
		double proj_out_vec = sin(theta_out);
		double out_vec_x = proj_out_vec*cos(fi_out);
		double out_vec_y = proj_out_vec*sin(fi_out);
		double out[3]= {out_vec_x,out_vec_y,out_vec_z};
		normalize(out);


		// compute halfway vector
		double half_x = (in_vec_x + out_vec_x)/2.0f;
		double half_y = (in_vec_y + out_vec_y)/2.0f;
		double half_z = (in_vec_z + out_vec_z)/2.0f;
		double half[3] = {half_x,half_y,half_z};
		normalize(half);

		// compute  theta_half, fi_half
		theta_half = acos(half[2]);
		fi_half = atan2(half[1], half[0]);


		double bi_normal[3] = {0.0, 1.0, 0.0};
		double normal[3] = { 0.0, 0.0, 1.0 };
		double temp[3];
		double diff[3];

		// compute diff vector
		rotate_vector(in, normal , -fi_half, temp);
		rotate_vector(temp, bi_normal, -theta_half, diff);

		// compute  theta_diff, fi_diff
		theta_diff = acos(diff[2]);
		fi_diff = atan2(diff[1], diff[0]);

	}


	// Lookup theta_half index
	// This is a non-linear mapping!
	// In:  [0 .. pi/2]
	// Out: [0 .. 89]
	inline int theta_half_index(double theta_half) const
	{
		if (theta_half <= 0.0)
			return 0;
		double theta_half_deg = ((theta_half / (M_PI/2.0))*BRDF_SAMPLING_RES_THETA_H);
		double temp = theta_half_deg*BRDF_SAMPLING_RES_THETA_H;
		temp = sqrt(temp);
		int ret_val = (int)temp;
		if (ret_val < 0) ret_val = 0;
		if (ret_val >= BRDF_SAMPLING_RES_THETA_H)
			ret_val = BRDF_SAMPLING_RES_THETA_H-1;
		return ret_val;
	}

	// Lookup theta_half from index
	// This is a non-linear mapping!
	// In:  [0 .. 89]
	// Out: [0 .. pi/2]
	inline double theta_half_from_index(int index) const
	{
	    const double temp = index+0.5;
	    const double theta_half_deg = (temp*temp) / BRDF_SAMPLING_RES_THETA_H;
		const double ret_val = theta_half_deg * (0.5*M_PI) / BRDF_SAMPLING_RES_THETA_H;
		return ret_val;
	}

	// Lookup theta_diff index
	// In:  [0 .. pi/2]
	// Out: [0 .. 89]
	inline int theta_diff_index(double theta_diff) const
	{
		int tmp = int(theta_diff / (M_PI * 0.5) * BRDF_SAMPLING_RES_THETA_D);
		if (tmp < 0)
			return 0;
		else if (tmp < BRDF_SAMPLING_RES_THETA_D - 1)
			return tmp;
		else
			return BRDF_SAMPLING_RES_THETA_D - 1;
	}

	// Lookup theta_diff from index
	// In:  [0 .. 89]
	// Out: [0 .. pi/2]
	inline double theta_diff_from_index(int index) const
	{
	    const double temp = index+0.5;
		const double theta_diff = temp * (0.5*M_PI) / BRDF_SAMPLING_RES_THETA_D;
		return theta_diff;
	}


	// Lookup phi_diff index
	inline int phi_diff_index(double phi_diff) const
	{
		// Because of reciprocity, the BRDF is unchanged under
		// phi_diff -> phi_diff + M_PI
		if (phi_diff < 0.0)
			phi_diff += M_PI;

		// In: phi_diff in [0 .. pi]
		// Out: tmp in [0 .. 179]
		int tmp = int(phi_diff / M_PI * BRDF_SAMPLING_RES_PHI_D / 2);
		if (tmp < 0)
			return 0;
		else if (tmp < BRDF_SAMPLING_RES_PHI_D / 2 - 1)
			return tmp;
		else
			return BRDF_SAMPLING_RES_PHI_D / 2 - 1;
	}

	// Lookup phi_diff from index
	//
	inline double phi_diff_from_index(int index) const
	{
	    const double phi_diff = (index+0.5) * M_PI / (BRDF_SAMPLING_RES_PHI_D/2);
		return phi_diff;
	}


	// Given a pair of incoming/outgoing angles, look up the BRDF.
	void lookup_brdf_val(double* brdf, double theta_half,
				  double theta_diff, double fi_diff,
				  double& red_val,double& green_val,double& blue_val) const
	{
		// The phi index needs to be positive.
		if(fi_diff < 0.0) {
			fi_diff = - fi_diff;
		}

		// The data is symmetric on fi_diff with respect to PI
		if(fi_diff > M_PI) {
			fi_diff = 2.0*M_PI - fi_diff;
		}

	    // Testing the input domain to avoid indexing outside of the
		 // allocated memory.

		 //ROMAIN PAC:  fi_diff IS ALWAYS >= 0.0 ACCORDING TO WHAT IS DONE BEFORE!!!
		if(theta_half < 0.0 || theta_half > 0.5*M_PI ||
		   theta_diff < 0.0 || theta_diff > 0.5*M_PI ||
		   fi_diff > M_PI || fi_diff < 0.0) {
			red_val   = 0.0;
			green_val = 0.0;
			blue_val  = 0.0;
			return;
		}

		// Find index.
		// Note that phi_half is ignored, since isotropic BRDFs are assumed
		int ind = phi_diff_index(fi_diff) +
			  theta_diff_index(theta_diff) * BRDF_SAMPLING_RES_PHI_D / 2 +
			  theta_half_index(theta_half) * BRDF_SAMPLING_RES_PHI_D / 2 *
						         BRDF_SAMPLING_RES_THETA_D;

		red_val = brdf[ind] * RED_SCALE;
		green_val = brdf[ind + BRDF_SAMPLING_RES_THETA_H*BRDF_SAMPLING_RES_THETA_D*BRDF_SAMPLING_RES_PHI_D/2] * GREEN_SCALE;
		blue_val = brdf[ind + BRDF_SAMPLING_RES_THETA_H*BRDF_SAMPLING_RES_THETA_D*BRDF_SAMPLING_RES_PHI_D] * BLUE_SCALE;

	#ifdef DEBUG
		if (red_val < 0.0 || green_val < 0.0 || blue_val < 0.0)
		{
			fprintf(stderr, "Negative value [%f, %f, %f].\n", theta_half, theta_diff, fi_diff);
	    std::cerr << " red_val = " << red_val << " green_val = " << green_val << " blue_val = " << blue_val << std::endl;
	    std::cerr << " AT INDEX = " << ind << std::endl;
		}

	#endif
	}

	// Read BRDF data
  friend data* load_data(std::istream&, const arguments&);
};

static bool read_brdf(std::istream& input, double* &brdf)
{
		int dims[3];
    input.read((char *) &dims, sizeof dims);
		int n = dims[0] * dims[1] * dims[2];
		if (n != BRDF_SAMPLING_RES_THETA_H *
        BRDF_SAMPLING_RES_THETA_D *
        BRDF_SAMPLING_RES_PHI_D / 2)
		{
        fprintf(stderr, "Dimensions don't match\n");
        return false;
		}

    input.read((char *) brdf, 3 * n * sizeof(double));

		return true;
}


ALTA_DLL_EXPORT data* provide_data(size_t size, const parameters& params,
                                   const arguments& args)
{
    return new MERL(params);
}

ALTA_DLL_EXPORT data* load_data(std::istream& input, const arguments& args)
{
    MERL* result = new MERL();

    if(!read_brdf(input, result->brdf))
		{
        std::cerr << "<<ERROR>> unable to load the data as a MERL file" << std::endl ;
        throw;
		}

    return result;
}
