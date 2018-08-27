/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2017 Unity Technologies

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */
#pragma once

// STL includes
#include <stdexcept>
#include <cassert>

// TinyEXR includes
#define TINYEXR_IMPLEMENTATION
#include <tinyexr/tinyexr.h>


/*! \class EXR_IO
 *
 * \details
 * EXR_IO provides static method to either load an RGB EXR file into a floatting
 * point array or save a RGB floatting point array to disk._data
 *
 * \author Laurent Belcour
 */
template<typename FType>
class t_EXR_IO
{
	public:
		/*! \brief Load an EXR file from an input file stream.
		 */
		static bool LoadEXR(std::istream& input, int& W,int& H, FType *& pix,int nC=3)
		{
			/* Convert the input std::istream into an unsigned char array */
			std::vector<unsigned char> _memory;
			size_t _n = 0, _m = 0;
			while(input.good()) {
				if(_n == _m) {
					_m += 256;
					_memory.resize(_m);
				}

				_memory[_n++] = (unsigned char)input.get();
			}

			/* Check the different part of loading */
			int _r = -1;
			const char* _err;

			/* Load the EXR version using TinyEXR */
			EXRVersion _version;
			_r = ParseEXRVersionFromMemory(&_version, &_memory[0], _memory.size());
			if(_r != TINYEXR_SUCCESS) {
				std::cerr << "<<ERROR>> Could not load EXR version from stream" << std::endl;
				return false;
			} else {
				std::cout << "<<DEBUG>> Loading a v" << _version.version << " EXR file" << std::endl;
			}

			/* Load the EXR header */
			EXRHeader _header;
  			InitEXRHeader(&_header);
			_r = ParseEXRHeaderFromMemory(&_header, &_version, &_memory[0], _memory.size(), &_err);
			if(_r != TINYEXR_SUCCESS) {
				std::cerr << "<<ERROR>> Could not load EXR header from stream" << std::endl;
				std::cerr << "<<ERROR>> " << _err << std::endl;
				return false;
			} else {
				std::cout << "<<DEBUG>> Loading a " << _header.num_channels << " channels EXR file" << std::endl;
			}

			/* Set the requested pixel types  */
			for(int i=0; i<_header.num_channels; ++i) {
				_header.requested_pixel_types[i] = TINYEXR_PIXELTYPE_FLOAT;
			}

			/* Load the EXR image */
			EXRImage _image;
  			InitEXRImage(&_image);

			_r = LoadEXRImageFromMemory(&_image, &_header, &_memory[0], _memory.size(), &_err);
			if (_r != TINYEXR_SUCCESS) {
				std::cerr << "<<ERROR>> Could not load EXR image from stream" << std::endl;
				std::cerr << "<<ERROR>> " << _err << std::endl;
				return false;
			} else {
				std::cout << "<<DEBUG>> Loading a " << _image.width << "x" << _image.height << " EXR file" << std::endl;
			}

			/*! \todo Check the different channels if they match RGB */

			/*! \todo  handle VS data using multi-channel */

			/* Recopy the image into the provided pixel array */
			W = _image.width;
			H = _image.height;
			pix = new FType[W*H*3];
			for(int k=0; k<_image.num_channels; ++k) {
				for(int i=0; i<W*H; ++i) {
					const FType val = (FType)(reinterpret_cast<float **>(_image.images)[k][i]);
					switch(_header.channels[k].name[0]) {
						case 'R':
							pix[3*i + 0] = val;
							break;
						case 'G':
							pix[3*i + 1] = val;
							break;
						case 'B':
							pix[3*i + 2] = val;
							break;
						default:
							std::cout << "<<DEBUG>> Unknow EXR channel \'" << _header.channels[k].name[0] << "\'" << std::endl;
					}
				}
			}

			/* Free TinyEXR memory */
			_r = FreeEXRHeader(&_header);
			_r = FreeEXRImage(&_image);
			return true ;
		}

		/*! \brief Save a RGB image into and OpenEXR file using TinyEXR.
		 *  This code uses TinyEXR's reference implementation of saving a file as
		 *  we do not manipulate streams here.
		 */
		static bool SaveEXR(const char *filename,int W,int H, const FType *pix,int nC=3)
		{
			EXRHeader header;
			InitEXRHeader(&header);

			EXRImage image;
			InitEXRImage(&image);

			image.num_channels = 3;

			std::vector<float> images[3];
			images[0].resize(W*H);
			images[1].resize(W*H);
			images[2].resize(W*H);

			// Split RGBRGBRGB... into R, G and B layer
			for (int i = 0; i < W*H; i++) {
				images[0][i] = (float)pix[3*i+0];
				images[1][i] = (float)pix[3*i+1];
				images[2][i] = (float)pix[3*i+2];
			}

			float* image_ptr[3];
			image_ptr[0] = &(images[2].at(0)); // B
			image_ptr[1] = &(images[1].at(0)); // G
			image_ptr[2] = &(images[0].at(0)); // R

			image.images = (unsigned char**)image_ptr;
			image.width  = W;
			image.height = H;

			header.num_channels = 3;
			header.channels = (EXRChannelInfo *)malloc(sizeof(EXRChannelInfo) * header.num_channels);
			// Must be (A)BGR order, since most of EXR viewers expect this channel order.
			strncpy(header.channels[0].name, "B", 255); header.channels[0].name[strlen("B")] = '\0';
			strncpy(header.channels[1].name, "G", 255); header.channels[1].name[strlen("G")] = '\0';
			strncpy(header.channels[2].name, "R", 255); header.channels[2].name[strlen("R")] = '\0';

			header.pixel_types = (int *)malloc(sizeof(int) * header.num_channels);
			header.requested_pixel_types = (int *)malloc(sizeof(int) * header.num_channels);
			for (int i = 0; i < header.num_channels; i++) {
				header.pixel_types[i] = TINYEXR_PIXELTYPE_FLOAT; // pixel type of input image
				header.requested_pixel_types[i] = TINYEXR_PIXELTYPE_FLOAT; // pixel type of output image to be stored in .EXR
			}

			const char* err;
			int ret = SaveEXRImageToFile(&image, &header, filename, &err);
			if (ret != TINYEXR_SUCCESS) {
				std::cerr << "<<ERROR>> Unable to save to EXR file" << std::endl;
				std::cerr << "<<ERROR>> " << err << std::endl;
				return ret;
			}

			free(header.channels);
			free(header.pixel_types);
			free(header.requested_pixel_types);
			return true ;
		}
};

typedef t_EXR_IO<float> EXR_IO ;
