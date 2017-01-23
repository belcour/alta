#pragma once
/* EXR_IO provides static method to either load an RGB EXR file into a floatting
 * point array or save a RGB floatting point array to disk._data
 *
 * Original author: Cyril Soler
 * Modified by: Laurent Belcour
 */

#include <stdexcept>
#include <cassert>

#ifdef USE_OPENEXR
	#include <ImfRgbaFile.h>
	#include <ImfArray.h>

	// XXX: This header is not installed as of version 2.2.0 of OpenEXR, see
	// <https://lists.nongnu.org/archive/html/openexr-devel/2016-06/msg00001.html>
	// and <https://github.com/openexr/openexr/pull/184>.
	#include <ImfStdIO.h>
#else
	#define TINYEXR_IMPLEMENTATION
	#include <tinyexr/tinyexr.h>
#endif

template<typename FType>
class t_EXR_IO
{
	public:
		static bool LoadEXR(std::istream& input, int& W,int& H, FType *& pix,int nC=3)
		{
#ifdef USE_OPENEXR
			/* XXX: OpenEXR implements its own IStream and OStream classes, but
			 * they are completely independent from those of libstdc++.  The
			 * closest thing it has is 'StdIFStream', hence this hack.
			 */
			std::ifstream* ifstream = dynamic_cast<std::ifstream*>(&input);
			assert(ifstream != NULL);

			Imf::StdIFStream iifstream(*ifstream, "unknown file name");

			Imf::RgbaInputFile file(iifstream);
			Imath::Box2i dw = file.dataWindow();

			W = dw.max.x - dw.min.x + 1;
			H = dw.max.y - dw.min.y + 1;

			Imf::Array2D<Imf::Rgba> pixels;
			pixels.resizeErase(H, W);

			file.setFrameBuffer (&pixels[0][0] - dw.min.x - dw.min.y * W, 1, W);
			file.readPixels (dw.min.y, dw.max.y);

			pix = new FType[W*H*3] ;

			switch(nC)
			{
				case 3: for(int i=0;i<H;++i)
							  for(int j=0;j<W;++j)
							  {
								  pix[3*(j+i*W)+0] = pixels[H-i-1][j].r ;
								  pix[3*(j+i*W)+1] = pixels[H-i-1][j].g ;
								  pix[3*(j+i*W)+2] = pixels[H-i-1][j].b ;
							  }
						  break ;

				case 1: for(int i=0;i<H;++i)
							  for(int j=0;j<W;++j)
								  pix[j+i*W] = 0.3*pixels[H-i-1][j].r + 0.59*pixels[H-i-1][j].g + 0.11*pixels[H-i-1][j].b ;
						  break ;
				default:
						  throw std::runtime_error("Unexpected case in EXR_IO.") ;
			}
#else
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
			_r = ParseEXRVersionFromMemory(&_version, &_memory[0]);
			if(_r != TINYEXR_SUCCESS) {
				std::cerr << "<<ERROR>> Could not load EXR version from stream" << std::endl;
				return false;
			} else {
				std::cout << "<<DEBUG>> Loading a v" << _version.version << " EXR file" << std::endl;
			}

			/* Load the EXR header */
			EXRHeader _header;
  			InitEXRHeader(&_header);
			_r = ParseEXRHeaderFromMemory(&_header, &_version, &_memory[0], &_err);
			if(_r != TINYEXR_SUCCESS) {
				std::cerr << "<<ERROR>> Could not load EXR header from stream" << std::endl;
				std::cerr << "<<ERROR>> " << _err << std::endl;
				return false;
			} else {
				std::cout << "<<DEBUG>> Loading a " << _header.num_channels << " channels EXR file" << std::endl;
			}

			/* Set the requested pixel types  */
			for(int i=0; i<_header.num_channels; ++i) {
				_header.pixel_types[i] = TINYEXR_PIXELTYPE_FLOAT;
			}

			/* Load the EXR image */
			EXRImage _image;
  			InitEXRImage(&_image);

			_r = LoadEXRImageFromMemory(&_image, &_header, &_memory[0], &_err);
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
			for(int k=0; k<_header.num_channels; ++k) {
				float* img = reinterpret_cast<float *>(_image.images[k]);

				for(int i=0; i<W*H; ++i) {
					const double val = (FType)(img[i]);
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

			/*! \todo Free TinyEXR memory : _image, _header and _version */
#endif
			return true ;
		}

		static bool SaveEXR(const char *filename,int W,int H, const FType *pix,int nC=3)
		{
#ifdef USE_OPENEXR
			Imf::Array2D<Imf::Rgba> pixels(H,W);

			/* Convert separated channel representation to per pixel representation */
			switch(nC)
			{
				case 3:
					for (int row=0;row<H;row++)
						for(int i=0;i<W;i++)
						{
							Imf::Rgba &p = pixels[H-row-1][i];

							p.r = pix[3*(i+row*W)+0] ;
							p.g = pix[3*(i+row*W)+1] ;
							p.b = pix[3*(i+row*W)+2] ;
							p.a = 1.0 ;
						}
					break ;

				case 1:
					for (int row=0;row<H;row++)
						for(int i=0;i<W;i++)
						{
							Imf::Rgba &p = pixels[H-row-1][i];

							p.r = pix[i+row*W] ;
							p.g = pix[i+row*W] ;
							p.b = pix[i+row*W] ;
							p.a = FType(1.0) ;
						}
					break ;
				default:
					throw std::runtime_error("Unexpected case in EXR_IO.") ;
			}

			Imf::RgbaOutputFile file(filename, W, H, Imf::WRITE_RGBA);
			file.setFrameBuffer(&pixels[0][0], 1, W);
			file.writePixels(H);

#else
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
				images[0][i] = pix[3*i+0];
				images[1][i] = pix[3*i+1];
				images[2][i] = pix[3*i+2];
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
				header.requested_pixel_types[i] = TINYEXR_PIXELTYPE_HALF; // pixel type of output image to be stored in .EXR
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
#endif
			return true ;
		}
};

typedef t_EXR_IO<float> EXR_IO ;