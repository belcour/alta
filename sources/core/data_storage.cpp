/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2015, 2016, 2017 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include "data.h"
#include "data_storage.h"
#include "vertical_segment.h"

#include <iostream>
#include <limits>
#include <iomanip>
#include <cassert>

#ifdef __GLIBC__
# include <endian.h>
#endif

using namespace alta;
using namespace Eigen;

// A deleter for arrays, to work around the lack of array support in C++11's
// 'shared_ptr':
// <http://stackoverflow.com/questions/8947579/why-isnt-there-a-stdshared-ptrt-specialisation#8947700>.
static void delete_array(double *thing)
{
    delete[] thing;
}


static bool cosine_correction(vecref v, unsigned int dimX, unsigned int dimY,
                              alta::params::input in_param)
{
    using namespace alta;

    assert(v.size() == dimX + dimY);

    // Correction of the data by 1/cosine(theta_L)
    double factor = 1.0;
    double cart[6];
    params::convert(&v(0), in_param, params::CARTESIAN, cart);
    if(cart[5] > 0.0 && cart[2] > 0.0)
    {
        factor = 1.0/cart[5]*cart[2];
        for(unsigned int i = 0; i < dimY; ++i)
        {
            std::cout << "i = " << i << " dimx = " << dimX
                      << " + = " << i + dimX << " | " << v.size()
                      << std::endl;
            v(i + dimX) /= factor;
        }
        return true;
    }
    else return false;
}

// Return true if V is between MIN and MAX.
static bool within_bounds(const vecref v,
                          const vec& min, const vec& max)
{
    return (v.array() >= min.array()).all()
        && (v.array() <= max.array()).all();
}

// Return the 'ci_kind' value corresponding to VS_VALUE, an integer found in
// a '#VS' header.
static vertical_segment::ci_kind ci_kind_from_number(int vs_value)
{
    return vs_value == 2 ? vertical_segment::ASYMMETRICAL_CONFIDENCE_INTERVAL
        : (vs_value == 1 ? vertical_segment::SYMMETRICAL_CONFIDENCE_INTERVAL
           : vertical_segment::NO_CONFIDENCE_INTERVAL);
}

// Return the number used to represent KIND in '#VS' headers.
static int number_from_ci_kind(vertical_segment::ci_kind kind)
{
    return kind == vertical_segment::ASYMMETRICAL_CONFIDENCE_INTERVAL
        ? 2 : (kind == vertical_segment::SYMMETRICAL_CONFIDENCE_INTERVAL
               ? 1 : 0);
}

// Read a confidence interval on the output parameters from INPUT into V.
static void read_confidence_interval(std::istream& input,
                                     vecref v,
                                     vertical_segment::ci_kind kind,
                                     unsigned int dimX,
                                     unsigned int dimY,
                                     const alta::arguments& args)
{
    assert(v.size() == dimX + 3 * dimY);

    for(int i=0; i<dimY; ++i)
    {
        double min_dt = 0.0, max_dt = 0.0;

        if(i == 0 && kind == vertical_segment::ASYMMETRICAL_CONFIDENCE_INTERVAL)
        {
            input >> min_dt ;
            input >> max_dt ;
            min_dt = min_dt-v(dimX + i);
            max_dt = max_dt-v(dimX + i);
        }
        else if(i == 0 && kind == vertical_segment::SYMMETRICAL_CONFIDENCE_INTERVAL)
        {
            double dt ;
            input >> dt ;
            min_dt = -dt;
            max_dt =  dt;
        }
        else
        {
            // Confidence interval data not provided in INPUT.
            double dt = args.get_double("dt", 0.1d);
            min_dt = -dt;
            max_dt =  dt;
        }

        if(args.is_defined("dt-relative"))
        {
            v(dimX +   dimY+i) = v(dimX + i) * (1.0 + min_dt) ;
            v(dimX + 2*dimY+i) = v(dimX + i) * (1.0 + max_dt) ;
        }
        else if(args.is_defined("dt-max"))
        {
            v(dimX +   dimY+i) = v(dimX + i) + std::max(v(dimX + i) * min_dt, min_dt);
            v(dimX + 2*dimY+i) = v(dimX + i) + std::max(v(dimX + i) * max_dt, max_dt);
        }
        else
        {
            v(dimX +   dimY+i) = v(dimX + i) + min_dt ;
            v(dimX + 2*dimY+i) = v(dimX + i) + max_dt ;
        }

        // You can enforce the vertical segment to stay in the positive
        // region using the --data-positive command line argument. Note
        // that the data point is also clamped to zero if negative.
        if(args.is_defined("dt-positive"))
        {
            v(dimX +        i) = std::max(v(dimX +        i), 0.0);
            v(dimX +   dimY+i) = std::max(v(dimX +   dimY+i), 0.0);
            v(dimX + 2*dimY+i) = std::max(v(dimX + 2*dimY+i), 0.0);
        }

#ifdef DEBUG
        std::cout << "<<DEBUG>> vs = [" << v[dimX +   dimY+i] << ", " << v[dimX + 2*dimY+i] << "]" << std::endl;
#endif
    }
}

alta::data* alta::load_data_from_text(std::istream& input,
                                      const alta::arguments& header,
                                      const alta::arguments& args)
{
  vec min, max ;
  vec ymin, ymax;

  if(! header.is_defined("DIM")) {
    std::cerr << "<<ERROR>> Undefined dimensions ! ";
    std::cerr << "Please add DIM [int] [int] into the file header." << std::endl;
    throw;
  }

  params::input  in_param  = params::parse_input(header.get_string("PARAM_IN", "UNKNOWN_INPUT"));
  params::output out_param = params::parse_output(header.get_string("PARAM_OUT", "UNKNOWN_OUTPUT"));

  std::pair<int, int> dim = header.get_pair<int>("DIM");
  size_t row_count = dim.first + 3 * dim.second;

  min = args.get_vec("min", dim.first, -std::numeric_limits<float>::max()) ;
  max = args.get_vec("max", dim.first,  std::numeric_limits<float>::max()) ;
#ifdef DEBUG
  std::cout << "<<DEBUG>> data will remove outside of " << min << " -> " << max << " x-interval" << std::endl;
#endif

  ymin = args.get_vec("ymin", dim.second, -std::numeric_limits<float>::max()) ;
  ymax = args.get_vec("ymax", dim.second,  std::numeric_limits<float>::max()) ;
#ifdef DEBUG
  std::cout << "<<DEBUG>> data will remove outside of " << ymin << " -> " << ymax << " y-interval" << std::endl;
#endif

  auto kind = ci_kind_from_number(header.get_int("VS"));
  std::vector<double> content;

  // Now read the body.
  while(input.good())
  {
    std::string line ;
    std::getline(input, line) ;
    std::stringstream linestream(line) ;

    // Discard comments and empty lines.
    if(line.empty() || linestream.peek() == '#')
    {
      continue ;
    }
    else
    {
      auto start = content.size();

      // Read the data point x and y coordinates
      for(int i=0; i < dim.first + dim.second; ++i)
      {
          double item;
          linestream >> item;
          content.push_back(item);
      }

      // Save room for the confidence interval.
      for (int i = 0; i < 2 * dim.second; i++)
          content.push_back(0.);

      Map<VectorXd> v(&content[start], row_count);

      // Read the confidence interval data if available.
      read_confidence_interval(linestream, v, kind,
                               dim.first, dim.second, args);

      // Check if we need to filter out what we just read according to ARGS.
      // TODO: Move filtering to a post-parsing operation on 'data'.
      if (!(within_bounds(v.segment(0, dim.first), min, max)
            && within_bounds(v.segment(dim.first, dim.second),
                             ymin, ymax)))
      {
          content.resize(start);
      }
      else if (args.is_defined("data-correct-cosine"))
      {
          if (!cosine_correction(v.segment(0, dim.first + dim.second),
                                 dim.first, dim.second, in_param))
              content.resize(start);
      }
    }
  }

  std::cout << "<<INFO>> loaded input stream" << std::endl ;
  std::cout << "<<INFO>> loading data input of R^"
            << dim.first
            << " -> R^" << dim.second << std::endl ;
  std::cout << "<<INFO>> " << content.size() << " elements to fit" << std::endl ;


  double *raw_content = new double[content.size()];
  memcpy(raw_content, content.data(), content.size() * sizeof(double));

  parameters param(dim.first, dim.second, in_param, out_param);
  size_t element_count = content.size() / row_count;
  data* result = new vertical_segment(param, element_count,
                                      std::shared_ptr<double>(raw_content,
                                                                            delete_array));
  if(args.is_defined("data-correct-cosine"))
      result->save("/tmp/data-corrected.dat");

  return result;
}

void alta::save_data_as_text(std::ostream& out, const alta::data &data)
{
        using namespace alta;

    out << "#DIM " << data.parametrization().dimX() << " " << data.parametrization().dimY() << std::endl;
    out << "#PARAM_IN  "
        << params::get_name(data.parametrization().input_parametrization())
        << std::endl;
    out << "#PARAM_OUT "
        << params::get_name(data.parametrization().output_parametrization())
        << std::endl;

    for(int i=0; i < data.size(); ++i)
    {
        vec x = data.get(i);
        for(int j=0; j< data.parametrization().dimX() + data.parametrization().dimY(); ++j)
        {
                    out << std::setprecision(std::numeric_limits<double>::digits10)
                        << x[j] << "\t";
        }
        out << std::endl;
    }
}

void alta::save_data_as_binary(std::ostream &out, const alta::data& data)
{
    using namespace alta;

    auto maybe_vs = dynamic_cast<const vertical_segment*>(&data);
    auto kind = maybe_vs != NULL
        ? maybe_vs->confidence_interval_kind()
        : vertical_segment::NO_CONFIDENCE_INTERVAL;

    out << "#DIM " << data.parametrization().dimX() << " " << data.parametrization().dimY() << std::endl;
    out << "#PARAM_IN  "
        << params::get_name(data.parametrization().input_parametrization())
        << std::endl;
    out << "#PARAM_OUT "
        << params::get_name(data.parametrization().output_parametrization())
        << std::endl;
    out << "#FORMAT binary" << std::endl;
    out << "#VERSION 0" << std::endl;
    out << "#PRECISION ieee754-double" << std::endl;
    out << "#SAMPLE_COUNT " << data.size() << std::endl;
    out << "#VS " << number_from_ci_kind(kind) << std::endl;

    // FIXME: Note: on non-glibc systems, both macros may be undefined, so
    // the conditional is equivalent to "#if 0 == 0", which is usually what
    // we want.
#if __BYTE_ORDER == __LITTLE_ENDIAN
    out << "#ENDIAN little" << std::endl;
#else
    out << "#ENDIAN big" << std::endl;
#endif

    out << "#BEGIN_STREAM" << std::endl;

    if (kind == vertical_segment::NO_CONFIDENCE_INTERVAL)
    {
        // No confidence interval to be saved.
        for(int i=0; i < data.size(); ++i)
        {
            vec sample = data.get(i);
            const double *numbers = sample.data();

            assert(sample.size() == data.parametrization().dimX() + data.parametrization().dimY());
            out.write((const char *)numbers, sample.size() * sizeof(*numbers));
        }
    }
    else
    {
        // Saving data along with the confidence interval.
        auto matrix = maybe_vs->matrix_view();
        auto byte_count = matrix.cols() * matrix.rows() * sizeof(double);

        // MATRIX has no stride so we can access the underlying storage
        // directly.
        out.write((const char *)matrix.data(), byte_count);
    }

    out << std::endl << "#END_STREAM" << std::endl;
}

alta::data* alta::load_data_from_binary(std::istream& in, const alta::arguments& header)
{
    using namespace alta;

    // FIXME: For now we make a number of assumptions.
    assert(header["FORMAT"] == "binary");
    assert(header.get_int("VERSION") == 0);
    assert(header["PRECISION"] == "ieee754-double");
#if __BYTE_ORDER == __LITTLE_ENDIAN
    assert(header["ENDIAN"] == "little");
#else
    assert(header["ENDIAN"] == "big");
#endif

    std::pair<int, int> dim = header.get_pair<int>("DIM");
    assert(dim.first > 0 && dim.second > 0);

    auto kind = ci_kind_from_number(header.get_int("VS"));

    in.exceptions(std::ios_base::failbit);

    int sample_count = header.get_int("SAMPLE_COUNT");
      if(sample_count <= 0) {
         std::cerr << "<<ERROR>> Uncorrect or not samples count in the header, please check \'SAMPLE_COUNT\'" << std::endl;
      }

      size_t ci_rows =
          kind == vertical_segment::ASYMMETRICAL_CONFIDENCE_INTERVAL
          ? 2 : (kind == vertical_segment::SYMMETRICAL_CONFIDENCE_INTERVAL
                 ? 1 : 0);

      size_t rows = dim.first + dim.second + ci_rows * dim.second;
      size_t element_count = sample_count * rows;
      double *content = new double[element_count];
      size_t byte_count = element_count * sizeof *content;

      // TODO: Arrange to use mmap and make it zero-alloc and zero-copy.
      for (std::streamsize total = 0;
           total < byte_count && !in.eof();
           total += in.gcount())
      {
          in.read((char *) content + total, byte_count - total);
      }

    parameters param(dim.first, dim.second,
                     params::parse_input(header["PARAM_IN"]),
                     params::parse_output(header["PARAM_OUT"]));

    return new alta::vertical_segment(param, sample_count,
                                      std::shared_ptr<double>(content,
                                                              delete_array),
                                      kind);
}
