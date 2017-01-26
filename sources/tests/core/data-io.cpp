/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2015, 2016, 2017 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */


/* Test whether the 'load' and 'save' methods of 'data' work as expected.  We
 * do that by bitwise comparison of the files produced by 'save', which is
 * not ideal (we may want to compare the in-memory representations
 * instead.)  */

#include <core/args.h>
#include <core/data.h>
#include <core/data_storage.h>
#include <core/vertical_segment.h>
#include <core/plugins_manager.h>
#include <tests.h>

#include <string>
#include <iostream>
#include <sstream>

#include <cstring>
#include <cstdlib>

// Get the 'unlink' declaration.
#ifdef _WIN32
# include <io.h>
# define unlink _unlink
#else
# include <unistd.h>
#endif

using namespace alta;

// Try loading a simple example from a text-format stream.
static void test_simple_load_from_text()
{
    static const char example[] = "\
#DIM 1 3f\n\
#VS 0\n\
#PARAM_IN COS_TH\n\
#PARAM_OUT RGB_COLOR\n\
0 1 2 3\n\
4 5 6 7\n\
8 9 10 11\n";

    std::istringstream input(example);

    auto data = dynamic_pointer_cast<vertical_segment>(
        plugins_manager::load_data("vertical_segment", input));
    TEST_ASSERT(data != NULL);

    TEST_ASSERT(data->size() == 3);
    TEST_ASSERT(data->parametrization().input_parametrization()
                == params::COS_TH);
    TEST_ASSERT(data->parametrization().output_parametrization()
                == params::RGB_COLOR);

    // Currently we always get ASYMMETRICAL_CONFIDENCE_INTERVAL for backward
    // compatibility reasons.
    TEST_ASSERT(data->confidence_interval_kind()
                == vertical_segment::ASYMMETRICAL_CONFIDENCE_INTERVAL);

    TEST_ASSERT(data->get(0) == Eigen::Vector4d(0., 1., 2., 3.));
    TEST_ASSERT(data->get(1) == Eigen::Vector4d(4., 5., 6., 7.));
    TEST_ASSERT(data->get(2) == Eigen::Vector4d(8., 9., 10., 11.));

    // Check the result of the other 'get' methods.
    vec x, y_lower, y_upper;
    data->get(1, x, y_lower, y_upper);
    TEST_ASSERT(x == Eigen::VectorXd::Constant(1, 4.));
    TEST_ASSERT(y_lower == Eigen::Vector3d(5., 6., 7.));
    TEST_ASSERT(y_lower == y_upper);
    data->get(1, y_lower, y_upper);
    TEST_ASSERT(y_lower == Eigen::Vector3d(5., 6., 7.));
    TEST_ASSERT(y_lower == y_upper);

    // The "matrix view" includes confidence interval data.
    auto view = data->matrix_view();
    TEST_ASSERT(view.cols() == 1 + 3 + 2 * 3);
    TEST_ASSERT(view.rows() == 3);
    TEST_ASSERT(view.col(0) == Eigen::Vector3d(0., 4., 8.));
    TEST_ASSERT(view.col(1) == Eigen::Vector3d(1., 5., 9.));
    TEST_ASSERT(view.col(2) == Eigen::Vector3d(2., 6., 10.));
    TEST_ASSERT(view.col(3) == Eigen::Vector3d(3., 7., 11.));
}

// Files that are automatically deleted upon destruction.
class temporary_file
{
public:
    temporary_file()
    {
        static std::string suffix;
        _name = "t-data-io" + suffix;
        suffix += "-";
    }

    ~temporary_file()
    {
        ::unlink(_name.c_str());
    }

    const std::string& name() { return _name; };
    operator const std::string& () { return _name; }

private:
    std::string _name;
};

// Try hard to read N bytes from INPUT into BUF.
static std::streamsize read_exactly(std::istream &input, char *buf, std::streamsize n)
{
    std::streamsize total;

    for (total = 0;
         total < n && !input.eof();
         total += input.gcount())
    {
        input.read(&buf[total], n - total);
    }

    return total;
}

// Return true if the contents of FILE1 are identical to the contents of
// FILE2.
static bool files_are_equal(const std::string &file1, const std::string &file2)
{
    std::ifstream i1, i2;

    i1.exceptions(std::ios_base::failbit);
    i1.open(file1.c_str());
    i2.open(file2.c_str());
    i1.exceptions(std::ios_base::goodbit);

    while (!i1.eof() && !i2.eof())
    {
        char buf1[16384], buf2[16384];
        std::streamsize len1 = read_exactly(i1, buf1, sizeof buf1);
        std::streamsize len2 = read_exactly(i2, buf2, sizeof buf2);
        if (len1 != len2 || memcmp(buf1, buf2, len1) != 0) {
            return false;
        }
    }

    return i1.eof() && (i1.eof() == i2.eof());
}

int main(int argc, char** argv)
{
    std::string input_file;
    temporary_file temp_file1, temp_file2, temp_file3;

    if (argc > 1)
        // Process the user-specified file.
        input_file = argv[1];
    else
    {
        // Process the default test file.
        static const std::string data_file = "Kirby2.dat";
        std::string data_dir = getenv("TEST_DATA_DIRECTORY") != NULL
            ? getenv("TEST_DATA_DIRECTORY") : ".";
        input_file = data_dir + "/" + data_file;
    }

    // Simple test first.
    test_simple_load_from_text();

    // Try a sequence of loads and saves.
    try
    {
        std::ifstream input;
        input.open(input_file);

        // Use the standard load/save methods.
        auto sample1 = plugins_manager::load_data("vertical_segment", input);
        sample1->save(temp_file1);

        std::ifstream temp1;
        temp1.open(temp_file1);

        auto sample2 = plugins_manager::load_data("vertical_segment", temp1);
        sample2->save(temp_file2);

        // Now use the binary output format.
        std::ofstream out;
        out.open(temp_file3.name().c_str(), std::ofstream::binary);
        save_data_as_binary(out, *sample2);
        out.close();

        // This should automatically load using the binary format loader.
        std::ifstream temp3;
        temp3.open(temp_file3);
        auto sample3 = plugins_manager::load_data("vertical_segment", temp3);

        TEST_ASSERT(sample1->equals(*sample2));
        TEST_ASSERT(files_are_equal(temp_file1, temp_file2));
        TEST_ASSERT(sample2->equals(*sample3));

        TEST_ASSERT(sample1->min().size() == sample1->parametrization().dimX());
        TEST_ASSERT(sample1->max().size() == sample1->parametrization().dimX());

        TEST_ASSERT(sample1->min() == sample2->min());
        TEST_ASSERT(sample1->max() == sample2->max());
        TEST_ASSERT(sample1->min() == sample3->min());
        TEST_ASSERT(sample1->max() == sample3->max());

        // Make sure confidence interval data was preserved.
        auto vs_sample1 = dynamic_pointer_cast<vertical_segment>(sample1);
        auto vs_sample2 = dynamic_pointer_cast<vertical_segment>(sample2);
        auto vs_sample3 = dynamic_pointer_cast<vertical_segment>(sample3);
        TEST_ASSERT(vs_sample1 && vs_sample2 && vs_sample3);

        TEST_ASSERT(vs_sample1->confidence_interval_kind()
                    == vertical_segment::ASYMMETRICAL_CONFIDENCE_INTERVAL);
        TEST_ASSERT(vs_sample2->confidence_interval_kind()
                    == vertical_segment::ASYMMETRICAL_CONFIDENCE_INTERVAL);
        TEST_ASSERT(vs_sample3->confidence_interval_kind()
                    == vertical_segment::ASYMMETRICAL_CONFIDENCE_INTERVAL);

        // Compare the full "matrix views", which contains both X, Y, and the
        // confidence interval on Y.
        TEST_ASSERT(vs_sample1->matrix_view() == vs_sample2->matrix_view());
        TEST_ASSERT(vs_sample1->matrix_view() == vs_sample3->matrix_view());
    }
    CATCH_FILE_IO_ERROR(input_file);

    return EXIT_SUCCESS;
}
