/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2018 Unity

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#pragma once

// Pybind11 includes
#include <pybind11/pybind11.h>
namespace py = pybind11;

// ALTA include
#include <core/args.h>
#include <core/ptr.h>
#include <core/data.h>

using namespace alta;


/* This class is a wrapper to ALTA's arguments class to add Python specific
 * behaviour such as dictionnary initialization.
 *
 * Right now it does not handle automatic conversion in function call for
 * example. The following code is not possible:
 *
 *    import alta
 *    alta.get_data('data_merl', {'params' : 'STARK_2D'})
 *
 * Instead, one has to construct the arguments object from the ALTA library
 * to use it afterwards:
 *
 *    import alta
 *    args = alta.arguments({'params' : 'STARK_2D'})
 *    alta.get_data('data_merl', args)
 */
struct python_arguments : public arguments {
	python_arguments() : arguments() {}
	python_arguments(py::dict d) : arguments() {
		for(auto item : d) {
			const std::string s_key = std::string(py::str(item.first));
			const std::string s_val = std::string(py::str(item.second));
			this->update(s_key, s_val);
		}
	}
};


/* Register the `arguments` class to python */
inline void register_args(py::module& m) {
    py::class_<python_arguments>(m, "arguments")
        .def(py::init<>())
		.def(py::init<py::dict>())
        .def("__getitem__", &arguments::operator[])
		.def("update", &arguments::update);
}