/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2015 Université de Montréal
   Copyright (C) 2018 Unity

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#pragma once

// Pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

// ALTA include
#include <core/common.h>

// STL include
#include <iostream>

/* Specific function call to set a vector's element
 */
static void vec_set_item(vec& x, int i, double a) {
	x[i] = a;
}

/* Operators on vec
 */
static vec vec_add(const vec& a, const vec& b) {
   return a + b;
}
static vec vec_sub(const vec& a, const vec& b) {
   return a - b;
}

/* Specific convert a vec to a string
 */
static std::string vec_str(const vec& x) {
   std::stringstream stream;
   stream << "[";
   for(int i=0; i<x.size(); ++i) {
      if(i > 0) { stream << ", "; }
      stream << x[i];
   }
   stream << "]";
   return stream.str();
}

/* Register the vector class to python
 */
inline void register_vec(py::module& m) {

	py::class_<vec>(m, "vec")
        .def(py::init<vec>())
        .def(py::init([]( const std::vector<double>& a) {
            vec v(a.size());
            memcpy((void*)v.data(), (const void*)a.data(), a.size()*sizeof(double));
            return v;
        }))
		.def("__add__", &vec_add)
		.def("__sub__", &vec_sub)
		// .def("__len__", &vec::size)
        .def("__getitem__", [](const vec &s, unsigned int i) {
            if (i >= s.size()) throw py::index_error();
            return s[i];
        })
		.def("__setitem__", &vec_set_item)
		.def("__str__", &vec_str);
}
