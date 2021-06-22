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
#include <core/common.h>
#include <core/ptr.h>
#include <core/data.h>

using namespace alta;

/* Create a data object from a plugin's name and the data filename. This
 * function is here to accelerate the loading of data file.
 */
static ptr<data> load_data(const std::string& plugin_name,
                           const std::string& filename)
{
    return plugins_manager::load_data(filename, plugin_name);
}

static ptr<data> get_data_with_args(const std::string& plugin_name,
                                    size_t size,
                                    const parameters& params,
                                    const arguments& args)
{
  return plugins_manager::get_data(plugin_name, size, params, args);
}

static ptr<data> get_data(const std::string& plugin_name, size_t size,
                          const parameters& params)
{
  return plugins_manager::get_data(plugin_name, size, params);
}

/* Register the `data` class to python
 */
inline void register_data(py::module& m) {
   py::class_<data, ptr<data>>(m, "data")
      .def("__getitem__", &data::get)
      .def("__setitem__", &data::set)
      .def("__len__", &data::size)
      .def("size",  &data::size)
      .def("get",   &data::get)
      .def("set",   &data::set)
      .def("value", &data::value)
      .def("save",  &data::save)
      .def("parametrization", &data::parametrization);
   m.def("get_data",  get_data);
   m.def("get_data",  get_data_with_args);
   m.def("load_data", load_data);
}
