/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2014, 2015, 2016 Inria
   Copyright (C) 2015 Université de Montréal
   Copyright (C) 2018 Unity

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

// Pybind11 includes
#include <pybind11/pybind11.h>
namespace py = pybind11;

// ALTA include
#include <core/common.h>
#include <core/ptr.h>
#include <core/function.h>
#include <core/rational_function.h>
#include <core/plugins_manager.h>
#include <core/vertical_segment.h>
#include <core/metrics.h>

// STL include
#include <iostream>

// Local includes
#include "wrapper_args.h"
#include "wrapper_vec.h"
#include "wrapper_data.h"

using namespace alta;



/* Creating functions for the plugins_manager calls
 *
 * TODO: Throw python exceptions if the function is not correctly created.
 *       Those function should disapear when the return type of get_Function
 *       in the plugin_manager will be ptr<function>.
 */
static ptr<function> get_function(const std::string& plugin_name,
                                  const parameters& params)
{
    ptr<function> func(plugins_manager::get_function(plugin_name,
                                                     params));
    if(!func) {
    	std::cerr << "<<ERROR>> no function created" << std::endl;
    }
    return func;
}
static ptr<function> get_function_from_args(const python_arguments& args,
                                            const parameters& params) {
    ptr<function> func(plugins_manager::get_function(args, params));
    if(!func) {
    	std::cerr << "<<ERROR>> no function created" << std::endl;
    }
    return func;
}

/* Load a function object from a file. The arguments object is never used here.
 * The file is supposed to load the function correctly.
 */
static ptr<function> load_function(const std::string& filename) {
    return ptr<function>(plugins_manager::load_function(filename));
}
static ptr<function> load_function_with_args(const std::string& filename, const python_arguments&) {
    return ptr<function>(plugins_manager::load_function(filename));
}

/* Loading a function object from file
 *
 * TODO: Add exceptions
 */
static void load_from_file_with_args(const ptr<function>& func,
                                     const std::string& filename,
                                     const python_arguments& args) {

   // Open a stream
   std::ifstream file;
   file.open(filename.c_str()) ;
   if(!file.is_open()) {
      std::cerr << "<<ERROR>> unable to open file \"" << filename << "\"" << std::endl ;
      return;
   }

   // Parse the associated header
   arguments header = arguments::parse_header(file);

   // Load the function stream
   func->load(file);
}
static void load_from_file(const ptr<function>& func, const std::string& filename) {

   python_arguments args;
   load_from_file_with_args(func, filename, args);
}

/* Operators on function object. This provide the ability to create compounds
 * and product in the command line. This is only possible for nonlinear_functions
 *
 * TODO: The compound and product function should store the shared pointer to the
 * function objects. They might stay in memory longer than the input functions.
 */
static ptr<function> add_function(const ptr<function>& f1, const ptr<function>& f2) {
	// Convert to a nonlinear function
	ptr<nonlinear_function> nf1 = dynamic_pointer_cast<nonlinear_function>(f1);
	ptr<nonlinear_function> nf2 = dynamic_pointer_cast<nonlinear_function>(f2);

	if(nf1 && nf2) {
    std::vector<ptr<nonlinear_function> > functions;
    functions.push_back(nf1);
    functions.push_back(nf2);
    std::vector<arguments> args;
    args.push_back(alta::arguments());
    args.push_back(alta::arguments());
    auto cf = new compound_function(functions, args);
    return ptr<function>(cf);

	// Failure case, one of the function is a NULL ptr.
	} else {
		std::cerr << "<<ERROR>> One of the input functions is NULL" << std::endl;
		return ptr<function>(NULL);
	}
}
static ptr<function> mult_function(const ptr<function>& f1, const ptr<function>& f2) {
	// Convert to a nonlinear function
	ptr<nonlinear_function> nf1 = dynamic_pointer_cast<nonlinear_function>(f1);
	ptr<nonlinear_function> nf2 = dynamic_pointer_cast<nonlinear_function>(f2);

	// Check
	if(!nf1 || !nf2) {
		std::cerr << "<<ERROR>> One of the input function of the product is NULL" << std::endl;
		return ptr<function>(NULL);
	}

	product_function* pf = new product_function(nf1, nf2);
	return ptr<function>(pf);
}

/* Setting/Get the parameters of a function object
 *
 * TODO: Add the rational function interface
 */
static void set_function_params(ptr<function>& f, const vec& v) {

	// Try to set the parameter as a nonlinear function
	ptr<nonlinear_function> nf = dynamic_pointer_cast<nonlinear_function>(f);
	if(nf) {
		if(nf->nbParameters() == v.size()) {
			nf->setParameters(v);
		} else {
			std::cerr << "<<ERROR>> Vector of params and function have different sizes." << std::endl;
		}
		return;
	}

	// ptr<rational_function> rf = dynamic_pointer_cast<rational_function>(f);
	// if(rf) {
	// 	int np, nq;
	// 	rf->size(np, nq);

	// 	if(np+nq == v.size()) {
	// 		rf->setParameters(v);
	// 	} else {
	// 		std::cerr << "<<ERROR>> Vector of parameters has different size that the functions number of parameters" << std::endl;
	// 	}
	// 	return;
	// }
}

static vec get_function_params(ptr<function>& f) {
		// Try to set the parameter as a nonlinear function
	ptr<nonlinear_function> nf = dynamic_pointer_cast<nonlinear_function>(f);
	if(nf) {
		return nf->parameters();
	}

	std::cerr << "<<ERROR>> Parameters cannot be retrieved" << std::endl;
	vec res(1);
	return res;
}

/* Save a function object to a file, without any argument option. This will save the function
 * object in ALTA's format.
 */
static void save_function_without_args(const ptr<function>& f, const std::string& filename) {
   arguments args;
   f->save(filename, args);
}

/* Fitter interface to allow to launch fit with and without providing an
 * arguments object.
 */
static bool fit_data_without_args(const ptr<fitter>& _fitter,
                                  const ptr<data>& _data,
                                  ptr<function>& _func) {
   arguments args;
   return _fitter->fit_data(_data, _func, args);
}

static bool fit_data_with_args(ptr<fitter>& _fitter, const ptr<data>& _data,
                               ptr<function>& _func,
                               const python_arguments& args) {
   _fitter->set_parameters(args);
   return _fitter->fit_data(_data, _func, args);
}


/* Softs functions. Those function recopy the softs main function, without
 * the command line arguments.
 * TODO: Add the command line arguments in the parameters
 */
static ptr<data> data2data(const data* d_in, const parameters& target)
{
    auto delete_array = [](double* array)
    {
        delete[] array;
    };

    // For each sample we store the X and Y values, followed by the
    // confidence interval on Y.
    size_t sample_size = target.dimX() + 3 * target.dimY();

    double* content = new double[d_in->size() * sample_size];

    for (auto sample = 0; sample < d_in->size(); sample++)
    {
        // Copy the input vector
        vec x = d_in->get(sample);
        params::convert(&x[0],
                        d_in->parametrization().input_parametrization(),
                        target.input_parametrization(),
                        &content[sample * sample_size]);
        params::convert(&x[d_in->parametrization().dimX()],
                        d_in->parametrization().output_parametrization(),
                        d_in->parametrization().dimY(),
                        target.output_parametrization(),
                        target.dimY(),
                        &content[sample * sample_size + target.dimX()]);

        // Set the confidence interval on Y to zero.
        for (auto i = 0; i < 2 * target.dimY(); i++)
        {
            content[sample * sample_size + target.dimX() + target.dimY()
                    + i] = 0.;
        }
    }

    data* result = new vertical_segment(target, d_in->size(),
                                        std::shared_ptr<double>(content,
                                                                delete_array));
    return ptr<data>(result);
}

/* This function provides a similar behaviour that the brdf2data function.
 * An input function object is evaluated on the input position of a data
 * object. To do so, the data object must contains some positions. It can
 * be so when the data object is a laoded data sample or when the data type
 * has predefined sample sets.
 */
static void brdf2data(const ptr<function>& f, ptr<data>& d) {
	if(d->size() == 0) {
		std::cerr << "<<ERROR>> Please provide a data object with a sample structure or load a data file with defined positions." << std::endl;
		return;
	}

	vec temp(f->parametrization().dimX());
	for(int i=0; i<d->size(); ++i) {
		// Copy the input vector
		vec x = d->get(i);

		// Convert the data to the function's input space.
      if(f->parametrization().input_parametrization() == params::UNKNOWN_INPUT) {
	    	memcpy(&temp[0], &x[0], f->parametrization().dimX()*sizeof(double));
	    } else {
			params::convert(&x[0],
                      d->parametrization().input_parametrization(),
                      f->parametrization().input_parametrization(),
                      &temp[0]);
		}
		vec y = f->value(temp);

		for(int j=0; j<d->parametrization().dimY(); ++j) {
			x[d->parametrization().dimX() + j] = y[j];
		}

		d->set(i, x);
	}
}

/* Compute distance metric between 'in' and 'ref'.
 */
static py::dict data2stats(const ptr<data>& in, const ptr<data>& ref) {
   // Compute the metrics
   errors::metrics res;
   errors::compute(in.get(), ref.get(), nullptr, res);

   // Fill the resulting Python vector
   py::dict py_res;
   for(auto rpair : res) {
      py_res[py::str(rpair.first)] = rpair.second;
   }
   return py_res;
}

inline void register_function(py::module& m) {
	py::class_<function, ptr<function>>(m, "function")
		.def("__add__",  &add_function)
		.def("__mul__",  &mult_function)
		.def("__rmul__", &mult_function)
		.def("value",    &function::value)
		.def("load",     &load_from_file)
		.def("load",     &load_from_file_with_args)
		.def("save",     &function::save)
        .def("save",     &save_function_without_args)
		.def("set",      &set_function_params)
		.def("get",      &get_function_params);
	m.def("get_function", get_function, py::arg("name") = "nonlinear_function_diffuse",
                                            py::arg("param") = parameters(6, 3, params::CARTESIAN, params::RGB_COLOR));
	m.def("get_function", get_function_from_args);
	m.def("load_function", load_function);
	m.def("load_function", load_function_with_args);
}

inline void register_fitter(py::module& m) {
	py::class_<fitter, ptr<fitter>>(m, "fitter")
		.def("fit_data", &fit_data_with_args)
                .def("fit_data", &fit_data_without_args);
	m.def("get_fitter",  plugins_manager::get_fitter);
}

#define STRINGIFY_(x) #x
#define STRINGIFY(x)  STRINGIFY_(x)

PYBIND11_MODULE(alta, m) {
    m.doc() = "ALTA python bindinds";
    m.def("norm", &norm, "Compute the norm of vector 'v'", py::arg("v"));

    /* Register the 'parameters' class and the enums */
    py::class_<parameters>(m, "parameters")
      .def(py::init<unsigned int, unsigned int, params::input, params::output>());

#define PARAM_VALUE(name)                       \
   .value(STRINGIFY(name), params:: name)

   py::enum_<params::input>(m, "input_parametrization")
        PARAM_VALUE(RUSIN_TH_PH_TD_PD)
        PARAM_VALUE(RUSIN_TH_PH_TD)
        PARAM_VALUE(RUSIN_TH_TD_PD)
        PARAM_VALUE(RUSIN_TH_TD)
        PARAM_VALUE(RUSIN_VH_VD)
        PARAM_VALUE(RUSIN_VH)
        PARAM_VALUE(COS_TH_TD)
        PARAM_VALUE(COS_TH)
        PARAM_VALUE(SCHLICK_TK_PK)
        PARAM_VALUE(SCHLICK_VK)
        PARAM_VALUE(SCHLICK_TL_TK_PROJ_DPHI)
        PARAM_VALUE(COS_TK)
        PARAM_VALUE(RETRO_TL_TVL_PROJ_DPHI)
        PARAM_VALUE(STEREOGRAPHIC)
        PARAM_VALUE(SPHERICAL_TL_PL_TV_PV)
        PARAM_VALUE(COS_TLV)
        PARAM_VALUE(COS_TLR)
        PARAM_VALUE(ISOTROPIC_TV_TL)
        PARAM_VALUE(ISOTROPIC_TV_TL_DPHI)
        PARAM_VALUE(ISOTROPIC_TV_PROJ_DPHI)
        PARAM_VALUE(ISOTROPIC_TL_TV_PROJ_DPHI)
        PARAM_VALUE(ISOTROPIC_TD_PD)
        PARAM_VALUE(STARK_2D)
        PARAM_VALUE(STARK_3D)
        PARAM_VALUE(NEUMANN_2D)
        PARAM_VALUE(NEUMANN_3D)
        PARAM_VALUE(CARTESIAN)
        PARAM_VALUE(UNKNOWN_INPUT);

    py::enum_<params::output>(m, "output_parametrization")
        PARAM_VALUE(INV_STERADIAN)
        PARAM_VALUE(INV_STERADIAN_COSINE_FACTOR)
        PARAM_VALUE(ENERGY)
        PARAM_VALUE(RGB_COLOR)
        PARAM_VALUE(XYZ_COLOR)
        PARAM_VALUE(UNKNOWN_OUTPUT);

#undef PARAM_VALUE

    /* Register the `args` and `vec` types */
    register_args(m);
    register_vec(m);

    /* register the `data`, `function` and `fitter` objects */
    register_data(m);
    register_function(m);
    register_fitter(m);

	/* register `softs` */
	m.def("data2data",  data2data);
	m.def("data2stats", data2stats);
	m.def("brdf2data",  brdf2data);
}
