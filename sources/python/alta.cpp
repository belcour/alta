// Boost includes
#include <boost/python.hpp>

// ALTA include
#include <core/common.h>
#include <core/ptr.h>
#include <core/function.h>
#include <core/plugins_manager.h>

// STL include
#include <iostream>

#define bp boost::python


/* The following code register ALTA's shared pointer as a valid shared ptr
 * to be used by boost::python .
 */
template <typename T> 
T* get_pointer(ptr<T> const& p) {
  return const_cast<T*>(p.get());
}

namespace boost {
	namespace python {
   		template <typename T>
    	struct pointee< ::ptr<T> > {
        	typedef T type;
    	};
	}
}


/* Wrapper to ALTA's vec class. This is only here to allow init with Python's
 * list.
 *
 * TODO: Make sure that the value passed to this vector are floatting point
 *       convertible.
 */
struct my_vec : public vec {
    my_vec() : vec() {}

    my_vec(const bp::list& l) : vec(bp::len(l)) {
        for(size_t i=0; i<bp::len(l); ++i) {
            (*this)[i] = bp::extract<double>(l[i]);
        }
    }
};


/* This class is a wrapper to ALTA's arguments class to add Python specific
 * behaviour such as dictionnary initialization.
 */
struct python_arguments : public arguments {
	python_arguments() : arguments() {}
	python_arguments(bp::dict d) : arguments() {
		bp::list keys = d.keys();
		for(int i=0; i<bp::len(keys); ++i) {
			const std::string s_key = bp::extract<std::string>(keys[i]);
			const std::string s_val = bp::extract<std::string>(d[keys[i]]);
			this->update(s_key, s_val);
		}
	}
};


/* Create a data object from a plugin's name and the data filename. This 
 * function is here to accelerate the loading of data file.
 */
ptr<data> load_data(const std::string& plugin_name, const std::string& filename) {
	ptr<data> d = plugins_manager::get_data(plugin_name);
	d->load(filename);
	return d;
}


/* Creating functions for the plugins_manager calls
 * 
 * TODO: Throw python exceptions if the function is not correctly created.
 *       Those function should disapear when the return type of get_Function
 *       in the plugin_manager will be ptr<function>.
 */
ptr<function> get_function(const std::string& filename) {
    ptr<function> func(plugins_manager::get_function(filename));
    if(!func) {
    	std::cerr << "<<ERROR>> no function created" << std::endl;
    }
	return func;
}
ptr<function> get_function_from_args(const arguments& args) {
    ptr<function> func(plugins_manager::get_function(args));
    if(!func) {
    	std::cerr << "<<ERROR>> no function created" << std::endl;
    }
    return func;
}

/* Exporting the ALTA module */
BOOST_PYTHON_MODULE(alta)
{
	// Argument class
	//
	bp::class_<python_arguments>("arguments")
		.def(bp::init<>())
		.def(bp::init<bp::dict>())
		.def("update", &arguments::update);

	bp::class_<my_vec>("vec")
		.def(bp::init<bp::list>())
		/*        .def(" __setitem__", &my_vec::operator[])*/;

	// Function interface
	//
	bp::class_<function, ptr<function>, boost::noncopyable>("function", bp::no_init)
		.def("value", &function::value)
		.def("load",  &function::load)
		.def("save",  &function::save);
	bp::def("get_function", get_function);
	bp::def("get_function", get_function_from_args);

	// Data interface
	//
	bp::class_<data, ptr<data>, boost::noncopyable>("data", bp::no_init)
		.def("load", static_cast< void(data::*)(const std::string&)>(&data::load))
		.def("size", &data::size)
		.def("save", &data::save);
	bp::def("get_data",  plugins_manager::get_data);
	bp::def("load_data", load_data);

	// Fitter interface
	//
	bp::class_<fitter, ptr<fitter>, boost::noncopyable>("fitter", bp::no_init)
		.def("fit_data", &fitter::fit_data);
	bp::def("get_fitter", plugins_manager::get_fitter);
}
