/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2015 CNRS
   Copyright (C) 2013, 2014, 2016 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include "function.h"		

#include "common.h"
#include "params.h"
#include "plugins_manager.h"

#include <cassert>

using namespace alta;

/*--- Functions implementation ----*/

void function::bootstrap(const ptr<data>, const arguments& args)
{
  #ifdef BOOTSTRAP_DEBUG
  std::cout << __FILE__ << " " << __LINE__ << " args = " << args << std::endl;
  #endif

  // If the bootstrap option contains a filename, load it
  if(args.is_defined("bootstrap") && !args["bootstrap"].empty())
  {
    std::ifstream in(args["bootstrap"].c_str());
    if(in.is_open())
    {
      if(!load(in))
			{
				std::cerr << "<<ERROR>> cannot bootstrap from file \"" << args["bootstrap"] << "\"" << std::endl;
			}
    }
    else // No file for the boostrap arguments
    {
      std::cerr << " Bootstrapping from command line is not implemented" << std::endl;
			NOT_IMPLEMENTED();
    }
  }//end of if bootstrap is defined and not empty
    
}

void function::save(const std::string& filename, const arguments& args) const
{
	bool const is_alta     = !args.is_defined("export") || args["export"] == "alta";
	bool const is_cpp      = args["export"] == "C++" || args["export"] == "c++" || args["export"] == "cpp";
	bool const is_explorer = args["export"] == "explorer";
	bool const is_shader   = args["export"] == "shader" || is_explorer;
	bool const is_matlab   = args["export"] == "matlab";

	// Open the file
	std::ofstream file(filename.c_str());
	if(!file.is_open())
	{
		std::cerr << "<<ERROR>> unable to open output file for writing" << std::endl;
	}

	// If the export is the alta format, use the maximum precision formatting
	if(is_alta)
	{
		file.precision(10);
	}

	if(is_explorer)
	{
		file << "analytic" << std::endl;
		file << std::endl;
		file << "::begin shader" << std::endl;
	}

	// Save common header
	save_header(file, args);

	// Save function definition
	save_body(file, args);

	if(is_cpp)
	{
		file << "#include <Eigen/Core>" << std::endl;
		file << "typedef Eigen::VectorXd vec;" << std::endl;
		file << std::endl;
		file << "vec brdf(const vec& x)" << std::endl;
		file << "{" << std::endl;
		file << "\tvec res(" << parametrization().dimY() << ");" << std::endl;
	}
	else if(is_matlab)
	{
		file << "function res = brdf(x)" << std::endl;
		file << "{" << std::endl;
		file << "\tres = zeros(" << parametrization().dimY() << ");" << std::endl;
		file << "\t";
	}
	else if(is_shader)
	{
		file << "vec3 BRDF(vec3 L, vec3 V, vec3 N, vec3 X, vec3 Y)" << std::endl;
		file << "{" << std::endl;
		file << "\tvec3 res = ";
	}

	// Save fit data
	save_call(file, args);
	if(is_cpp || is_shader)
	{
		file << "\treturn res;" << std::endl;
		file << "}" << std::endl;
	}
	else if(is_matlab)
	{
		file << ";" << std::endl;
		file << "\treturn res;" << std::endl;
		file << "endfunction" << std::endl;
	}
	file << std::endl;

	if(is_explorer)
	{
		file << "::end shader" << std::endl;
	}
}
		
//! \brief save the header of the output function file. The header should
//! store general information about the fit such as the command line used
//! the dimension of the fit. L2 and L_inf distance could be added here.
void function::save_header(std::ostream& out, const arguments& args) const
{
	if(!args.is_defined("export"))
	{
		out << "#ALTA FUNC HEADER" << std::endl;
		out << "#CMD " << args.get_cmd() << std::endl;
		out << "#DIM " << parametrization().dimX() << " " << parametrization().dimY() << std::endl;
		out << "#PARAM_IN  " << params::get_name(parametrization().input_parametrization()) << std::endl;
		//out << "#PARAM_OUT " << params::get_name(output_parametrization()) << std::endl;*
		if(args.is_defined("export-append")) 
		{
			out << args["export-append"] << std::endl;
		}
		out << "#ALTA HEADER END" << std::endl;
		out << std::endl;
	}
}

//! \brief save function specific data. This has no use for ALTA export
//! but allows to factorize the code in the C++ or matlab export by
//! defining function calls that are common to all the plugins.
void function::save_body(std::ostream&, const arguments&) const
{


}

//! \brief save object specific information. For an ALTA export the
//! coefficients will be exported. For a C++ or matlab export, the call
//! to the associated function will be done.
void function::save_call(std::ostream&, const arguments&) const
{

}

//! \brief L2 norm to data.
double function::L2_distance(const ptr<data>& d) const
{
	double l2_dist = 0.0;
	for(int i=0; i<d->size(); ++i)
	{
		vec dat = d->get(i);
    vec x(parametrization().dimX());
    vec y(parametrization().dimY());

		if(parametrization().input_parametrization() == params::UNKNOWN_INPUT)
		{
			memcpy(&x[0], &dat[0], parametrization().dimX()*sizeof(double));
		}
		else
		{
      params::convert(&dat[0],
                      d->parametrization().input_parametrization(),
                      parametrization().input_parametrization(), &x[0]);
		}
    memcpy(&y[0], &dat[d->parametrization().dimX()],
           parametrization().dimY()*sizeof(double));

    l2_dist += std::pow(norm(y-value(x)), 2);
	}

  assert(d->size() > 0.);
	l2_dist = std::sqrt(l2_dist / d->size());

	return l2_dist;
}

//! \brief Linf norm to data.
double function::Linf_distance(const ptr<data>& d) const
{
	vec mean = vec::Zero(parametrization().dimY());
	vec var  = vec::Zero(parametrization().dimY());

#ifdef DEBUG
	std::cout << "<<DEBUG>> input param here = " << params::get_name(parametrization().input_parametrization()) << std::endl;
#endif

	double linf_dist = 0.0;
	for(int i=0; i<d->size(); ++i)
	{
		vec dat = d->get(i);
		vec x(parametrization().dimX()), y(parametrization().dimY());

    // Convert the position of the data sample to the parametrization
    // of the function.
    if(parametrization().input_parametrization() == params::UNKNOWN_INPUT)
    {
        memcpy(&x[0], &dat[0], parametrization().dimX()*sizeof(double));
    }
    else
    {
        params::convert(&dat[0],
                        d->parametrization().input_parametrization(),
                        parametrization().input_parametrization(), &x[0]);
    }

		// Copy the value part of the data vector in a vector to perform vector
		// operations on it (used in the computation of the mean).
    memcpy(&y[0], &dat[d->parametrization().dimX()],
           d->parametrization().dimY() * sizeof(double));

		// Take the componentwise-max of the two vectors.
		const vec v = value(x);
    for(int j=0; j<d->parametrization().dimY(); ++j)
		{
			linf_dist = std::max<double>(linf_dist, std::abs(y[j]-v[j]));
		}

		// Compute the mean
		mean += (y-v) / static_cast<double>(d->size());
	}
	
	// Compute the standard deviation with respect to the mean error
	for(int i=0; i<d->size(); ++i)
	{
		vec dat = d->get(i);
    vec x(parametrization().dimX()), y(d->parametrization().dimY()), val(parametrization().dimY());
		
        // Convert the position of the data sample to the parametrization
        // of the function.
        if(parametrization().input_parametrization() == params::UNKNOWN_INPUT)
        {
            memcpy(&x[0], &dat[0], parametrization().dimX()*sizeof(double));
        }
        else
        {
            params::convert(&dat[0],
                            d->parametrization().input_parametrization(),
                            parametrization().input_parametrization(),
                            &x[0]);
        }

		// Copy the value part of the data vector in a vector to perform vector
		// operations on it (used in the computation of the mean).
    memcpy(&y[0], &dat[d->parametrization().dimX()], parametrization().dimY()*sizeof(double));

		val = value(x);
    for(int j=0; j<d->parametrization().dimY(); ++j)
		{ 
      y[j] = dat[d->parametrization().dimX()+j];
			var[j] += pow(mean[j] - (val[j]-y[j]), 2) / double(d->size());
		}
	}

	std::cout << "<<INFO>> Mean = " << mean << ", Var = " << var << std::endl;

	return linf_dist;
}


/*--- Non-linear functions implementation ----*/

nonlinear_function::nonlinear_function()
{
}

bool nonlinear_function::load(std::istream& in)
{
	// Parse line until the next comment
	while(in.peek() != '#')
	{
		char line[256];
		in.getline(line, 256);

		// If we cross the end of the file, or the badbit is
		// set, the file cannot be loaded
		if(!in.good())
			return false;
	}

	// Checking for the comment line #FUNC nonlinear_function_phong
	std::string token;
	in >> token;
	if(token.compare("#FUNC") != 0) 
	{ 
		std::cerr << "<<ERROR>> parsing the stream. The #FUNC is not the next line defined." << std::endl; 
#ifdef DEBUG
		std::cout << "<<DEBUG>> got: \"" << token << "\"" << std::endl;
#endif
		return false;
	}

	in >> token;
	if(token.compare("nonlinear_function") != 0)
	{
		std::cerr << "<<ERROR>> parsing the stream. A function name is defined." << std::endl;
		std::cerr << "<<ERROR>> did you forget to specify the plugin used to export?" << std::endl;
		return false;
	}

	int nb_params = nbParameters();
	vec p(nb_params);
	for(int i=0; i<nb_params; ++i)
	{
		in >> token >> p[i];
	}

	setParameters(p);
	return true;
}

void nonlinear_function::save_call(std::ostream& out, const arguments& args) const
{
	if(!args.is_defined("export"))
	{
		// Dump a #FUNC nonlinear
		out << "#FUNC nonlinear_function" << std::endl;

		// Dump the parameters in order
		vec p = parameters();
		for(int i=0; i<p.size(); ++i)
		{
			out << "param_" << i+1 << "\t" << p[i] << std::endl;
		}
		out << std::endl;
	}

	function::save_call(out, args);
}

void nonlinear_function::bootstrap(const ptr<data> d, const arguments& args)
{
  #ifdef BOOTSTRAP_DEBUG
  std::cout << __FILE__ << " " << __LINE__ << " args = " << args << std::endl;
  #endif

  if(args.is_vec("bootstrap"))
  {
    vec p = args.get_vec("bootstrap", nbParameters());

    #ifdef BOOTSTRAP_DEBUG
    std::cout << __FILE__ << " " << __LINE__ << " args = " << args << std::endl;
    std::cout << "BOOTSTRAPPING VALUE = " << p << std::endl
              << " with nbParameters = " << nbParameters() << std::endl;
    #endif

    setParameters(p);
  }
  else
  {
    #ifdef BOOTSTRAP_DEBUG
    std::cout << __FILE__ << " " << __LINE__ << " args = " << args << std::endl;
    #endif

    function::bootstrap(d, args);
  }
}
		
vec nonlinear_function::getParametersMax() const
{
	vec M(nbParameters());
	for(int i=0; i<nbParameters(); ++i) 
	{
		M[i] = std::numeric_limits<double>::max();
	}
	return M;
}

vec nonlinear_function::getParametersMin() const
{
	vec m(nbParameters());
	for(int i=0; i<nbParameters(); ++i) 
	{
		m[i] = -std::numeric_limits<double>::max();
	}
	return m;
}



/*--- Compound functions implementation ----*/
compound_function::~compound_function()
{
}

vec compound_function::operator()(const vec& x) const
{
	return value(x);
}
vec compound_function::value(const vec& x) const
{
	vec res = vec::Zero(parametrization().dimY());
	for(unsigned int i=0; i<fs.size(); ++i)
	{
		vec temp_x(fs[i]->parametrization().dimX());
    params::convert(&x[0], parametrization().input_parametrization(),
                    fs[i]->parametrization().input_parametrization(), &temp_x[0]);
		res = res + fs[i]->value(temp_x);
	}
	return res;
}

vec compound_function::parametersJacobian(const vec& x) const
{
	int nb_params = nbParameters();
	vec jac = vec::Zero(nb_params*parametrization().dimY());

	int start_i = 0;

	// Export the sub-Jacobian for each function
	for(unsigned int f=0; f<fs.size(); ++f)
	{
		const ptr<nonlinear_function>& func = fs[f];
		int nb_f_params = func->nbParameters(); 

		// Only export Jacobian if there are non-linear parameters
		if(nb_f_params > 0 && !is_fixed[f])
		{

			vec temp_x(func->parametrization().dimX());
      params::convert(&x[0], parametrization().input_parametrization(),
                      func->parametrization().input_parametrization(), &temp_x[0]);
			vec func_jac = func->parametersJacobian(temp_x);

			for(int i=0; i<nb_f_params; ++i)
			{
				for(int y=0; y<parametrization().dimY(); ++y)
				{
					jac[y*nb_params + (i+start_i)] = func_jac[y*nb_f_params + i];
				}
			}

			start_i += nb_f_params;
		}
	}

	return jac;
}

// Return the parameters of the composition of FUNCTIONS.
static const parameters
compound_parameters(const std::vector<ptr<nonlinear_function> >& functions)
{
    assert(functions.size() > 0);

    parameters result = functions[0]->parametrization();

    for (auto f: functions)
    {
        auto input = f->parametrization().input_parametrization();
        auto output = f->parametrization().output_parametrization();

        if (result.input_parametrization() == params::UNKNOWN_INPUT)
            result = alta::parameters(result.dimX(), result.dimY(),
                                      input,
                                      result.output_parametrization());
        else if (result.input_parametrization() != input)
            // Parametrization mismatch: fall back to Cartesian.
            result = result.set_input(6, params::CARTESIAN);

        if (result.output_parametrization() == params::UNKNOWN_OUTPUT)
            result = alta::parameters(result.dimX(), result.dimY(),
                                      result.input_parametrization(),
                                      output);
        else if (result.output_parametrization() != output)
            abort();
    }

    return result;
}

compound_function::compound_function(const std::vector<ptr<nonlinear_function> >& functions,
                                     const std::vector<arguments> args)
    : nonlinear_function(compound_parameters(functions)),
      fs(functions), fs_args(args)
{
    assert(functions.size() == args.size());

    is_fixed.resize(functions.size());
    for (unsigned int x = 0; x < functions.size(); x++)
    {
        is_fixed[x] = args[x].is_defined("fixed");
    }
}


nonlinear_function* compound_function::operator[](int i) const
{
#ifdef DEBUG
	assert(i >= 0 && i < fs.size());
#endif
	return fs[i].get();
}
		
unsigned int compound_function::size() const
{
	return fs.size();
}

bool compound_function::load(std::istream& in)
{
	int nb_good = 0; // Number of correctly openned functions
	for(unsigned int i=0; i<fs.size(); ++i)
	{
		std::streampos pos = in.tellg();
		if(!fs[i]->load(in))
		{
			in.seekg(pos);
		}
		else
		{
			nb_good++;
		}
	}
	std::cout << "<<DEBUG>> number of correct function loaded in compound: " << nb_good  << " / " << fs.size() << std::endl;
	return nb_good > 0;
}

void compound_function::bootstrap(const ::ptr<data> d, const arguments& args)
{
	bool const global_bootstrap = args.is_defined("bootstrap");

	#ifdef BOOTSTRAP_DEBUG
	std::cout << __FILE__ << " " << __LINE__ << " args = " << args << std::endl;
	#endif

	// Check if the bootstrap is global
	if(global_bootstrap)
	{
		if(args.is_vec("bootstrap"))
		{
			vec p = args.get_vec("bootstrap", nbParameters());
			setParameters(p);
		
			for(unsigned int i=0; i<fs.size(); ++i)
			{
				if(fs[i]->nbParameters() == 0)
				{
					fs[i]->bootstrap(d, fs_args[i]);
				}
			}
      std::cout << "<<DEBUG>> Number of parameters for compound function: " << nbParameters() << std::endl;
      std::cout << "<<INFO>> Will use " << p << " as a bootstrap vector for the compound non-linear function" << std::endl;
		}
		else
		{
			std::ifstream file;
			file.open(args["bootstrap"].c_str()) ;
			if(file.is_open())
			{
				// Skip the header
				std::string line ;
				do
				{
					std::getline(file, line) ;
				}
				while(line != "#ALTA HEADER END");

				// Load the file
				// For each individual element, the input file is parsed
				// in order. if the function cannot parse it, the file is
				// restored to the previous state and the loading continues
				for(unsigned int i=0; i<fs.size(); ++i)
				{
					std::streampos pos = file.tellg();
					
					ptr<product_function> p = dynamic_pointer_cast<product_function>(fs[i]);
					if(p) {

						nonlinear_function* f1 = p->first();
						nonlinear_function* f2 = p->second();

						pos = file.tellg();
						if(!f1->load(file))
						{
							file.seekg(pos);

							// Bootstrap the function as if it was not loaded
							f1->bootstrap(d, args);
							std::cout << "<<DEBUG>> Unable to load first function of product, regular bootstraping" << std::endl;
						}

						// If the second function cannot be loaded, put the input stream
						// in the previous state and bootstrap normaly this function.
						pos = file.tellg();
						if(!f2->load(file))
						{
							file.seekg(pos);

							// Bootstrap the function as if it was not loaded
							f2->bootstrap(d, args);
							std::cout << "<<DEBUG>> Unable to load second function of product, regular bootstraping" << std::endl;
						}
					}
					else
					{	
						// If the function cannot be loaded, put the input stream
						// in the previous state and bootstrap normaly this function.
						if(!fs[i]->load(file))
						{
							file.seekg(pos);

							// Bootstrap the function as if it was not loaded
							fs[i]->bootstrap(d, fs_args[i]);
							std::cout << "<<DEBUG>> Unable to load one function of compound, regular bootstraping" << std::endl;
						}
					}//end of else
				}//end of for-loop
			}//end of if
			else
			{
				std::cerr << "<<ERROR>> you must provide a vector of parameters or a file to load with the bootstrap" << std::endl;
			}
		}
	}
	else
	{
		std::cout << "<<DEBUG>> No gobal bootstrap option. Bootstraping per function" << std::endl;
		for(unsigned int i=0; i<fs.size(); ++i)
		{
			fs[i]->bootstrap(d, fs_args[i]);
		}
	}
}

void compound_function::setMin(const vec& min) 
{
	function::setMin(min);
	for(unsigned int i=0; i<fs.size(); ++i)
	{
		fs[i]->setMin(min);
	}
}
void compound_function::setMax(const vec& max) 
{
	function::setMax(max);
	for(unsigned int i=0; i<fs.size(); ++i)
	{
		fs[i]->setMax(max);
	}
}

//! Number of parameters to this non-linear function
int compound_function::nbParameters() const
{
	int nb_params = 0;
	for(unsigned int i=0; i<fs.size(); ++i)
	{
		if(!is_fixed[i]) {
			nb_params += fs[i]->nbParameters();
		}
	}
	return nb_params;
}

//! Get the vector of parameters for the function
vec compound_function::parameters() const
{
	vec params(nbParameters());
	int current_i = 0;
	for(unsigned int f=0; f<fs.size(); ++f)
	{
		int f_size = fs[f]->nbParameters();

		// Handle when there is no parameters to include
		if(f_size > 0 && !is_fixed[f])
		{
			vec f_params = fs[f]->parameters();
			for(int i=0; i<f_size; ++i)
			{
				params[i + current_i] = f_params[i];
			}

			current_i += f_size;
		}
	}

	return params;
}

//! Get the vector of min parameters for the function
vec compound_function::getParametersMin() const
{
	vec params(nbParameters());
	int current_i = 0;
	for(unsigned int f=0; f<fs.size(); ++f)
	{
		int f_size = fs[f]->nbParameters();

		// Handle when there is no parameters to include
		if(f_size > 0 && !is_fixed[f])
		{
			vec f_params = fs[f]->getParametersMin();
			for(int i=0; i<f_size; ++i)
			{
				params[i + current_i] = f_params[i];
			}

			current_i += f_size;
		}
	}

	return params;
}

//! Get the vector of min parameters for the function
vec compound_function::getParametersMax() const
{
	vec params(nbParameters());
	int current_i = 0;
	for(unsigned int f=0; f<fs.size(); ++f)
	{
		int f_size = fs[f]->nbParameters();

		// Handle when there is no parameters to include
		if(f_size > 0 && !is_fixed[f])
		{
			vec f_params = fs[f]->getParametersMax();
			for(int i=0; i<f_size; ++i)
			{
				params[i + current_i] = f_params[i];
			}

			current_i += f_size;
		}
	}

	return params;
}

//! Update the vector of parameters for the function
void compound_function::setParameters(const vec& p) 
{
	int current_i = 0;
	for(unsigned int f=0; f<fs.size(); ++f)
	{
		int f_size = fs[f]->nbParameters();

		// Handle when there is no parameters to include
		if(f_size > 0 && !is_fixed[f])
		{
			vec f_params(f_size);
			for(int i=0; i<f_params.size(); ++i)
			{
				f_params[i] = p[i + current_i];
			}

			fs[f]->setParameters(f_params);
			current_i += f_size;
		}
	}
}
		
void compound_function::save_body(std::ostream& out, const arguments& args) const
{
	for(unsigned int i=0; i<fs.size(); ++i)
	{
		fs[i]->save_body(out, args);
		out << std::endl;
	}

	function::save_body(out, args);
}
		
void compound_function::save_call(std::ostream& out, const arguments& args) const
{
	bool is_cpp    = args["export"] == "C++";
	bool is_shader = args["export"] == "shader" || args["export"] == "explorer";
	bool is_matlab = args["export"] == "matlab";

	// This part is export specific. For ALTA, the coefficients are just
	// dumped as is with a #FUNC {plugin_name} header.
	//
	// For C++ export, the function call should be done before hand and
	// the line should look like:
	//   res += call_i(x);
	for(unsigned int i=0; i<fs.size(); ++i)
	{
		if(i != 0 && (is_cpp || is_matlab || is_shader))
		{
			out << "\tres += ";
		}

		fs[i]->save_call(out, args);

		if(is_cpp || is_matlab || is_shader)
		{
			out << ";" << std::endl;
		}
	}

	function::save_call(out, args);
}



// Return the parameters of the product of G1 and G2.
static const parameters
product_parameters(const ptr<nonlinear_function>& g1,
                   const ptr<nonlinear_function>& g2)
{
    parameters result = g1->parametrization();

    if(g1->parametrization().input_parametrization() !=
       g2->parametrization().input_parametrization())
    {
        result = alta::parameters(6, result.dimY(),
                                  params::CARTESIAN,
                                  result.output_parametrization());
    }

    return result;
}

/*--- Product functions implementation ----*/

product_function::product_function(const ptr<nonlinear_function>& g1,
                                   const ptr<nonlinear_function>& g2,
                                   bool is_g1_fixed, bool is_g2_fixed)
: nonlinear_function(product_parameters(g1, g2)),
  f1( g1 ),
	f2( g2 ),
	_is_fixed( std::pair<bool,bool>( is_g1_fixed, is_g2_fixed) )
{
}

product_function::~product_function()
{
}


vec product_function::operator()(const vec& x) const
{
	return value(x);
}
vec product_function::value(const vec& x) const
{
	// Convert input space to the 1rst function parametrization and compute the
	// output value
	vec xf1(f1->parametrization().dimX());
	params::convert(&x[0], parametrization().input_parametrization(),
                  f1->parametrization().input_parametrization(), &xf1[0]);
	vec f1res = f1->value(xf1);

	// Convert input space to the 2nd function parametrization and compute the
	// output value
	vec xf2(f2->parametrization().dimX());
	params::convert(&x[0], parametrization().input_parametrization(),
                  f2->parametrization().input_parametrization(), &xf2[0]);
	vec f2res = f2->value(xf2);

	vec res = product(f1res, f2res);

	//std::cout << f1res << " :::::::::: " << res << std::endl;
	//std::cout << f2res << std::endl;
	return res;
}
		
bool product_function::load(std::istream& in)
{
	bool loaded_f1 = false,loaded_f2 = false;
	std::streampos pos = in.tellg();

	// Load the first function
	if(f1) {
		loaded_f1 = f1->load(in);
		if(! loaded_f1)
		{
			in.seekg(pos);
		}
	}

	pos = in.tellg();
	// Load the second function
	if(f2) {
		loaded_f2 = f2->load(in);
		if(! loaded_f2)
		{
			in.seekg(pos);
		}
	}

	return loaded_f1 || loaded_f2;
}

void product_function::save_body(std::ostream& out, const arguments& args) const
{
	f1->save_body(out, args);
	out << std::endl;
	
	f2->save_body(out, args);
	out << std::endl;

	function::save_body(out, args);
}

void product_function::save_call(std::ostream& out, const arguments& args) const
{
	bool is_cpp    = args["export"] == "C++";
	bool is_shader = args["export"] == "shader" || args["export"] == "explorer";
	bool is_matlab = args["export"] == "matlab";

	// This part is export specific. For ALTA, the coefficients are just
	// dumped as is with a #FUNC {plugin_name} header.
	//
	// For C++ export, the function call should be done before hand and
	// the line should look like:
	//   res += call_i(x);
	if(is_cpp || is_matlab || is_shader)
	{
		out << "(";
	}

	f1->save_call(out, args);

	if(is_cpp || is_matlab || is_shader)
	{
		out << " * ";
	}

	f2->save_call(out, args);

	if(is_cpp || is_matlab || is_shader)
	{
		out << ")" << std::endl;
	}

	function::save_call(out, args);
}

void product_function::bootstrap(const ptr<data> d, const arguments& args)
{
	const bool global_bootstrap = args.is_defined("bootstrap");

	// Check if the bootstrap is global
	if(global_bootstrap)
	{
		if(args.is_vec("bootstrap"))
		{
			vec p = args.get_vec("bootstrap", nbParameters());
			setParameters(p);
		}
		else
		{
			std::ifstream file;
			file.open(args["bootstrap"].c_str()) ;
			if(file.is_open())
			{
				// Skip the header
				std::string line ;
				do
				{
					std::getline(file, line) ;
				}
				while(line != "#ALTA HEADER END");

				std::streampos pos = file.tellg();

				// If the first function cannot be loaded, put the input stream
				// in the previous state and bootstrap normaly this function.
				if(!f1->load(file))
				{
					file.seekg(pos);

					// Bootstrap the function as if it was not loaded
					f1->bootstrap(d, args);
					std::cout << "<<DEBUG>> Unable to load first function of product, regular bootstraping" << std::endl;
				}
				
				// If the second function cannot be loaded, put the input stream
				// in the previous state and bootstrap normaly this function.
				if(!f2->load(file))
				{
					file.seekg(pos);

					// Bootstrap the function as if it was not loaded
					f2->bootstrap(d, args);
					std::cout << "<<DEBUG>> Unable to load second function of product, regular bootstraping" << std::endl;
				}
			}
			else
			{
				std::cerr << "<<ERROR>> you must provide a vector of parameters or a file to load with the bootstrap" << std::endl;
			}
		}
	}
	else
	{
		f1->bootstrap(d, args);
		f2->bootstrap(d, args);
	}
}

void product_function::setMin(const vec& min) 
{
	function::setMin(min);
	f1->setMin(min);
	f2->setMin(min);
}

void product_function::setMax(const vec& max) 
{
	function::setMax(max);
	f1->setMax(max);
	f2->setMax(max);
}

int product_function::nbParameters() const
{
	int num_parameters = 0;
	
	num_parameters +=	(! _is_fixed.first) ? (f1->nbParameters()) : (0) ;

	num_parameters +=	(! _is_fixed.second) ? (f2->nbParameters()) : (0) ;

	return num_parameters;
}

vec product_function::parameters() const
{
	// int nb_f1_params = f1->nbParameters();
	// int nb_f2_params = f2->nbParameters();
	// int nb_params = nb_f1_params + nb_f2_params;
	int const nb_params = product_function::nbParameters();	

	int nb_f1_params = 0;
	
	vec params(nb_params);
	if( !_is_fixed.first )
	{
		//F1 is active save the number of paramters of f1 
		nb_f1_params = f1->nbParameters();

		vec f1_params = f1->parameters();
		for(int i=0; i<nb_f1_params; ++i)
		{
			params[i] = f1_params[i];
		}

	}

	if( ! _is_fixed.second )
	{
		int const nb_f2_params = f2->nbParameters();
		vec f2_params = f2->parameters();
		for(int i=0; i<nb_f2_params; ++i)
		{
			params[i+nb_f1_params] = f2_params[i];
		}
	}

	return params;
}

void product_function::setParameters(const vec& p)
{

	int nb_f1_params = 0;

	if( ! _is_fixed.first )
	{
		nb_f1_params = f1->nbParameters();
	
		vec f1_params(nb_f1_params);
		for(int i=0; i<nb_f1_params; ++i)
		{
			f1_params[i] = p[i];
		}

		f1->setParameters(f1_params);	    
	}

	if( ! _is_fixed.second )
	{
		int const nb_f2_params = f2->nbParameters();
		vec f2_params(nb_f2_params);
		for(int i=0; i<nb_f2_params; ++i)
		{
			f2_params[i] = p[i+nb_f1_params];
		}
		
		f2->setParameters(f2_params);			
	}


}

//! Get the vector of max parameters for the function
vec product_function::getParametersMax() const
{
	int const nb_params = product_function::nbParameters();	

	int nb_f1_params = 0;

	vec params(nb_params);

	if( !_is_fixed.first )
	{
		nb_f1_params = f1->nbParameters();
		vec f1_params = f1->getParametersMax();
		for(int i=0; i<nb_f1_params; ++i)
		{
			params[i] = f1_params[i];
		}
	}

	if( ! _is_fixed.second )
	{
		int const nb_f2_params = f2->nbParameters();
		vec f2_params = f2->getParametersMax();
		for(int i=0; i<nb_f2_params; ++i)
		{
			params[i+nb_f1_params] = f2_params[i];
		}

	}

	return params;
}

//! Get the vector of min parameters for the f1tion
vec product_function::getParametersMin() const
{
	int nb_f1_params = 0;

	int const nb_params = product_function::nbParameters();	

	vec params(nb_params);

	if( ! _is_fixed.first )
	{
		nb_f1_params = f1->nbParameters();

		vec f1_params = f1->getParametersMin();
		for(int i=0; i<nb_f1_params; ++i)
		{
			params[i] = f1_params[i];
		}
	}

	if( ! _is_fixed.second )
	{
		int const nb_f2_params = f2->nbParameters();
		vec f2_params = f2->getParametersMin();
		for(int i=0; i<nb_f2_params; ++i)
		{
			params[i+nb_f1_params] = f2_params[i];
		}	
	}


	

	return params;
}

vec product_function::parametersJacobian(const vec& x) const
{
//ß	std::cout << "ENTERING parametersJacobian " __FILE__ << " " << __LINE__ << std::endl;


	int const nb_f1_params = f1->nbParameters();
	int const nb_f2_params = f2->nbParameters();
	//int nb_params = nb_f1_params + nb_f2_params;
	int const nb_params = product_function::nbParameters();


	// Convert the input value x to the input space of the f1tion
	vec xf2(f2->parametrization().dimX());
	params::convert(&x[0], parametrization().input_parametrization(),
                  f2->parametrization().input_parametrization(), &xf2[0]);
	vec xf1(f1->parametrization().dimX());
	params::convert(&x[0], parametrization().input_parametrization(),
                  f1->parametrization().input_parametrization(), &xf1[0]);

	//Value of each function for the given x
	vec f1_value = f1->value(xf1);
	vec f2_value = f2->value(xf2);

	//Jacobien of each function for the given x
	vec f1_jacobian = f1->parametersJacobian(xf1);
	vec f2_jacobian = f2->parametersJacobian(xf2);

	// F = f2nel; f = f1tion
	// d(F * f)(x) /dp = F(x) df(x) /dp + f(x) dF(x) / dp
	//vec jac(nb_params*parametrization().dimY());
	
	//Forcing Zero by default
	vec jac= vec::Zero( nb_params * parametrization().dimY() );


	for(int y=0; y<parametrization().dimY(); ++y)
	{
		if( ! _is_fixed.first)
		{
			for(int i=0; i<nb_f1_params; ++i)
			{
				jac[y*nb_params + i] = f1_jacobian[y*nb_f1_params + i] * f2_value[y];
			}	
		}
		
		if( ! _is_fixed.second )
		{
			for(int i=0; i<nb_f2_params; ++i)
			{
				if( _is_fixed.first )
				{
					jac[y*nb_params + i ] = f2_jacobian[y*nb_f2_params + i] * f1_value[y];
				}
				else
				{
					jac[y*nb_params + (i+nb_f1_params)] = f2_jacobian[y*nb_f2_params + i] * f1_value[y];	
				}
				
			}				
		}
	}

	return jac;
}

nonlinear_function* product_function::first() const
{
	return f1.get();
}

nonlinear_function* product_function::second() const
{
	return f2.get();
}
