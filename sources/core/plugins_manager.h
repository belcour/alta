#pragma once 

#include <map>
#include <string>

#include "function.h"
#include "data.h"
#include "fitter.h"
#include "args.h"

/*! \brief This class handles the loading of plugins and insure that they can
 *  talk to each others through coordinates transforms.
 *
 *  \details
 *
 *  \todo Should we put this class as a singleton ? I am tempted by it so that
 *  each plugin can access some informations.
 */
class plugins_manager
{
	public: //functions

		//! \brief Create the object, parse the argument and load all the plugins
		plugins_manager(const arguments& args) ;

		//! \brief Get instances of the function, the data and the fitter. Select 
		//! the first in the map,
        static function* get_function() ;
		data*     get_data()     const ;
		fitter*   get_fitter()   const ;
		
		//! \brief Get instances of the function, the data and the fitter, select one 
		//! based on the name. Return null if no one exist.
        static function* get_function(const std::string& n) ;
        static data*     get_data(const std::string& n) ;
        fitter*   get_fitter(const std::string& n)   const ;

        //! \brief check if a data object and a function object are compatibles.
        //! this has to be done before fitting to ensure that the
        //! parametrizations spaces are the same.
        //! \todo specify an output parametrization for the function ?
        static void check_compatibility(data*& d, function*& f,
                                        const arguments& args)
        {
            if(f->parametrization() == params::UNKNOWN_INPUT)
            {
                std::cout << "<<DEBUG>> function will take the parametrization of the data" << std::endl;
                f->setParametrization(d->parametrization());
            }
            else if(d->parametrization() != f->parametrization())
            {
                std::cout << "<<INFO>> has to change the parametrization of the input data" << std::endl;
                data_params* dd = new data_params(d, f->parametrization());
                d = dd ;
            }
            else
            {
                std::cout << "<<DEBUG>> no change was made to the parametrization" << std::endl;
            }
        }

		//! \brief Provide a measure of how much memory there is on the system.
		//! \details It permits to know is one can allocate more memory for a fitting
		//! procedure for example.
		static size_t get_system_memory() ;

	private: //data

		std::map<std::string, function*> _functions ;
		std::map<std::string, data*>     _datas ;
		std::map<std::string, fitter*>   _fitters ;
} ;
