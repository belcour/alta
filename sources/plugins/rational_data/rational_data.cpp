#include "rational_data.h"

#include <boost/regex.hpp>
#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

void rational_data::load(const std::string& filename) 
{
	load(filename, -std::numeric_limits<double>::max(), std::numeric_limits<double>::max()) ;

}
		
void rational_data::load(const std::string& filename, double min, double max) 
{
	std::ifstream file(filename) ;
	if(!file.is_open())
	{
		std::cerr << "<<ERROR>> unable to open file \"" << filename << "\"" << std::endl ;
	}

	// N-Floats regexp
//	boost::regex e ("^([0-9]*\.?[0-9]+[\\t ]?)");

	_nX = 0 ; _nY = 0 ;
	std::vector<int> vs ;

	double x, y, dy ;
	while(file.good())
	{
		std::string line ;
		std::getline(file, line) ;
		std::stringstream linestream(line) ;
		
		// Discard incorrect lines
		if(linestream.peek() == '#')
		{
			linestream.ignore(1) ;

			std::string comment ;
			linestream >> comment ;
			std::cout << comment << std::endl ;

			if(comment == std::string("DIM"))
			{
				linestream >> _nX >> _nY ;

				vs.resize(dimY()) ;
				for(int k=0; k<dimY(); ++k)
				{
					vs[k] = 0 ;
				}

				_min.resize(dimX()) ;
				_max.resize(dimX()) ;
				for(int k=0; k<dimY(); ++k)
				{
					_min[k] =  std::numeric_limits<double>::max() ;
					_max[k] = -std::numeric_limits<double>::max() ;
				}
			}
			continue ;
		} 
		else if(line.empty()/*!boost::regex_match(line,e)*/)
		{
			continue ;
		}

		vec v ;
		v.resize(dimX() + 3*dimY()) ;
		for(int i=0; i<dimX(); ++i)
			linestream >> v[i] ;

		for(int i=0; i<dimY(); ++i)
			linestream >> v[dimX() + i] ;

		for(int i=0; i<dimY(); ++i)
		{
			// TODO, the firts case does not account for the
			// dimension of the ouput vector
			if(linestream.good()) 
			{
				linestream >> v[dimX() + dimY()+i] ;
			} 
			else 
			{
				// TODO Specify the delta in case
				// Handle multiple dim
				v[dimX() + dimY()+i] = v[dimX() + i] - 0.01f ;
				v[dimX() + 2*dimY()+i] = v[dimX() + i] + 0.01f ;
			}
		}
		
		// If data is not in the interval of fit
		// TODO: Update to more dims
		if(v[0] < min || v[0] > max)
		{
			continue ;
		}

		_data.push_back(v) ;

		// Update min and max
		for(int k=0; k<dimX(); ++k)
		{
			_min[k] = std::min(_min[k], v[k]) ;
			_max[k] = std::max(_max[k], v[k]) ;
		}
	}

	//TODO Test for small data
/*	std::vector<std::vector<double> > temp ;
	std::ofstream temp_out("temp.gnuplot", std::ios_base::trunc) ;
	for(int i=0; i<20; ++i)
	{
		int k = (i * _data.size()) / 20 ;
		temp_out << _data[k][0] << "\t" << _data[k][1] << "\t" << 0.5*(_data[k][2] - _data[k][1]) << std::endl ;
		temp.push_back(_data[k]) ;
	}
	_data = temp ;
*/
	// Sort the vector
//	std::sort(_data.begin(), _data.end(), [](const std::vector<double>& a, const std::vector<double>& b){return (a[0]<b[0]);});

	std::cout << "<<INFO>> loaded file \"" << filename << "\"" << std::endl ;
	std::cout << "<<INFO>> data inside " << _min << " ... " << _max << std::endl ;
	std::cout << "<<INFO>> loading data file of R^" << dimX() << " -> R^" << dimY() << std::endl ;
	std::cout << "<<INFO>> " << _data.size() << " elements to fit" << std::endl ;
}

bool rational_data::get(int i, double& x, double& yl, double& yu) const
{
	if(i >= (int)_data.size())
	{
		return false ;
	}

	x  = _data[i][0] ;
	yl = _data[i][1] ;
	yu = _data[i][2] ;

	return true ;
}
		
void rational_data::get(int i, vec& yl, vec& yu) const
{
	yl.resize(dimY()) ; yu.resize(dimY()) ;
	for(int j=0; j<dimY(); ++j)	
	{
		yl[j] = _data[i][dimX() + dimY() + j] ;
		yu[j] = _data[i][dimX() + 2*dimY() + j] ;
	}
}
		
const vec& rational_data::operator[](int i) const
{
	return _data[i] ;
}
const vec& rational_data::get(int i) const 
{
	return _data[i] ;
}

int rational_data::size() const
{
	return _data.size() ;
}

vec rational_data::min() const 
{
	return _min ;
}

vec rational_data::max() const 
{
	return _max ;
}

Q_EXPORT_PLUGIN2(rational_data, rational_data)
