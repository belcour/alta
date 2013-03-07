#include <iostream>
#include <fstream>
#include <cmath>

#include <core/args.h>

int main(int argc, char** argv)
{
	std::ofstream f("input.gnuplot") ;
	arguments args(argc, argv);	

	int nbx = 100;
	int nby = 100;
	int nbz = 100;
	if(args.is_defined("nbx"))
		nbx = args.get_int("nbx", 100) ;
	if(args.is_defined("nby"))
		nby = args.get_int("nby", 100) ;

	const int k = args.get_int("f", 1) ;
	if(k == 1)
	{
		f << "#DIM 1 1" << std::endl ;
		for(int i=0; i<nbx; ++i)
		{
			const float x = i / (float)nbx ;
			const float y = 100.0f * exp(-10.0 * x*x) * x*x - 0.01 *x*x*x ;

			f << x << "\t" << y << "\t" << 0.1f << std::endl ;
		}
	}
	else if(k == 2)
	{
		f << "#DIM 2 1" << std::endl ;
		for(int i=0; i<nbx; ++i)
			for(int j=0; j<nby; ++j)
			{
				const float x = i / (float)nbx ;				
				const float y = j / (float)nby ; 
				const float z = 1 + 0.1f*x;
			
				f << x << "\t" << y << "\t" << z << "\t" << 0.1f << std::endl ;
			}
	}
	else if(k == 3)
	{
		f << "#DIM 2 1" << std::endl ;
		for(int i=0; i<nbx; ++i)
			for(int j=0; j<nby; ++j)
			{
				const float x = i / (float)nbx ;
				const float y = j / (float)nby ; 
				const float z = exp(-10.0 * x*x) +  x*y ;
			
				f << x << "\t" << y << "\t" << z << "\t" << 0.1f << std::endl ;
			}
	}

	return 0 ;
}