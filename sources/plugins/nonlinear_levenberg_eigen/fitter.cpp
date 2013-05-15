#include "fitter.h"

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <unsupported/Eigen/LevenbergMarquardt>

#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

#include <QTime>

fitter* provide_fitter()
{
    return new nonlinear_fitter_eigen();
}

struct EigenFunctor: Eigen::DenseFunctor<double>
{
    EigenFunctor(nonlinear_function* f, const data* d) : 
		 DenseFunctor<double>(f->nbParameters(), d->size()), _f(f), _d(d)
	{
#ifndef DEBUG
		std::cout << "<<DEBUG>> constructing an EigenFunctor for n=" << inputs() << " parameters and m=" << values() << " points" << std::endl ;
#endif
	}

	int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& y) const
	{
		std::cout << "input vector = " << x << std::endl;


		// Update the parameters vector
		vec _p(inputs());
		for(int i=0; i<inputs(); ++i)
		{
			_p[i] = x(i);
		}
		_f->setParameters(_p);

		for(int s=0; s<values(); ++s)
		{
			vec _x  = _d->get(s);

			vec _di = vec(_f->dimY());
			for(int i=0; i<_f->dimY(); ++i)
				_di[i] = _x[_f->dimX() + i];

			// Should add the resulting vector completely
			vec _y = _di - (*_f)(_x);
			y(s) = _y[0];
		}
		
		return 0;
	}

	int df(const Eigen::VectorXd &x, Eigen::MatrixXd &fjac) const
	{
		vec _p(inputs());
		for(int i=0; i<inputs(); ++i) _p[i] = x(i);
		_f->setParameters(_p);

		for(int i=0; i<values(); ++i)
		{
			vec xi = _d->get(i);

			vec _jac = _f->parametersJacobian(xi);
			for(int j=0; j<_f->nbParameters(); ++j)
			{
				fjac(i,j) = -_jac[j];
			}
		}

		return 0;
	}


	nonlinear_function* _f;
    const data* _d;
};

nonlinear_fitter_eigen::nonlinear_fitter_eigen() : QObject()
{
}
nonlinear_fitter_eigen::~nonlinear_fitter_eigen() 
{
}

bool nonlinear_fitter_eigen::fit_data(const data* d, function* fit, const arguments &args)
{
    // I need to set the dimension of the resulting function to be equal
    // to the dimension of my fitting problem
    fit->setDimX(d->dimX()) ;
    fit->setDimY(d->dimY()) ;
    fit->setMin(d->min()) ;
    fit->setMax(d->max()) ;

    if(dynamic_cast<nonlinear_function*>(fit) == NULL)
    {
        std::cerr << "<<ERROR>> the function is not a non-linear function" << std::endl;
        return false;
    }
    nonlinear_function* nf = dynamic_cast<nonlinear_function*>(fit);

#ifndef DEBUG
	 std::cout << "<<DEBUG>> number of parameters: " << nf->nbParameters() << std::endl;
#endif

    /* the following starting values provide a rough fit. */
    int info;
    Eigen::VectorXd x;
    x.setConstant(nf->nbParameters(), 1.);

    EigenFunctor functor(nf, d);
    Eigen::LevenbergMarquardt<EigenFunctor> lm(functor);

	 info = lm.minimize(x);

    if(info == Eigen::LevenbergMarquardtSpace::ImproperInputParameters)
    {
        std::cerr << "<<ERROR>> incorrect parameters" << std::endl;
        return false;
    }


    vec p(nf->nbParameters());

    for(int i=0; i<p.size(); ++i)
    {
        p[i] = x(i);
    }
    std::cout << "<<INFO>> found parameters: " << p << std::endl;
#ifndef DEBUG
    std::cout << "<<DEBUG>> using " << lm.iterations() << " iterations" << std::endl;
#endif
    nf->setParameters(p);
    return true;

}

void nonlinear_fitter_eigen::set_parameters(const arguments& args)
{
}
		
data*     nonlinear_fitter_eigen::provide_data() const { return NULL; }
function* nonlinear_fitter_eigen::provide_function() const { return NULL; }

Q_EXPORT_PLUGIN2(nonlinear_fitter_eigen, nonlinear_fitter_eigen)
