import os, sys

##----------------------------------------------------------------##
## This file describes required and optional arguments to ALTA    ##
## compilation. If you want to manualy tune the use of an already ##
## present library, uncomment the according lines.                ##
##----------------------------------------------------------------##


## Compiler and build flags.
##
#CXX            = 'g++'


## Python and boost-python library
##
#PYTHON_INC    = ['/usr/include/python2.7']
#PYTHON_DIR    = ['/usr/lib/x86_64-linux-gnu']
#PYTHON_LIB    = ['boost_python-py27', 'python2.7']


## Eigen library
##
EIGEN_INC     = ['#external/build/include/Eigen']


## OpenEXR library
##
#OPENEXR_INC    = ['/usr/include/OpenEXR']
#OPENEXR_DIR    = ['/usr/lib/x86_64_linux-gnu']
#OPENEXR_LIB    = ['Half', 'IlmImf', 'IlmThread']


## FLANN library
##
FLANN_INC    = ['/usr/include/flann']
FLANN_DIR    = ['/usr/lib/x86_64_linux-gnu']
FLANN_LIB    = ['flann']


## QUADPROG library
##
## You have to specify the directory of the QuadProg library
##
QUADPROG_INC      = ['#external/quadprog++']
QUADPROG_DIR      = ['#external/build/lib']
QUADPROG_LIB      = ['quadprog++']


## CERES library
##
## There is no 'ceres' or 'glog' packages available on Debian. If they are not
## provided by the user, the automatic installation tool will be runnned.
##
CERES_INC      = ['#external/build/include']
CERES_DIR      = ['#external/build/lib']
CERES_LIB      = ['ceres', 'glog']
CERES_OPT_LIB  = ['gomp', 'lapack', 'blas', 'amd', 'camd', 'ccolamd', 'colamd', 'cholmod', 'cxsparse']


## NlOpt library
##
## NlOpt has a default pkgconfig configuration file. If you want to use you
## own NlOpt installation with no support of pkgconfig (not advised), please
## uncomment and update the following variables.
##
#NLOPT_INC      = ['#external/build/include']
#NLOPT_DIR      = ['#external/build/lib']
#NLOPT_LIB      = ['nlopt']
#NLOPT_OPT_LIB  = []


## coin IpOpt library
##
## IpOpt has a default pkgconfig configuration file. If you want to use you
## own IpOpt installation with no support of pkgconfig (not advised), please
## uncomment and update the following variables.
##
#IPOPT_INC      = ['#external/build/include']
#IPOPT_DIR      = ['#external/build/lib']
#IPOPT_LIB      = ['ipopt']
#IPOPT_OPT_LIB  = []


## MATLAB library and Engine
##
## The Directory where engine.h is located
MATLAB_INC  = []
## The Directory where libeng.so  libmex.so libmat.so (or .dll) libraries are located
MATLAB_DIR  = []
# Matlb Libraries needed at compilation and runtime
MATLAB_LIB  = ['eng', 'mex','mat']
