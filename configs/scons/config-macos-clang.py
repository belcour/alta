import os, sys

##----------------------------------------------------------------##
## This file describes required and optional arguments to ALTA    ##
## compilation. If you want to manualy tune the use of an already ##
## present library, update the according lines                    ##
##                                                                ##
## This configuration file is made of a MacOS X operating system  ##
## version 10.9, with clang and OpenEXR and FLANN libraries       ##
## installed using Ports.                                         ##
##----------------------------------------------------------------##


## Compiler and build flags.
##
#CXX            = 'clang++'
CCFLAGS        = '-std=c++11 -O1 -g -Wall'


## Python and boost-python library
##
## Be sure to install boost-python and python
## You can do this by using Macports (http://www.macports.org/) 
## or Homebrew (http://brew.sh/)
##
# PYTHON_INC    = ['/opt/local/include', '/usr/include/python2.7']
# PYTHON_DIR    = ['/opt/local/lib']
# PYTHON_LIB    = ['python2.7']


## Eigen library
##
EIGEN_INC     = ['#external/build/include']

## OpenEXR library
##
## OpenEXR in ports has a default pkgconfig file. If you want to specify your
## own compilation of OpenEXR, please uncomment and update the two variables.
##
#OPENEXR_INC    = ['/opt/local/include/OpenEXR']
#OPENEXR_DIR    = ['/opt/local/lib']


## FLANN library
##
FLANN_INC    = ['/opt/local/include']
FLANN_DIR    = ['/opt/local/lib']
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
## On OSX with ports, 'ceres' and 'glog' packages are available. By defaults
## we try to get them. If they are not provided by the user, the automatic
## installation tool will be runnned.
##
CERES_INC      = ['/opt/local/include']
CERES_DIR      = ['/opt/local/lib']
CERES_LIB      = ['ceres', 'glog']
CERES_OPT_LIB  = ['lapack', 'blas', 'amd', 'camd', 'ccolamd', 'colamd', 'cholmod', 'cxsparse', 'atlas']


## NlOpt library
##
## NlOpt has a default pkgconfig configuration file. If you want to use you
## own NlOpt installation with no support of pkgconfig (not advised), please
## uncomment and update the following variables.
##
# NLOPT_INC      = ['#external/build/include']
# NLOPT_DIR      = ['#external/build/lib']
# NLOPT_LIB      = ['nlopt']
# NLOPT_OPT_LIB  = []


## coin IpOpt library
##
## IpOpt has a default pkgconfig configuration file. If you want to use you
## own IpOpt installation with no support of pkgconfig (not advised), please
## uncomment and update the following variables.
##
# IPOPT_INC      = ['#external/build/include']
# IPOPT_DIR      = ['#external/build/lib']
# IPOPT_LIB      = ['ipopt']
# IPOPT_OPT_LIB  = []


## MATLAB library and Engine
##
#MATLAB_INC  = ['/Applications/MATLAB_R2015b.app/extern/include/']
#MATLAB_DIR  = ['/Applications/MATLAB_R2015b.app/bin/maci64/']
#MATLAB_LIB  = ['eng', 'mx','mat']
