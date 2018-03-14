import os, sys

##----------------------------------------------------------------##
## This file describes required and optional arguments to ALTA    ##
## compilation. If you want to manualy tune the use of an already ##
## present library, uncomment the according lines.                ##
##----------------------------------------------------------------##


## Compilators build flags
##
CXX            = 'cl'

# We used to have '/arch:AVX' here, but that generates invalid code if
# the underlying CPU does not support AVX (which is typically the case
# in VMs.)
CCFLAGS        = '/Zi /nologo /O2 /Ox /EHsc /MD'
#LINKFLAGS      = '/MACHINE:X86'


## Python and boost-python library
##
PYTHON_LIB    = ['python27', 'boost_python-vc140-mt-1_63']
PYTHON_INC    = ['C:/Python27/include/', 'C:/Libraries/boost_1_63_0']
PYTHON_DIR    = ['C:/Python27/libs/', 'C:/Libraries/boost_1_63_0/lib32-msvc-14.0']


## Eigen library
##
EIGEN_INC     = ['#external/build/include/Eigen']


## OpenMP flags
##
OPENMP_CFLAGS    = ' /openmp'
OPENMP_LINKFLAGS = ''


## OpenEXR library
##
OPENEXR_INC    = ['']
OPENEXR_DIR    = ['']


## FLANN library
##
FLANN_INC    = ['']
FLANN_DIR    = ['']
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
## You have to specify both the directory of the CERES library
## and the glog library
##
CERES_INC      = ['#external/build/include', '#external/build/include/ceres/internal/miniglog']
CERES_DIR      = ['#external/build/lib']
CERES_LIB      = ['ceres']
CERES_OPT_LIB  = ['glog', 'gomp', 'lapack', 'blas']


## NlOpt library
##
## You have to specify the directory of the NlOpt library
##
NLOPT_INC      = ['#external/build/include']
NLOPT_DIR      = ['#external/build/lib']
NLOPT_LIB      = ['nlopt']
NLOPT_OPT_LIB  = []


## coin IpOpt library
##
## You have to specify the directory of the IpOpt library
##
IPOPT_INC      = ['#external/build/include']
IPOPT_DIR      = ['#external/build/lib']
IPOPT_LIB      = ['ipopt']
IPOPT_OPT_LIB  = []


## MATLAB library and Engine
##
#MATLAB_INC  = ['']
#MATLAB_DIR  = ['']
#MATLAB_LIB  = ['eng', 'mex','mat']
