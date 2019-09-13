# Copyright (c) 2011-2018, The DART development contributors
# All rights reserved.
#
# The list of contributors can be found at:
#   https://github.com/dartsim/dart/blob/master/LICENSE
#
# This file is provided under the "BSD-style" License

# Find NLOPT
#
# This sets the following variables:
# NLOPT_FOUND
# NLOPT_INCLUDE_DIRS
# NLOPT_LIBRARIES
# NLOPT_DEFINITIONS
# NLOPT_VERSION

# 2018 : SMALL MODIFICATIONS FROM R.P. romain dot pacanowski @ institutoptique DOT fr
# NLOPT_LIB_DIRS : the paths to the nlopt libraries
# SETTING A version even if pkg_config is not available

find_package(PkgConfig QUIET)

# Check to see if pkgconfig is installed.
pkg_check_modules(PC_NLOPT nlopt QUIET)

# Definitions
set(NLOPT_DEFINITIONS ${PC_NLOPT_CFLAGS_OTHER})

# Include directories
find_path(NLOPT_INCLUDE_DIRS
    NAMES nlopt.h
    HINTS ${PC_NLOPT_INCLUDEDIR}
    PATHS "${CMAKE_INSTALL_PREFIX}/include" "C:/Program Files/NLOPT/include" )

# Libraries
find_library(NLOPT_LIBRARIES
    NAMES nlopt nlopt_cxx
    HINTS ${PC_NLOPT_LIBDIR} "C:/Program Files/NLOPT/lib")

# Add-on from RP.
# Libraries Path
find_path(NLOPT_LIB_DIRS
    NAMES nlopt.lib nlopt.dll nlopt-0.dll
    HINTS ${PC_NLOPT_LIBDIR} 
    PATHS "C:/Program Files/NLOPT/lib/" )


# Version
if(DEFINED ${PC_NLOPT_VERSION})
  set(NLOPT_VERSION ${PC_NLOPT_VERSION})
else()
  set(NLOPT_VERSION "VERSION NOT FOUND")
endif()


# Set (NAME)_FOUND if all the variables and the version are satisfied.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NLOPT
    FAIL_MESSAGE  DEFAULT_MSG
    REQUIRED_VARS NLOPT_INCLUDE_DIRS NLOPT_LIBRARIES
VERSION_VAR NLOPT_VERSION)