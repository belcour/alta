cmake_minimum_required (VERSION 3.2)
project (alta CXX)

#############################
#      Preconfiguration     #
#############################

# Set Python version number
set(Python_Version "3" CACHE STRING "Python interface version (default Python 3)")
set(Python_ADDITIONAL_VERSIONS $(Python_Version))

# Require C++11 options
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


#############################
#        Dependencies       #
#############################
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/configs/cmake/)

# Find required packages
find_package(Eigen3)
find_package(OpenMP)
find_package(PythonLibs)
find_package(Doxygen)
find_package(Ceres)
find_package(NLOPT)

add_definitions(-DALTA_PLUGIN_DIRECTORY="/usr/lib/alta_plugins")

# Add core and library include files
include_directories("external/build" "sources" ${EIGEN3_INCLUDE_DIR})
include_directories("external" "external/Catch/include" "external/quadprog++")
include_directories("external/pybind11/include")

# Look for header only dependencies
find_file(EIGEN_FOUND  "Eigen/Core"      HINTS ${EIGEN3_INCLUDE_DIR})
find_file(CATCH_FOUND  "catch.hpp"       HINTS external/Catch/include)
find_file(FLANN_FOUND  "flann/flann.hpp")
find_file(PYBIND_FOUND "pybind11/pybind11.h" HINTS "external/pybind11/include")

# Check if Eigen is found
if(NOT EIGEN_FOUND)
    message(FATAL_ERROR "Unable to find Eigen header")
endif()

# Update compilation option to include OpenMP if present
if(OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS ${OpenMP_CXX_FLAGS})
endif(OPENMP_FOUND)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)


#############################
#      Helper functions     #
#############################

function(alta_add_plugin name sources)
	add_library(${name} SHARED sources/plugins/${sources})
        target_link_libraries(${name} core)
endfunction(alta_add_plugin)

function(alta_add_soft name sources)
	add_executable(${name} sources/softs/${sources})
        target_link_libraries(${name} core)
endfunction(alta_add_soft)

function(alta_test)
    set(options OPTIONAL)
    set(oneValueArgs NAME)
    set(multiValueArgs COMMAND)
    cmake_parse_arguments(TEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_COMMAND}
             WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/tests")
    set_tests_properties("${TEST_NAME}"
                         PROPERTIES ENVIRONMENT "ALTA_PLUGIN_PATH=${CMAKE_BINARY_DIR}/plugins")  
endfunction(alta_test)

function(alta_test_python)
    set(options OPTIONAL)
    set(oneValueArgs NAME)
    set(multiValueArgs COMMAND)
    cmake_parse_arguments(TEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_COMMAND}
             WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/tests")
    set_tests_properties("${TEST_NAME}"
        PROPERTIES ENVIRONMENT "ALTA_PLUGIN_PATH=${CMAKE_BINARY_DIR}/plugins;PYTHONPATH=${CMAKE_BINARY_DIR}/python")  
endfunction(alta_test_function)

function(alta_test_unit name sources)
        add_executable(${name} sources/tests/${sources})
        target_link_libraries(${name} core)
        add_test(NAME ${name} COMMAND ${CMAKE_BINARY_DIR}/tests/${name} WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
        set_tests_properties(${name} PROPERTIES ENVIRONMENT "TEST_DATA_DIRECTORY=${CMAKE_SOURCE_DIR}/sources/tests;ALTA_PLUGIN_PATH=${CMAKE_BINARY_DIR}/plugins")
endfunction(alta_test_unit)

function(alta_test_soft name)
        add_test(NAME ${name} COMMAND ${CMAKE_BINARY_DIR}/softs/${name} --help)
endfunction(alta_test_soft)

function(alta_download_merl name)
    find_file(${name}_FOUND "${name}.binary" PATHS "${CMAKE_BINARY_DIR}/tests/")
    if(NOT ${name}_FOUND)
        message("Dowloading ${name} from MERL database:")
	file(DOWNLOAD
	     "http://people.csail.mit.edu/wojciech/BRDFDatabase/brdfs/${name}.binary"
	     "${CMAKE_BINARY_DIR}/tests/${name}.binary" SHOW_PROGRESS)
	set(${name}_FOUND TRUE PARENT_SCOPE)
    endif()
endfunction(alta_download_merl)


#############################
#         ALTA core         #
#############################

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

add_library(core STATIC
            sources/core/args.h
            sources/core/common.h
            sources/core/common.cpp
            sources/core/metrics.h
            sources/core/metrics.cpp
            sources/core/params.h
            sources/core/params.cpp
            sources/core/data.h
            sources/core/data.cpp
            sources/core/data_storage.h
            sources/core/data_storage.cpp
            sources/core/vertical_segment.h
            sources/core/vertical_segment.cpp
            sources/core/function.h
            sources/core/function.cpp
            sources/core/rational_function.h
            sources/core/rational_function.cpp
            sources/core/plugins_manager.h
            sources/core/plugins_manager.cpp)


add_library(quadprog STATIC
            external/quadprog++/QuadProg++.hh
            external/quadprog++/QuadProg++.cc)

target_link_libraries(core ${CMAKE_DL_LIBS})


#############################
#         Plugins           #
#############################

#    set(CMAKE_SHARED_LIBRARY_PREFIX "")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/plugins)

# Data
alta_add_plugin(data_astm       data_io/astm.cpp)
alta_add_plugin(data_merl       data_io/merl.cpp)
alta_add_plugin(data_utia       data_io/utia.cpp)
alta_add_plugin(data_brdf_slice data_io/slice.cpp)

alta_add_plugin(data_grid       data_interpolants/grid.cpp)
if(FLANN_FOUND)
    alta_add_plugin(data_rbf    data_interpolants/rbf.cpp)
endif()
#alta_add_plugin(data_matlab data_interpolants/matlab.cpp)

# Functions
alta_add_plugin(rational_function_legendre		rational_function_legendre/rational_function.cpp)
alta_add_plugin(rational_function_chebychev		rational_function_chebychev/rational_function.cpp)
alta_add_plugin(rational_function_cosine 		rational_function_cosine/rational_function.cpp)
alta_add_plugin(nonlinear_function_abc                  nonlinear_function_abc/function.cpp)
alta_add_plugin(nonlinear_function_beckmann             nonlinear_function_beckmann/function.cpp)
alta_add_plugin(nonlinear_function_blinn                nonlinear_function_blinn/function.cpp)
alta_add_plugin(nonlinear_function_diffuse              nonlinear_function_diffuse/function.cpp)
alta_add_plugin(nonlinear_function_lafortune            nonlinear_function_lafortune/function.cpp)
alta_add_plugin(nonlinear_function_sgd                  nonlinear_function_sgd/function.cpp)
alta_add_plugin(nonlinear_function_spherical_gaussian   nonlinear_function_spherical_gaussian/function.cpp)
alta_add_plugin(nonlinear_function_ward                 nonlinear_function_ward/function.cpp)

alta_add_plugin(nonlinear_function_retrobeckmann        nonlinear_function_retrobeckmann/function.cpp)
alta_add_plugin(nonlinear_function_retroblinn		nonlinear_function_retroblinn/function.cpp)
alta_add_plugin(nonlinear_function_retroyoo             nonlinear_function_retroyoo/function.cpp)
alta_add_plugin(nonlinear_fresnel_retroschlick		nonlinear_fresnel_retroschlick/function.cpp)

alta_add_plugin(nonlinear_fresnel_schlick               nonlinear_fresnel_schlick/function.cpp)
alta_add_plugin(nonlinear_fresnel_normalized_schlick    nonlinear_fresnel_normalized_schlick/function.cpp)

alta_add_plugin(nonlinear_shadowing_schlick             nonlinear_shadowing_schlick/function.cpp)
alta_add_plugin(nonlinear_shadowing_smith               nonlinear_shadowing_smith/function.cpp)
alta_add_plugin(nonlinear_shadowing_walter_smith        nonlinear_shadowing_walter_smith/function.cpp)

# Fitters
alta_add_plugin(rational_fitter_eigen 			rational_fitter_eigen/rational_fitter.cpp)
alta_add_plugin(rational_fitter_leastsquare		rational_fitter_leastsquare/rational_fitter.cpp)
alta_add_plugin(rational_fitter_quadprog		rational_fitter_quadprog/rational_fitter.cpp)
alta_add_plugin(rational_fitter_parallel		rational_fitter_parallel/rational_fitter.cpp)
alta_add_plugin(nonlinear_fitter_eigen                  nonlinear_fitter_eigen/fitter.cpp)
target_link_libraries(rational_fitter_quadprog quadprog)
target_link_libraries(rational_fitter_parallel quadprog)


# TODO: Add check before compiling  IPOPT
if (CERES_FOUND)
    include_directories( ${CERES_INCLUDE_DIR} )
    alta_add_plugin(nonlinear_fitter_ceres                  nonlinear_fitter_ceres/fitter.cpp)
    target_link_libraries(nonlinear_fitter_ceres ${CERES_LIBRARIES})
    list(APPEND nonlinear_fitters ceres)
endif()
if (NLOPT_FOUND)
    include_directories( ${NLOPT_INCLUDE_DIRS} )
    alta_add_plugin(nonlinear_fitter_nlopt                  nonlinear_fitter_nlopt/fitter.cpp)
    target_link_libraries(nonlinear_fitter_nlopt ${NLOPT_LIBRARIES})
    list(APPEND nonlinear_fitters nlopt)
endif()
#alta_add_plugin(nonlinear_fitter_ipopt                  nonlinear_fitter_ipopt/fitter.cpp)

# Python bindings
if(PYTHONLIBS_FOUND AND PYBIND_FOUND)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python)
    include_directories(${PYTHON_INCLUDE_DIRS})
    add_library(alta SHARED sources/python/alta.cpp)
    target_link_libraries(alta core ${PYTHON_LIBRARIES})
    set_target_properties(alta PROPERTIES PREFIX "")
    if(APPLE)
        set_target_properties(alta PROPERTIES SUFFIX ".so")
    endif()
endif()


#############################
#          Softs            #
#############################

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/softs)

alta_add_soft(brdf2brdf    brdf2brdf/main.cpp)
alta_add_soft(brdf2data    brdf2data/main.cpp)
alta_add_soft(brdf2gnuplot brdf2gnuplot/main.cpp)
#alta_add_soft(brdf2moments brdf2moments/main.cpp)
alta_add_soft(brdf2stats   fit2stat/fit2stat.cpp)
alta_add_soft(data2data    data2data/main.cpp)
alta_add_soft(data2brdf    data2brdf/main.cpp)
alta_add_soft(data2stats   data2stats/data2stats.cpp)
alta_add_soft(data2moments data2moments/main.cpp)


#############################
#          Tests            #
#############################

enable_testing()

# Test all softs with `--help` option
alta_test_soft(brdf2data)
alta_test_soft(data2data)
alta_test_soft(data2stats)
alta_test_soft(data2brdf)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests)
include_directories("sources/tests")

# Unit test
alta_test_unit(test_data_io  core/data-io.cpp)
alta_test_unit(half-test-1   core/half-test-1.cpp)
alta_test_unit(half-test-2   core/half-test-2.cpp)
alta_test_unit(half-test-3   core/half-test-3.cpp)
alta_test_unit(half-test-4   core/half-test-4.cpp)
alta_test_unit(nonlinear-fit core/nonlinear-fit.cpp)
alta_test_unit(params-test-1 core/params-test-1.cpp)
alta_test_unit(params-test-2 core/params-test-2.cpp)

if(CPPQUICKCHECK_FOUND)
    alta_test_unit(params-qc-1 core/params-qc-1.cpp)
endif()

# Integration test for rational function fitting
foreach(fitter IN ITEMS eigen leastsquare quadprog parallel)
    add_test(NAME "data2dbrdf_kirby_${fitter}"
             COMMAND "data2brdf" "--input"   "${CMAKE_SOURCE_DIR}/sources/tests/Kirby2.dat"
                                 "--output"  "Kirby2.func"
                                 "--fitter"  "rational_fitter_${fitter}"
             WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/tests")

    set_tests_properties("data2dbrdf_kirby_${fitter}"
                         PROPERTIES ENVIRONMENT "ALTA_PLUGIN_PATH=${CMAKE_BINARY_DIR}/plugins")
endforeach()

add_test(NAME "brdf2data_kirby"
         COMMAND "brdf2data" "--input"     "Kirby2.func"
                             "--output"    "Kirby2.dat"
                             "--data-file" "${CMAKE_SOURCE_DIR}/sources/tests/Kirby2.dat"
         WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/tests")

add_test(NAME "data2stats_kirby"
         COMMAND "data2stats" "--input"     "Kirby2.dat"
                              "--ref"       "${CMAKE_SOURCE_DIR}/sources/tests/Kirby2.dat"
         WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/tests")

foreach(fitter IN ITEMS ${nonlinear_fitters})
    alta_test(NAME "data2brdf_pinkfelt_${fitter}"
              COMMAND "data2brdf"
                      "--input"  "${CMAKE_SOURCE_DIR}/data/brdf/pink-felt-1d.alta"
                      "--output" "pink-felt-1d-${fitter}-blinn.func"
                      "--fitter" "nonlinear_fitter_${fitter}"
                      "--func"   "[nonlinear_function_diffuse, nonlinear_function_blinn]") 
endforeach()

if(PYTHONLIBS_FOUND AND PYBIND_FOUND)
    alta_test_python(NAME "python_test_arguments"
                     COMMAND "python3" "${CMAKE_SOURCE_DIR}/sources/tests/python/test-arguments.py")

    alta_test_python(NAME "python_test_vec"
                     COMMAND "python3" "${CMAKE_SOURCE_DIR}/sources/tests/python/test-vec.py")

    alta_test_python(NAME "python_test_function"
                     COMMAND "python3" "${CMAKE_SOURCE_DIR}/sources/tests/python/test-python-function.py")
endif()

# add a target to generate API documentation with Doxygen
if(DOXYGEN_FOUND)
	add_custom_target(doc
			  ${DOXYGEN_EXECUTABLE} ${CMAKE_SOURCE_DIR}/documents/doxygen.conf
			  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/documents/
			  COMMENT "Generating API documentation with Doxygen" VERBATIM)
endif(DOXYGEN_FOUND)
