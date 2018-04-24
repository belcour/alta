cmake_minimum_required (VERSION 3.2)
project (alta CXX)

# Preconfiguration
# Set Python version number
set(Python_ADDITIONAL_VERSIONS 2.7)

# Find required packages
find_package(Eigen3)
find_package(OpenMP)
find_package(PythonLibs)
find_package(Boost COMPONENTS python)

add_definitions(-DALTA_PLUGIN_DIRECTORY="/usr/lib/alta_plugins")

# Add core and library include files
include_directories("external/build" "sources" ${EIGEN3_INCLUDE_DIR})
include_directories("external" "external/Catch/include")

# Look for header only dependencies
find_file(EIGEN_FOUND Eigen     HINTS ${EIGEN3_INCLUDE_DIR})
find_file(CATCH_FOUND catch.hpp HINTS external/Catch/include)

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
	target_compile_features(${name} PRIVATE cxx_range_for)
endfunction(alta_add_plugin)

function(alta_add_soft name sources)
	add_executable(${name} sources/softs/${sources})
	target_link_libraries(${name} core)
	target_compile_features(${name} PRIVATE cxx_range_for)
endfunction(alta_add_soft)

function(alta_test_unit name sources)
        add_executable(${name} sources/tests/${sources})
        target_link_libraries(${name} core)
        target_compile_features(${name} PRIVATE cxx_range_for)
        add_test(NAME ${name} COMMAND ${CMAKE_BINARY_DIR}/tests/${name})
        set_tests_properties(${name} PROPERTIES ENVIRONMENT "TEST_DATA_DIRECTORY=${CMAKE_SOURCE_DIR}/sources/tests")
endfunction(alta_test_unit)

function(alta_test_soft name)
        add_test(NAME ${name} COMMAND ${CMAKE_BINARY_DIR}/softs/${name} --help)
endfunction(alta_test_soft)


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
target_compile_features(core PRIVATE cxx_range_for)

target_link_libraries(core ${CMAKE_DL_LIBS})


#############################
#         Plugins           #
#############################

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/plugins)
set(CMAKE_SHARED_LIBRARY_PREFIX "")

# Data
alta_add_plugin(data_merl data_io/merl.cpp)
alta_add_plugin(data_utia data_io/utia.cpp)

# Functions
alta_add_plugin(nonlinear_function_diffuse nonlinear_function_diffuse/function.cpp)

# Fitters
alta_add_plugin(nonlinear_fitter_eigen     nonlinear_fitter_eigen/fitter.cpp)

# Python bindings
if(PYTHONLIBS_FOUND AND Boost_FOUND)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python)
    include_directories(${PYTHON_INCLUDE_DIRS} ${Boost_INCLUDE_DIR})
    add_library(alta SHARED sources/python/alta.cpp)
    target_link_libraries(alta core ${PYTHON_LIBRARIES} ${Boost_LIBRARIES})
    target_compile_features(alta PRIVATE cxx_range_for)
    if(APPLE)
        set_target_properties(alta PROPERTIES SUFFIX ".so")
    endif()
endif()

#############################
#          Softs            #
#############################

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/softs)

alta_add_soft(brdf2data  brdf2data/main.cpp)
alta_add_soft(data2data  data2data/main.cpp)
alta_add_soft(data2brdf  data2brdf/main.cpp)
alta_add_soft(data2stats data2stats/data2stats.cpp)

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

if(CATCH_FOUND)
    alta_test_unit(conversion-catch-1 core/conversion-catch-1.cpp )
endif()
