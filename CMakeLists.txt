cmake_minimum_required (VERSION 3.2)
project (alta CXX)

# Find required packages
find_package(Eigen3 REQUIRED)
find_package(OpenMP)

add_definitions(-DALTA_PLUGIN_DIRECTORY="/usr/lib/alta_plugins")

# Add core and library include files
include_directories("external/build" "sources" ${EIGEN3_INCLUDE_DIR})
include_directories("external")

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


#############################
#         ALTA core         #
#############################

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

add_library(core STATIC sources/core/common.cpp
			sources/core/metrics.cpp
			sources/core/params.cpp
			sources/core/data.cpp
			sources/core/data_storage.cpp
			sources/core/function.cpp
			sources/core/plugins_manager.cpp
			sources/core/vertical_segment.cpp
			sources/core/rational_function.cpp)
target_compile_features(core PRIVATE cxx_range_for)

target_link_libraries(core ${CMAKE_DL_LIBS})


#############################
#         Plugins           #
#############################

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/plugins)
set(CMAKE_SHARED_MODULE_PREFIX "")

# Data
alta_add_plugin(data_merl data_io/merl.cpp)
alta_add_plugin(data_utia data_io/utia.cpp)

# Functions
alta_add_plugin(nonlinear_function_diffuse nonlinear_function_diffuse/function.cpp)

# Fitters 
alta_add_plugin(nonlinear_fitter_eigen     nonlinear_fitter_eigen/fitter.cpp)

#############################
#          Softs            #
#############################

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/softs)

add_executable(brdf2data sources/softs/brdf2data/main.cpp)
target_link_libraries(brdf2data core)
target_compile_features(brdf2data PRIVATE cxx_range_for)

alta_add_soft(data2data data2data/main.cpp)
alta_add_soft(data2brdf data2brdf/main.cpp)


#############################
#          Tests            #
#############################

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests)

include_directories("sources/tests")

add_executable(test_data_io sources/tests/core/data-io.cpp)
target_link_libraries(test_data_io core)
target_compile_features(test_data_io PRIVATE cxx_range_for)

enable_testing()
add_test(data_io test_data_io)
set_tests_properties(data_io PROPERTIES ENVIRONMENT "TEST_DATA_DIRECTORY=${CMAKE_SOURCE_DIR}/sources/tests")
