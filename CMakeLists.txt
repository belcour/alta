cmake_minimum_required (VERSION 3.2)
project (alta CXX)

# Find required packages
find_package(Eigen3)
find_package(OpenMP)

add_definitions(-DALTA_PLUGIN_DIRECTORY="/usr/lib/alta_plugins")

# Add core and library include files
include_directories("external/build" "sources" ${EIGEN3_INCLUDE_DIR})
include_directories("external" "external/Catch/include")

# TODO Look for Eigen library there (do not rely on find_package)

# Update compilation option to include OpenMP if present
if(OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS ${OpenMP_CXX_FLAGS})
endif(OPENMP_FOUND)

# TODO Look for catch library
set(CATCH_FOUND TRUE)

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

function(alta_add_test name sources)
        add_executable(${name} sources/tests/${sources})
        target_link_libraries(${name} core)
        target_compile_features(${name} PRIVATE cxx_range_for)
        add_test(NAME ${name} COMMAND ${CMAKE_BINARY_DIR}/tests/${name})
        set_tests_properties(${name} PROPERTIES ENVIRONMENT "TEST_DATA_DIRECTORY=${CMAKE_SOURCE_DIR}/sources/tests")
endfunction(alta_add_test)


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

alta_add_soft(brdf2data  brdf2data/main.cpp)
alta_add_soft(data2data  data2data/main.cpp)
alta_add_soft(data2brdf  data2brdf/main.cpp)
alta_add_soft(data2stats data2stats/data2stats.cpp)

#############################
#          Tests            #
#############################

enable_testing()
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests)

include_directories("sources/tests")

alta_add_test(test_data_io  core/data-io.cpp)
alta_add_test(half-test-1   core/half-test-1.cpp)
alta_add_test(half-test-2   core/half-test-2.cpp)
alta_add_test(half-test-3   core/half-test-3.cpp)
alta_add_test(half-test-4   core/half-test-4.cpp)
alta_add_test(nonlinear-fit core/nonlinear-fit.cpp)
alta_add_test(params-test-1 core/params-test-1.cpp)
alta_add_test(params-test-2 core/params-test-2.cpp)

if(CPPQUICKCHECK_FOUND)
    alta_add_test(params-qc-1 core/params-qc-1.cpp)
endif()

if(CATCH_FOUND)
    alta_add_test(conversion-catch-1 core/conversion-catch-1.cpp )
endif()
