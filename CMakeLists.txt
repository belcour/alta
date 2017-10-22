cmake_minimum_required (VERSION 3.2)
project (ALTA CXX)

# Find required packages
find_package(EIGEN3)

add_definitions(-DALTA_PLUGIN_DIRECTORY="/usr/lib/alta_plugins")

# Add core and library include files
include_directories ("external/build" "sources" ${EIGEN3_INCLUDE_DIR})

# Add ALTA core library
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

# Add main test suite
include_directories("sources/tests")

add_executable(test_data_io sources/tests/core/data-io.cpp)
target_link_libraries(test_data_io core)
target_compile_features(test_data_io PRIVATE cxx_range_for)

enable_testing()
add_test(data_io test_data_io)
