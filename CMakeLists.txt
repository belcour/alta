cmake_minimum_required (VERSION 3.2)
project (ALTA CXX)

find_package(EIGEN3)

add_definitions(-DALTA_PLUGIN_DIRECTORY="/usr/lib/alta_plugins")

include_directories ("external/build" "." ${EIGEN3_INCLUDE_DIR})

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

# Add dependencies

# Add main test suite
#add_executable (TestCovariance4D    tests/Covariance4D.cpp)
#add_executable (TestInvCovariance4D tests/InvCovariance4D.cpp)
#target_compile_features(TestCovariance4D    PRIVATE cxx_range_for)
#target_compile_features(TestInvCovariance4D PRIVATE cxx_range_for)

#enable_testing()
#add_test(TestCovariance4D    TestCovariance4D)
#add_test(TestInvCovariance4D TestInvCovariance4D)

#add_executable (Tutorial1 tutorials/tutorial1.cpp)
#target_compile_features(Tutorial1 PRIVATE cxx_range_for)
#add_executable (Tutorial2 tutorials/tutorial2.cpp)
#target_compile_features(Tutorial2 PRIVATE cxx_range_for)
