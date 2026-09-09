cmake_minimum_required(VERSION 3.25)

#cmake_policy(SET CMP0167 NEW)

project(Transfer_entropy CXX)

set(CMAKE_CXX_STANDARD 23)
cmake_policy(SET CMP0167 NEW)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++2b" COMPILER_SUPPORTS_CXX2B)

#set(Boost_USE_STATIC_LIBS   OFF)

find_package(GTest)
find_package (Eigen3 REQUIRED)
set(boost_system_DIR /usr/lib/cmake/)
find_package(Boost 1.90 COMPONENTS filesystem serialization log program_options REQUIRED)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
find_package(MPFR REQUIRED)
find_package(GMP REQUIRED)
find_package(absl)

include_directories(${EIGEN_INCLUDE_DIR})

add_subdirectory(third-party)
add_subdirectory(src_cpp)

include(CPack)
