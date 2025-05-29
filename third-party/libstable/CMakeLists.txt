cmake_minimum_required(VERSION 3.20)

#cmake_policy(SET CMP0167 NEW)

project(libstable LANGUAGES C VERSION 1.0.0)
set(CMAKE_CXX_STANDARD 23)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++2b" COMPILER_SUPPORTS_CXX2B)

find_package(GSL REQUIRED)

add_subdirectory(stable/src)
