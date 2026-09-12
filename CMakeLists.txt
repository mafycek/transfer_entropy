cmake_minimum_required(VERSION 3.25)

cmake_policy(SET CMP0175 NEW)

project(Transfer_entropy CXX)

set(CMAKE_CXX_STANDARD 26)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
list(APPEND PKG_CONFIG_PATH /usr/lib/x86_64-linux-gnu/pkgconfig )

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++2b" COMPILER_SUPPORTS_CXX2B)

set(Boost_USE_STATIC_LIBS   ON)

find_package(ICU COMPONENTS data i18n uc REQUIRED)
find_package(ZLIB REQUIRED)
find_package(zstd REQUIRED)
find_package(GTest)
find_package (Eigen3 REQUIRED)
find_package(Boost 1.70 COMPONENTS filesystem serialization log program_options iostreams REQUIRED) #system
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
find_package(MPFR REQUIRED)
find_package(GMP REQUIRED)
find_package(TBB REQUIRED)
find_package(absl REQUIRED)
find_package(Python3 3.14.0 COMPONENTS Interpreter Development NumPy REQUIRED)

MESSAGE(STATUS ${Boost_LIBRARIES} YES)

# mffm fftw3
include_directories("/usr/include/mffm")
link_libraries(2Dfft fft)

include_directories(${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS})
include_directories(${EIGEN_INCLUDE_DIR})
link_libraries(${Python3_LIBRARIES})
link_libraries(TBB::tbb)

add_subdirectory(third-party)
add_subdirectory(src_cpp)

include(CPack)
