cmake_minimum_required(VERSION 3.7)
project(Pydfa)

INCLUDE(CheckFunctionExists)
FIND_PACKAGE(PythonLibs REQUIRED)
FIND_PACKAGE(Boost COMPONENTS python system thread REQUIRED)
FIND_PACKAGE(GSL REQUIRED)
FIND_PACKAGE(PkgConfig REQUIRED)
PKG_CHECK_MODULES(MPFR REQUIRED mpfr)
PKG_CHECK_MODULES(LAPACK REQUIRED lapack) #-netlib

check_library_exists(m pow "" LIBM)

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR} "/usr/include/python3.11")

link_directories( "${PYTHON_LIBRARIES}" )

# detrended_fluctuation_analysis.cpp
ADD_LIBRARY ( pydfa SHARED pydfa.cpp regression.cpp  fluctuation_analyzer.cpp stl_wrapper_for_lapack.cpp math_functions.cpp debug_support.cpp auxiliary_functions.cpp )
# math_functions.cpp

TARGET_LINK_LIBRARIES(pydfa LINK_PUBLIC ${PYTHON_LIBRARIES} Boost::python Boost::system Boost::thread ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${MPFR_LIBRARIES} ${LAPACK_LIBRARIES} m)
# -lboost_python-py37 -lboost_system -lpython3.6m   -llapack -lblas -lgfortran 

set_target_properties ( pydfa PROPERTIES PREFIX "" )
