/*! \file math_functions.h
\brief Mathematical functions declaration file

Mathematical functions that are used in other projects.
\ingroup zarja
*/
/***************************************************************************
 *   Copyright (C) 2005 by Hynek Lavička   *
 *   h.lavicka@email.cz   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <boost/numeric/ublas/io.hpp>

#include <boost/numeric/ublas/hermitian.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

#include <boost/math/special_functions.hpp>

#ifndef MATH_FUNCTIONS_H
#define MATH_FUNCTIONS_H

/*
// Matrices

#include "diagonal_matrix.h"
#include "hermitian_matrix.h"
#include "symmetric_matrix.h"
#include "matrix.h"

// Vectors

#include "general_vector.h"
#include "column_vector.h"
#include "row_vector.h"
*/

/**
\brief Namespace encapulates all mathematical functions.
\ingroup zarja
*/
namespace math_functions
{

/// A function which takes double and returns double
/**
The function computes absolut value of variable
*/
double absolut_value ( double a );

int absolut_value ( int a );

/// A function which takes integer and returns double
/**
The function computes factorial
\param i parameter of factorial
\return factorial
*/
double factorial ( int i );

double pythagoras_hypotenuse ( double a , double b );

double conditional_sign ( double a , double b );

double & INDEX_MATRIX ( double * A , int i , int j , int N );

std :: complex<double> & INDEX_MATRIX ( std :: complex<double> * A , int i , int j , int N );

int Jakobi_eigenvalue_method_with_eigenvectors ( const int N , double * A , double * eigenvalue , double * eigenvectors , int * number_of_iterations );

int Jakobi_eigenvalue_method ( const int N , double * A , double * eigenvalue , int * number_of_iterations );

int QL_algorithm_and_household_reduction_with_eigenvectors ( int N , double * A , double * eigenvalues , double * eigenvectors  );

int RealMatrixMultiplication ( double * A , int columns_A , int rows_A , double * B , int columns_B , int rows_B , double * C , int columns_C , int rows_C );

int ComplexMatrixMultiplication ( std :: complex<double> * A , int columns_A , int rows_A , std :: complex<double> * B , int columns_B , int rows_B , std :: complex<double> * C , int columns_C , int rows_C );

int ComplexMatrixTransposition ( std :: complex<double> * A , int columns_A , int rows_A , std :: complex<double> * B , int columns_B , int rows_B );

int ComplexMatrixConjugation ( std :: complex<double> * A , int columns_A , int rows_A , std :: complex<double> * B , int columns_B , int rows_B );

int ComplexMatrixHermitianConjugation ( std :: complex<double> * A , int columns_A , int rows_A , std :: complex<double> * B , int columns_B , int rows_B );

int ComplexErrorPrintMatrix ( std :: complex<double> * A , int columns_A , int rows_A );

long GetDigit ( long number , long decimals , long rank );

double sqr ( double x );

double cube ( double x );

double sqrt_inverse ( double x );

int minimum ( int a , int b );

double Signum ( double x );

double HeavysideFunction ( double x );

double Factorial ( int i );

double LogarithmusFactorial ( int i );

double Combination ( int n , int m );

double LogarithmusCombination ( int n , int m );

double GeneralizedHarmonicNumber ( unsigned int N , double s );

double & ScalarMultiplication ( boost::numeric::ublas::vector<double> & vector1 , boost::numeric::ublas::vector<double> & vector2 );

std :: complex<double> & ScalarMultiplication ( boost::numeric::ublas::vector<std :: complex<double> > & vector1 , boost::numeric::ublas::vector<std :: complex<double> > & vector2 );

typedef boost::numeric::ublas::hermitian_matrix<std :: complex<double>, boost::numeric::ublas::upper> complex_hermitian_matrix;
typedef boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::upper> real_symmetric_matrix;
typedef boost::numeric::ublas::matrix<double> real_matrix;
typedef boost::numeric::ublas::matrix<std :: complex<double> > complex_matrix;
typedef boost::numeric::ublas::matrix<double> real_sparse_matrix;
typedef boost::numeric::ublas::mapped_matrix<std :: complex<double> > complex_sparse_matrix;

typedef boost::numeric::ublas::vector<double> real_column_vector;
typedef boost::numeric::ublas::vector<std :: complex<double> > complex_column_vector;


namespace multiprecision
{

template <class TYPE>
inline TYPE power ( TYPE a , TYPE b )
{
	return boost::multiprecision::pow ( a , b );
}

template<>
inline double power ( double a , double b )
{
	return pow ( a , b );
}

template<>
inline float power ( float a , float b )
{
	return pow ( a , b );
}

template <class TYPE>
inline TYPE logarithm ( TYPE a )
{
	return boost::math::log1p  ( a + static_cast<TYPE> ( 1 ) );
}

template<>
inline double logarithm ( double a )
{
	return log ( a  );
}

template<>
inline float logarithm ( float a )
{
	return log ( a  );
}

template <class TYPE>
inline TYPE exponential ( TYPE a )
{
	return boost::math::expm1 ( a ) + static_cast<TYPE> ( 1 );
}

template<>
inline double exponential ( double a )
{
	return exp ( a  );
}

template<>
inline float exponential ( float a )
{
	return exp ( a  );
}

template <class TYPE>
inline TYPE absolute_value ( TYPE a )
{
	return ( a >= 0 ? a : -a );
}

template<>
inline double absolute_value ( double a )
{
	return std::abs ( a );
}

template<>
inline float absolute_value ( float a )
{
	return std::abs ( a );
}

template <class TYPE>
inline bool isnan ( TYPE a )
{
	return boost::math::isnan (a);
}

template<>
inline bool isnan ( double a )
{
	return std::isnan ( a );
}

template<>
inline bool isnan ( float a )
{
	return std::isnan ( a );
}


}

}

#endif
