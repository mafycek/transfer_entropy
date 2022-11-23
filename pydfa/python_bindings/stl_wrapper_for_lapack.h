/***************************************************************************
 *   Copyright (C) 2010 by Hynek Lavička   *
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

#ifndef STL_WRAPPER_FOR_LAPACK_H
#define STL_WRAPPER_FOR_LAPACK_H

#include <gsl/gsl_sf.h>

#include <vector>
#include <list>
#include "math_functions.h"

#include <boost/numeric/ublas/hermitian.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <lpp/lapack.hh>


namespace math_tools
{

double TrOperator ( boost::numeric::ublas::hermitian_adaptor<math_functions::complex_matrix, boost::numeric::ublas::upper> & input );

double TrOperator ( math_functions::complex_hermitian_matrix & input );

double TrOperator ( math_functions::complex_matrix & input );

int Eigenvalues ( math_functions::complex_hermitian_matrix & base , std::vector<double> & eigenvalues );

int Eigenvalues ( math_functions::complex_matrix & base , std::vector<std :: complex<double> > & eigenvalues );

int EigenvaluesAndEigenvectors ( math_functions::complex_hermitian_matrix & base , std::vector<double> & eigenvalues , std::vector<boost::numeric::ublas::vector<std :: complex<double> > * > & eigenvectors );

int EigenvaluesAndRightEigenvectors ( math_functions::complex_matrix & base , std::vector<std :: complex<double> > & eigenvalues , std::vector<math_functions::complex_column_vector * > & eigenvectors );

int EigenvaluesAndLeftEigenvectors ( math_functions::complex_matrix & base , std::vector<std :: complex<double> > & eigenvalues , std::vector<math_functions::complex_column_vector * > & eigenvectors );

int EigenvaluesAndLeftAndRightEigenvectors ( math_functions::complex_matrix & base , std::vector<std :: complex<double> > & eigenvalues , std::vector<math_functions::complex_column_vector * > & left_eigenvectors , std::vector<math_functions::complex_column_vector * > & right_eigenvectors );

void DestroyEigenvectors ( std::vector<boost::numeric::ublas::vector<std :: complex<double> > * > & eigenvectors );

int GeneralFunction ( math_functions::complex_hermitian_matrix & result , math_functions::complex_hermitian_matrix & base , double ( * function ) ( double ) );

int GeneralFunction ( math_functions::complex_matrix & result , math_functions::complex_hermitian_matrix & base , std :: complex<double> ( * function ) ( const std :: complex<double> & ) );

int GeneralFunction ( math_functions::complex_matrix & result , math_functions::complex_matrix & base , std :: complex<double> ( * function ) ( const std :: complex<double> & ) );

int Logarithm ( math_functions::complex_hermitian_matrix & result , math_functions::complex_hermitian_matrix & base );

template<class TYPE1, class TYPE2>
int Exp ( TYPE1 & result , TYPE2 & base )
{
    return GeneralFunction ( result , base , & std::exp );
}

int Sqrt ( math_functions::complex_hermitian_matrix & result , math_functions::complex_hermitian_matrix & base );

int Sqrt ( math_functions::complex_matrix & result , math_functions::complex_matrix & base );

int Sqr ( math_functions::complex_hermitian_matrix & result , math_functions::complex_hermitian_matrix & base );

int Sqr ( math_functions::complex_matrix & result , math_functions::complex_matrix & base );

int Inversion ( math_functions::complex_hermitian_matrix & result , math_functions::complex_hermitian_matrix & base );

int Inversion ( math_functions::complex_matrix & result , math_functions::complex_matrix & base );

int Inversion ( math_functions::real_matrix & result , math_functions::real_matrix & base );

template<class MATRIX_TYPE>
MATRIX_TYPE Inversion ( MATRIX_TYPE & mat )
{
	typedef typename MATRIX_TYPE::value_type TYPE_NAME;

	if ( mat.size1() == mat.size2 () )
	{
// inversion matrix construction
		MATRIX_TYPE inversion ( mat.size1() , mat.size2() );
		for ( unsigned int row = 0 ; row <  mat.size1() ; ++ row )
		{
			for ( unsigned int column = 0 ; column <  mat.size2() ; ++ column )
			{
				if ( row == column )
				{
					inversion ( row , column ) = static_cast <TYPE_NAME> ( 1.0 );
				}
				else
				{
					inversion ( row , column ) = static_cast <TYPE_NAME> ( 0.0 );
				}
			}
		}

// working copy of the matrix to invert
		MATRIX_TYPE working_copy ( mat );

// lower diagonal elimination
		for ( unsigned int pivot = 0 ; pivot < mat.size1 () ; ++ pivot )
		{
			unsigned int row = working_copy.size1 () - 1;
			if ( working_copy ( pivot , pivot ) == static_cast<TYPE_NAME> ( 0 ) )
			{
				// find a non-zero substitute for a pivot
				bool unable_to_find_subtitute = true;
				while ( row > pivot )
				{
					if ( working_copy ( row , pivot ) != static_cast<TYPE_NAME> ( 0 ) )
					{
						// exchange some row and pivot row
						for ( unsigned int column = pivot ; column < working_copy.size1 () ; ++ column )
						{
							TYPE_NAME swap = working_copy ( pivot , column );
							working_copy ( pivot , column ) = working_copy ( row , column );
							working_copy ( row , column ) = swap;
						}

						for ( unsigned int column = 0 ; column < working_copy.size1 () ; ++ column )
						{
							TYPE_NAME swap2 = inversion ( pivot , column );
							inversion ( pivot , column ) = inversion ( row , column );
							inversion ( row , column ) = swap2;
						}

						unable_to_find_subtitute = false;
						break;
					}
					-- row;
				}

				// singular matrix detected
				if ( unable_to_find_subtitute )
				{
					std::cerr << "Inversion: Singular matrix detected." << std::endl;
				}
			}

// setting pivot to 1
			for ( unsigned int column = pivot + 1 ; column < working_copy.size1 () ; ++ column )
			{
				working_copy ( pivot , column ) = working_copy ( pivot , column ) / working_copy ( pivot , pivot );
			}

			for ( unsigned int column = 0 ; column < working_copy.size1 () ; ++ column )
			{
				inversion ( pivot , column ) = inversion ( pivot , column ) / working_copy ( pivot , pivot );
			}
			working_copy ( pivot , pivot ) = static_cast<TYPE_NAME> ( 1 );

// eliminate non-zeros in pilot column
			for ( unsigned int zeroing_row = pivot + 1 ; zeroing_row <= row ; ++ zeroing_row )
			{
				for ( unsigned int column = pivot + 1 ; column <  working_copy.size1 () ; ++ column )
				{
					working_copy ( zeroing_row , column ) = working_copy ( zeroing_row , column ) * working_copy ( pivot , pivot ) - working_copy ( pivot , column ) * working_copy ( zeroing_row , pivot );
				}

// preforming transformation of the inversion of the matrix
				for ( unsigned int column = 0 ; column < working_copy.size1 () ; ++ column )
				{
					inversion ( zeroing_row , column ) = inversion ( zeroing_row , column ) * working_copy ( pivot , pivot ) - inversion ( pivot , column ) * working_copy ( zeroing_row , pivot );
				}

				working_copy ( zeroing_row , pivot ) = working_copy ( zeroing_row , pivot ) * working_copy ( pivot , pivot ) - working_copy ( zeroing_row , pivot ) * mat ( pivot , pivot );
			}
		}

		for ( unsigned int pivot = mat.size1 () - 1 ; pivot > 0 ; -- pivot )
		{
// eliminate non-zeros in pilot column
			for ( unsigned int zeroing_row = 0 ; zeroing_row < pivot ; ++ zeroing_row )
			{
// preforming transformation of the inversion of the matrix
				for ( unsigned int column = 0 ; column < working_copy.size1 () ; ++ column )
				{
					inversion ( zeroing_row , column ) = inversion ( zeroing_row , column ) * working_copy ( pivot , pivot ) - inversion ( pivot , column ) * working_copy ( zeroing_row , pivot );
				}

				working_copy ( zeroing_row , pivot ) = static_cast<TYPE_NAME> ( 0 ) ;
			}
		}

		return inversion;
	}
	else
	{
		std::cout << "Inversion: It is not a square matrix." << std::endl;
		return mat;
	}
}

double QuantumEntropy ( math_functions::complex_hermitian_matrix & input );

void Join ( math_functions::complex_matrix & result , math_functions::complex_matrix & base );

void HamiltonianContinuousQuantumRandomWalk ( math_functions::complex_hermitian_matrix & result  , double phase = 0 );

}

#endif
