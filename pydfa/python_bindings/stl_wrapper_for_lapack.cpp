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

#include "stl_wrapper_for_lapack.h"

namespace math_tools
{

double TrOperator ( boost::numeric::ublas::hermitian_adaptor<math_functions::complex_matrix, boost::numeric::ublas::upper> & input )
{
    if ( input . size1 ( ) == input . size2 ( ) )
    {
        double sum = 0;
        for ( unsigned int i = 0 ; i < input . size1 ( ) ; i ++ )
        {
            std::complex <double> diagonal_element = input ( i , i );
            sum += diagonal_element . real ( );
        }

        return sum;
    }
    else
    {
// 		Dout ( dc::warning, "TrOperator: Input is not a square matrix." );
// 		cout << "Problem" << std :: endl;

        return -1;
    }
}


double TrOperator ( math_functions::complex_hermitian_matrix & input )
{
    if ( input . size1 ( ) == input . size2 ( ) )
    {
        double sum = 0;
        for ( unsigned int i = 0 ; i < input . size1 ( ) ; i ++ )
        {
            std::complex <double> diagonal_element = input ( i , i );
            sum += diagonal_element . real ( );
        }

        return sum;
    }
    else
    {
// 		Dout ( dc::warning, "TrOperator: Input is not a square matrix." );
// 		cout << "Problem" << std :: endl;

        return -1;
    }
}


double TrOperator ( math_functions::complex_matrix & input )
{
    if ( input . size1 ( ) == input . size2 ( ) )
    {
        double sum = 0;
        for ( unsigned int i = 0 ; i < input . size1 ( ) ; i ++ )
        {
            std::complex <double> diagonal_element = input ( i , i );
            sum += diagonal_element . real ( );
        }

        return sum;
    }
    else
    {
// 		Dout ( dc::warning, "TrOperator: Input is not a square matrix." );
// 		cout << "Problem" << std :: endl;

        return -1;
    }
}


int Eigenvalues ( math_functions::complex_hermitian_matrix & base , std::vector<double> & eigenvalues )
{
    eigenvalues . clear ( );

    long info;
    char jobs = 'N';
    char uplo = 'U';
    long rows = base . size1 ( );
    long columns = base . size2 ( );
    boost::numeric::ublas::zero_matrix<std :: complex<double> > zero ( rows , columns );

    double * eigenvalues_double = new double [ columns ];
    std :: complex<double> * temporary_matrix = new std :: complex<double> [ rows * columns ];

    for ( int i = 0 ; i < columns ; i ++ )
    {
        temporary_matrix [ i * rows + i ] = base ( i , i );

        for ( int j = i + 1 ; j < rows ; j ++ )
        {
            std :: complex<double> element = base ( j , i );

            temporary_matrix [ i * columns + j ] = element;
            temporary_matrix [ j * rows + i ] = conj ( element );
        }
    }

    lpp :: heev ( & jobs , & uplo , & rows , temporary_matrix , & columns , eigenvalues_double , & info );

// 	for ( int i = 0 ; i < columns ; i ++ )
// 	{
// 		eigenvalues . push_back ( eigenvalues_double [ i ] );
// 	}

    if ( info != 0 )
    {
        if ( info < 0 )
        {
            std::cerr << "Heev got wrong argument at position " << - info << std::endl;
        }
        else
        {
            std::cerr << "Heev did not converge for " << info << " dimension and no eigenvectors produced." << std::endl;

            for ( int i = 0 ; i < columns - info ; i ++ )
            {
                eigenvalues . push_back ( eigenvalues_double [ i ] );

#ifdef DEBUG
                std::cerr << eigenvalues_double [ i ] << " ";
#endif
            }
#ifdef DEBUG
            std::cerr << std::endl;
#endif
        }

        std::cerr << "Heev finished with value: " << info << std::endl;
    }
    else
    {
        for ( int i = 0 ; i < columns ; i ++ )
        {
            eigenvalues . push_back ( eigenvalues_double [ i ] );
        }
    }

    delete [] eigenvalues_double;
    delete [] temporary_matrix;

    return info;
}


int Eigenvalues ( math_functions::complex_matrix & base , std::vector<std :: complex<double> > & eigenvalues )
{
    long info = 0;
    char jobs = 'V';
    char jobsN = 'N';
    const long rows = base . size1 ( );
    const long columns = base . size2 ( );
    boost::numeric::ublas::zero_matrix<std :: complex<double> > zero ( rows , columns );

    std :: complex<double> * eigenvalues_complex = new std :: complex<double> [ columns ];
    std :: complex<double> * temporary_matrix = new std :: complex<double> [ rows * columns ];
// 	std :: complex<double> * eigenvectors = new std :: complex<double> [ 1 * 1 ];

#ifdef DEBUG
    std::cout << "Temporary matrix:" << std::endl;
#endif

    for ( int i = 0 ; i < columns ; i ++ )
    {
        for ( int j = 0 ; j < rows ; j ++ )
        {
            std :: complex<double> element = base ( j , i );

            temporary_matrix [ i * columns + j ] = element;

#ifdef DEBUG
            std::cout << element << " ";
#endif
        }
        eigenvalues_complex [ i ] = 0;

#ifdef DEBUG
        std::cout << std::endl;
#endif
    }

#ifdef DEBUG
    std::cout << "Calculating eigenvalues for matrix:" << std::endl;
    std::cout << base << std::endl;
    std::cout << "Eigenvalues: ";
#endif

    lpp :: geev ( & jobsN , & jobsN , & rows , temporary_matrix , & columns , eigenvalues_complex , nullptr , & rows , nullptr , & columns , & info );

    if ( info != 0 )
    {
        if ( info < 0 )
        {
            std::cerr << "Geev got wrong argument at position " << - info << std::endl;
        }
        else
        {
            std::cerr << "Geev did not converge for " << info << " dimension and no eigenvectors produced." << std::endl;

            for ( int i = 0 ; i < columns - info ; i ++ )
            {
                eigenvalues . push_back ( eigenvalues_complex [ i ] );

#ifdef DEBUG
                std::cerr << eigenvalues_complex [ i ] << " ";
#endif
            }
#ifdef DEBUG
            std::cerr << std::endl;
#endif
        }

        std::cerr << "Geev finished with value: " << info << std::endl;
    }
    else
    {
        for ( int i = 0 ; i < columns ; i ++ )
        {
            eigenvalues . push_back ( eigenvalues_complex [ i ] );
        }
    }

    delete [] eigenvalues_complex;
    delete [] temporary_matrix;
// 	delete eigenvectors;

    return info;
}


int EigenvaluesAndEigenvectors ( math_functions::complex_hermitian_matrix & base , std::vector<double> & eigenvalues , std::vector<math_functions::complex_column_vector * > & eigenvectors )
{
    DestroyEigenvectors ( eigenvectors );

    long info;
    char jobs = 'V';
    char uplo = 'U';
    long rows = base . size1 ( );
    long columns = base . size2 ( );
    boost::numeric::ublas::zero_matrix<std :: complex<double> > zero ( rows , columns );

    double * eigenvalues_double = new double [ columns ];
    std :: complex<double> * temporary_matrix = new std :: complex<double> [ rows * columns ];

#ifdef DEBUG
    std::cout << "Temporary matrix:" << std::endl;
#endif

    for ( int i = 0 ; i < columns ; i ++ )
    {
        temporary_matrix [ i * rows + i ] = base ( i , i );

        for ( int j = i + 1 ; j < rows ; j ++ )
        {
            std :: complex<double> element = base ( j , i );

            temporary_matrix [ i * columns + j ] = element;
            temporary_matrix [ j * rows + i ] = conj ( element );

#ifdef DEBUG
            std::cout << element << " ";
#endif
        }
        eigenvalues_double [ i ] = 0;

#ifdef DEBUG
        std::cout << std::endl;
#endif
    }

#ifdef DEBUG
    std::cout << "Calculating eigenvalues for matrix:" << std::endl;
    std::cout << base << std::endl;
    std::cout << "Eigenvalues: ";
#endif

    lpp :: heev ( & jobs , & uplo , & rows , temporary_matrix , & columns , eigenvalues_double , & info );

    if ( info == 0 )
    {
        for ( int i = 0 ; i < columns ; i ++ )
        {
            eigenvalues . push_back ( eigenvalues_double [ i ] );

#ifdef DEBUG
            std::cout << eigenvalues_double [ i ] << " ";
#endif
        }

#ifdef DEBUG
        std::cout << std::endl;
#endif

#ifdef DEBUG
        std::cout << "Eigenvectors: " << std::endl;
#endif

        for ( int i = 0 ; i < rows ; i ++ )
        {
            boost::numeric::ublas::vector<std :: complex<double> > * eigenvector = new boost::numeric::ublas::vector<std :: complex<double> > ( columns );
            eigenvectors.push_back ( eigenvector );
            for ( int j = 0 ; j < columns ; j ++ )
            {
                std :: complex<double> value = temporary_matrix [ i * rows + j ];
                ( * eigenvector ) [ j ]  = value;
            }
#ifdef DEBUG
            std::cout << ( * eigenvector ) <<std::endl;
#endif
        }

#ifdef DEBUG
        std::cout << std::endl;
#endif
        return info;
    }
    else
    {
        std::cout << "ZHEEV finished with value: " << info << std::endl;
        return info;
    }
}


int EigenvaluesAndRightEigenvectors ( math_functions::complex_matrix & base , std::vector<std :: complex<double> > & eigenvalues , std::vector<math_functions::complex_column_vector * > & eigenvectors )
{
    DestroyEigenvectors ( eigenvectors );

    long info;
    char jobs = 'V';
    char jobsN = 'N';
    const long rows = base . size1 ( );
    const long columns = base . size2 ( );
    boost::numeric::ublas::zero_matrix<std :: complex<double> > zero ( rows , columns );

    std :: complex<double> * eigenvalues_complex = new std :: complex<double> [ columns ];
    std :: complex<double> * temporary_matrix = new std :: complex<double> [ rows * columns ];
    std :: complex<double> * right_eigenvectors = new std :: complex<double> [ rows * columns ];

#ifdef DEBUG
    std::cout << "Temporary matrix:" << std::endl;
#endif

    for ( int i = 0 ; i < columns ; i ++ )
    {
        for ( int j = 0 ; j < rows ; j ++ )
        {
            std :: complex<double> element = base ( j , i );

            temporary_matrix [ i * columns + j ] = element;
            right_eigenvectors [ i * columns + j ] = 0;

#ifdef DEBUG
            std::cout << element << " ";
#endif
        }
        eigenvalues_complex [ i ] = 0;

#ifdef DEBUG
        std::cout << std::endl;
#endif
    }

#ifdef DEBUG
    std::cout << "Calculating eigenvalues for matrix:" << std::endl;
    std::cout << base << std::endl;
    std::cout << "Eigenvalues: ";
#endif

    lpp :: geev ( & jobsN , & jobs , & rows , temporary_matrix , & columns , eigenvalues_complex , right_eigenvectors , & rows , right_eigenvectors , & columns , & info );

    if ( info != 0 )
    {
        for ( int i = 0 ; i < columns ; i ++ )
        {
            eigenvalues . push_back ( eigenvalues_complex [ i ] );

#ifdef DEBUG
            std::cout << eigenvalues_complex [ i ] << " ";
#endif
        }

#ifdef DEBUG
        std::cout << std::endl;
#endif

#ifdef DEBUG
        std::cout << "Right eigenvectors: " << std::endl;
#endif

        for ( int i = 0 ; i < rows ; i ++ )
        {
            boost::numeric::ublas::vector<std :: complex<double> > * eigenvector = new boost::numeric::ublas::vector<std :: complex<double> > ( columns );
            eigenvectors.push_back ( eigenvector );
            for ( int j = 0 ; j < columns ; j ++ )
            {
                std :: complex<double> value = right_eigenvectors [ i * rows + j ];
                ( * eigenvector ) [ j ]  = value;
            }
#ifdef DEBUG
            std::cout << ( * eigenvector ) <<std::endl;
#endif
        }

#ifdef DEBUG
        std::cout << std::endl;
#endif

        delete [] eigenvalues_complex;
        delete [] right_eigenvectors;
        delete [] temporary_matrix;

        return info;
    }
    else
    {
        std::cout << "ZGEEV finished with value: " << info << std::endl;
        return info;
    }
}


int EigenvaluesAndLeftEigenvectors ( math_functions::complex_matrix & base , std::vector<std :: complex<double> > & eigenvalues , std::vector<math_functions::complex_column_vector * > & eigenvectors )
{
    DestroyEigenvectors ( eigenvectors );

    long info;
    char jobs = 'V';
    char jobsN = 'N';
    const long rows = base . size1 ( );
    const long columns = base . size2 ( );
    boost::numeric::ublas::zero_matrix<std :: complex<double> > zero ( rows , columns );

    std :: complex<double> * eigenvalues_complex = new std :: complex<double> [ columns ];
    std :: complex<double> * temporary_matrix = new std :: complex<double> [ rows * columns ];
    std :: complex<double> * left_eigenvectors = new std :: complex<double> [ rows * columns ];

#ifdef DEBUG
    std::cout << "Temporary matrix:" << std::endl;
#endif

    for ( int i = 0 ; i < columns ; i ++ )
    {
        for ( int j = 0 ; j < rows ; j ++ )
        {
            std :: complex<double> element = base ( j , i );

            temporary_matrix [ i * columns + j ] = element;
            left_eigenvectors [ i * columns + j ] = 0;

#ifdef DEBUG
            std::cout << element << " ";
#endif
        }
        eigenvalues_complex [ i ] = 0;

#ifdef DEBUG
        std::cout << std::endl;
#endif
    }

#ifdef DEBUG
    std::cout << "Calculating eigenvalues for matrix:" << std::endl;
    std::cout << base << std::endl;
    std::cout << "Eigenvalues: ";
#endif

    lpp :: geev ( & jobsN , & jobs , & rows , temporary_matrix , & columns , eigenvalues_complex , left_eigenvectors , & rows , left_eigenvectors , & columns , & info );

    if ( info != 0 )
    {
        for ( int i = 0 ; i < columns ; i ++ )
        {
            eigenvalues . push_back ( eigenvalues_complex [ i ] );

#ifdef DEBUG
            std::cout << eigenvalues_complex [ i ] << " ";
#endif
        }

#ifdef DEBUG
        std::cout << std::endl;
#endif

#ifdef DEBUG
        std::cout << "Right eigenvectors: " << std::endl;
#endif

        for ( int i = 0 ; i < rows ; i ++ )
        {
            boost::numeric::ublas::vector<std :: complex<double> > * eigenvector = new boost::numeric::ublas::vector<std :: complex<double> > ( columns );
            eigenvectors.push_back ( eigenvector );
            for ( int j = 0 ; j < columns ; j ++ )
            {
                std :: complex<double> value = left_eigenvectors [ i * rows + j ];
                ( * eigenvector ) [ j ]  = value;
            }
#ifdef DEBUG
            std::cout << ( * eigenvector ) <<std::endl;
#endif
        }

#ifdef DEBUG
        std::cout << std::endl;
#endif

        delete [] eigenvalues_complex;
        delete [] left_eigenvectors;
        delete [] temporary_matrix;

        return info;
    }
    else
    {
        std::cout << "ZGEEV finished with value: " << info << std::endl;
        return info;
    }
}


int EigenvaluesAndLeftAndRightEigenvectors ( math_functions::complex_matrix & base , std::vector<std :: complex<double> > & eigenvalues , std::vector<math_functions::complex_column_vector * > & left_eigenvectors , std::vector<math_functions::complex_column_vector * > & right_eigenvectors )
{
    DestroyEigenvectors ( left_eigenvectors );
    DestroyEigenvectors ( right_eigenvectors );

    long info;
    char jobs = 'V';
    const long rows = base . size1 ( );
    const long columns = base . size2 ( );
    boost::numeric::ublas::zero_matrix<std :: complex<double> > zero ( rows , columns );

    std :: complex<double> * eigenvalues_complex = new std :: complex<double> [ columns ];
    std :: complex<double> * temporary_matrix = new std :: complex<double> [ rows * columns ];
    std :: complex<double> * geev_left_eigenvectors = new std :: complex<double> [ rows * columns ];
    std :: complex<double> * geev_right_eigenvectors = new std :: complex<double> [ rows * columns ];

#ifdef DEBUG
    std::cout << "Temporary matrix:" << std::endl;
#endif

    for ( int i = 0 ; i < columns ; i ++ )
    {
        for ( int j = 0 ; j < rows ; j ++ )
        {
            std :: complex<double> element = base ( j , i );

            temporary_matrix [ i * columns + j ] = element;
            geev_left_eigenvectors [ i * columns + j ] = 0;
            geev_right_eigenvectors [ i * columns + j ] = 0;

#ifdef DEBUG
            std::cout << element << " ";
#endif
        }
        eigenvalues_complex [ i ] = 0;

#ifdef DEBUG
        std::cout << std::endl;
#endif
    }

#ifdef DEBUG
    std::cout << "Calculating eigenvalues for matrix:" << std::endl;
    std::cout << base << std::endl;
    std::cout << "Eigenvalues: ";
#endif

    lpp :: geev ( & jobs , & jobs , & rows , temporary_matrix , & columns , eigenvalues_complex , geev_left_eigenvectors , & rows , geev_right_eigenvectors , & columns , & info );

    if ( info != 0 )
    {
        for ( int i = 0 ; i < columns ; i ++ )
        {
            eigenvalues . push_back ( eigenvalues_complex [ i ] );

#ifdef DEBUG
            std::cout << eigenvalues_complex [ i ] << " ";
#endif
        }

#ifdef DEBUG
        std::cout << std::endl;
#endif

#ifdef DEBUG
        std::cout << "Right eigenvectors: " << std::endl;
#endif

        for ( int i = 0 ; i < rows ; i ++ )
        {
            math_functions::complex_column_vector * left_eigenvector = new math_functions::complex_column_vector ( columns );
            math_functions::complex_column_vector * right_eigenvector = new math_functions::complex_column_vector ( columns );
            left_eigenvectors.push_back ( left_eigenvector );
            right_eigenvectors.push_back( right_eigenvector);

            for ( int j = 0 ; j < columns ; j ++ )
            {
                std :: complex<double> value = geev_left_eigenvectors [ i * rows + j ];
                ( * left_eigenvector ) [ j ]  = value;

                std :: complex<double> value2 = geev_right_eigenvectors [ i * rows + j ];
                ( * right_eigenvector ) [ j ]  = value2;
            }
#ifdef DEBUG
            std::cout << ( * left_eigenvector ) <<std::endl;
            std::cout << ( * right_eigenvector ) <<std::endl;
#endif
        }

#ifdef DEBUG
        std::cout << std::endl;
#endif

        delete [] eigenvalues_complex;
        delete [] geev_left_eigenvectors;
        delete [] geev_right_eigenvectors;
        delete [] temporary_matrix;

        return info;
    }
    else
    {
        std::cout << "ZGEEV finished with value: " << info << std::endl;
        return info;
    }
}


void DestroyEigenvectors ( std::vector<boost::numeric::ublas::vector<std :: complex<double> > * > & eigenvectors )
{
    for ( std::vector<boost::numeric::ublas::vector<std :: complex<double> > * > :: iterator iter = eigenvectors . begin ( ) ; iter != eigenvectors . end ( ) ; ++ iter )
    {
        boost::numeric::ublas::vector<std :: complex<double> > * eigenvector = ( * iter );

        delete eigenvector;
    }

    eigenvectors . clear ( );
}


int Logarithm ( math_functions::complex_hermitian_matrix & result , math_functions::complex_hermitian_matrix & base )
{
    long info;
    char jobs = 'V';
    char uplo = 'U';
    long rows = base . size1 ( );
    long columns = base . size2 ( );
    boost::numeric::ublas::zero_matrix<std :: complex<double> > zero ( rows , columns );

    double * eigenvalues_double = new double [ columns ];
    std :: complex<double> * temporary_matrix = new std :: complex<double> [ rows * columns ];

    for ( int i = 0 ; i < columns ; i ++ )
    {
        temporary_matrix [ i * columns + i ] = base ( i , i );

        for ( int j = i ; j < rows ; j ++ )
        {
            std :: complex<double> element = base ( i , j );

            temporary_matrix [ i * columns + j ] = element;
            temporary_matrix [ j * rows + i ] = conj ( element );
        }
    }

    lpp :: heev ( & jobs , & uplo , & rows , temporary_matrix , & columns , eigenvalues_double , & info );

    math_functions::complex_matrix diagonal_matrix ( rows , columns );
    diagonal_matrix = zero;

// 	Dout ( dc::notice, "Logarithm" );
    for ( int i = 0 ; i < columns ; i ++ )
    {
// 		Dout ( dc::notice, "Eigenvalue: " << eigenvalues_double [ i ] );

        if ( eigenvalues_double [ i ] > 0 )
        {
            diagonal_matrix ( i , i ) = log ( eigenvalues_double [ i ] );
        }
        else
        {
            diagonal_matrix ( i , i ) = 0;
        }
    }

    math_functions::complex_matrix transformation ( rows , columns );
    for ( int i = 0 ; i < columns ; i ++ )
    {
        for ( int j = 0 ; j < rows ; j ++ )
        {
            transformation ( i , j ) = temporary_matrix [ j * columns + i ];
        }
    }

    long int pivots [ columns ];

    lpp :: getrf ( & rows , & columns , temporary_matrix , & columns , pivots , & info );
    lpp :: getri ( & rows , temporary_matrix , & columns , pivots , & info );

    math_functions::complex_matrix inverse_transformation ( rows , columns );
    for ( int i = 0 ; i < columns ; i ++ )
    {
        for ( int j = 0 ; j < rows ; j ++ )
        {
            inverse_transformation ( i , j ) = temporary_matrix [ j * columns + i ];
        }
    }

    math_functions::complex_matrix semiresult ( rows , columns );
    semiresult = prec_prod ( transformation , diagonal_matrix );
    result = prec_prod ( semiresult , inverse_transformation );

    if ( eigenvalues_double != nullptr )
    {
        delete [] eigenvalues_double;
    }

    if ( temporary_matrix != nullptr )
    {
        delete [] temporary_matrix;
    }

    return 0;
}


int Inversion ( math_functions::complex_matrix & result , math_functions::complex_matrix & base )
{
    long info;
    long rows = base . size1 ( );
    long columns = base . size2 ( );
    long int pivots [ columns ];
    std :: complex<double> * temporary_matrix = new std :: complex<double> [ rows * columns ];
    boost::numeric::ublas::zero_matrix<std :: complex<double> > zero ( rows , columns );
    result = zero;

    for ( int i = 0 ; i < columns ; i ++ )
    {
        for ( int j = 0 ; j < rows ; j ++ )
        {
            std :: complex<double> element = base ( i , j );

            temporary_matrix [ i + rows * j ] = element;
        }
    }

    lpp :: getrf ( & rows , & columns , temporary_matrix , & columns , pivots , & info );
// 	if ( info != 0 )
// 	{
// 		return info;
// 	}

    lpp :: getri ( & rows , temporary_matrix , & columns , pivots , & info );
// 	if ( info != 0 )
// 	{
// 		return info;
// 	}

    for ( int i = 0 ; i < columns ; i ++ )
    {
        for ( int j = 0 ; j < rows ; j ++ )
        {
            std :: complex<double> element = temporary_matrix [ i + rows * j ];

            result ( i , j ) = element;
        }
    }

    return info;
}


int Inversion ( math_functions::complex_hermitian_matrix & result , math_functions::complex_hermitian_matrix & base )
{
    long info;
// 	const char jobs = 'V';
    const char uplo = 'U';
    long rows = base . size1 ( );
    long columns = base . size2 ( );
    std :: complex<double> * temporary_matrix = new std :: complex<double> [ rows * columns ];

    int k = 0;
    for ( int i = 0 ; i < columns ; i ++ )
    {
        temporary_matrix [ i * columns + i ] = base ( i , i );

        for ( int j = i ; j < rows ; j ++ )
        {
            std :: complex<double> element = base ( i , j );

            temporary_matrix [ i * columns + j ] = element;
            temporary_matrix [ j * rows + i ] = conj ( element );
        }
    }

    lpp :: potrf ( & uplo , & rows , temporary_matrix , & columns , & info );
    if ( info != 0 )
    {
        return info;
    }

    lpp :: potri ( & uplo , & rows , temporary_matrix , & columns , & info );

    k = 0;
    for ( int i = 0 ; i < columns ; i ++ )
    {
        for ( int j = i ; j < rows ; j ++ )
        {
            std :: complex<double> element = temporary_matrix [ k ];

            if ( i == j )
            {
                element =  element . real ( );
            }

            result ( i , j ) = element;
            result ( j , i ) = element;

            k ++;
        }
    }

    return info;
}


int Inversion ( math_functions::real_matrix & result , math_functions::real_matrix & base )
{
    long info;
    long rows = base . size1 ( );
    long columns = base . size2 ( );
    long int pivots [ columns ];
    double * temporary_matrix = new double [ rows * columns ];
    boost::numeric::ublas::zero_matrix<double> zero ( rows , columns );
    result = zero;

    for ( int i = 0 ; i < columns ; i ++ )
    {
        for ( int j = 0 ; j < rows ; j ++ )
        {
            double element = base ( i , j );

            temporary_matrix [ i + rows * j ] = element;
        }
    }

    lpp :: getrf ( & rows , & columns , temporary_matrix , & columns , pivots , & info );
// 	if ( info != 0 )
// 	{
// 		return info;
// 	}

    lpp :: getri ( & rows , temporary_matrix , & columns , pivots , & info );
// 	if ( info != 0 )
// 	{
// 		return info;
// 	}

    for ( int i = 0 ; i < columns ; i ++ )
    {
        for ( int j = 0 ; j < rows ; j ++ )
        {
            double element = temporary_matrix [ i + rows * j ];

            result ( i , j ) = element;
        }
    }

    return info;
}


int GeneralFunction ( math_functions::complex_hermitian_matrix & result , math_functions::complex_hermitian_matrix & base , double ( * function ) ( double ) )
{
    std::vector<double> eigenvalues;
    std::vector<boost::numeric::ublas::vector<std :: complex<double> > * > eigenvectors;
    math_tools :: EigenvaluesAndEigenvectors ( base , eigenvalues , eigenvectors );

    long rows = base . size1 ( );
    long columns = base . size2 ( );
    boost::numeric::ublas::zero_matrix<std :: complex<double> > zero ( rows , columns );

    math_functions::complex_matrix diagonal_matrix ( rows , columns );
    diagonal_matrix = zero;

    for ( int i = 0 ; i < columns ; i ++ )
    {
        double eigenvalue = eigenvalues [ i ];
        double value = function ( eigenvalue );
        diagonal_matrix ( i , i ) = value;
    }

    math_functions::complex_matrix transformation ( rows , columns );
    math_functions::complex_matrix inversion = zero;
    for ( int i = 0 ; i < columns ; i ++ )
    {
        for ( int j = 0 ; j < rows ; j ++ )
        {
            transformation ( i , j ) = ( * eigenvectors [ i ] ) [ j ];
        }
    }

    math_tools :: Inversion ( inversion , transformation );

    math_functions::complex_matrix semiresult ( rows , columns );
    noalias ( semiresult ) = prec_prod ( inversion , diagonal_matrix );
    noalias ( result ) = prec_prod ( semiresult , transformation );

    return 0;
}


int GeneralFunction ( math_functions::complex_matrix & result , math_functions::complex_hermitian_matrix & base , std :: complex<double> ( * function ) ( const std :: complex<double> & ) )
{
    std::vector<double> eigenvalues;
    std::vector<boost::numeric::ublas::vector<std :: complex<double> > * > eigenvectors;
    math_tools :: EigenvaluesAndEigenvectors ( base , eigenvalues , eigenvectors );

    long rows = base . size1 ( );
    long columns = base . size2 ( );
    boost::numeric::ublas::zero_matrix<std :: complex<double> > zero ( rows , columns );

    math_functions::complex_matrix diagonal_matrix ( rows , columns );
    diagonal_matrix = zero;

    for ( int i = 0 ; i < columns ; i ++ )
    {
        double eigenvalue = eigenvalues [ i ];
        std::complex<double> value = function ( eigenvalue );
        diagonal_matrix ( i , i ) = value;
    }

    math_functions::complex_matrix transformation ( rows , columns );
    math_functions::complex_matrix inversion = zero;
    for ( int i = 0 ; i < columns ; i ++ )
    {
        for ( int j = 0 ; j < rows ; j ++ )
        {
            transformation ( i , j ) = ( * eigenvectors [ i ] ) [ j ];
        }
    }

    math_tools :: Inversion ( inversion , transformation );

    math_functions::complex_matrix semiresult ( rows , columns );
    noalias ( semiresult ) = prec_prod ( inversion , diagonal_matrix );
    noalias ( result ) = prec_prod ( semiresult , transformation );

    return 0;
}


int GeneralFunction ( math_functions::complex_matrix & result , math_functions::complex_matrix & base , std :: complex<double> ( * function ) ( const std :: complex<double> & ) )
{
    std::vector<std :: complex<double> > eigenvalues;
    std::vector<boost::numeric::ublas::vector<std :: complex<double> > * > eigenvectors;
    math_tools :: EigenvaluesAndRightEigenvectors ( base , eigenvalues , eigenvectors );

    long rows = base . size1 ( );
    long columns = base . size2 ( );
    boost::numeric::ublas::zero_matrix<std :: complex<double> > zero ( rows , columns );

    math_functions::complex_matrix diagonal_matrix ( rows , columns );
    diagonal_matrix = zero;

    for ( int i = 0 ; i < columns ; i ++ )
    {
        std :: complex<double> eigenvalue = eigenvalues [ i ];
        std :: complex<double> value = function ( eigenvalue );
        diagonal_matrix ( i , i ) = value;
    }

    math_functions::complex_matrix transformation ( rows , columns );
    math_functions::complex_matrix inversion = zero;
    for ( int i = 0 ; i < columns ; i ++ )
    {
        for ( int j = 0 ; j < rows ; j ++ )
        {
            transformation ( i , j ) = ( * eigenvectors [ i ] ) [ j ];
        }
    }

    math_tools :: Inversion ( inversion , transformation );

    math_functions::complex_matrix semiresult ( rows , columns );
    noalias ( semiresult ) = prec_prod ( inversion , diagonal_matrix );
    noalias ( result ) = prec_prod ( semiresult , transformation );

    return 0;
}


int Sqrt ( math_functions::complex_hermitian_matrix & result , math_functions::complex_hermitian_matrix & base )
{
    return GeneralFunction ( result , base , & sqrt );
    /*
    	long info;
    	char jobs = 'V';
    	char uplo = 'U';
    	long rows = base . size1 ( );
    	long columns = base . size2 ( );
    	boost::numeric::ublas::zero_matrix<std :: complex<double> > zero ( rows , columns );

    	double * eigenvalues_double = new double [ columns ];
    	std :: complex<double> * temporary_matrix = new std :: complex<double> [ rows * columns ];

    	for ( int i = 0 ; i < columns ; i ++ )
    	{
    		temporary_matrix [ i * columns + i ] = base ( i , i );

    		for ( int j = i ; j < rows ; j ++ )
    		{
    			std :: complex<double> element = base ( i , j );

    			temporary_matrix [ i * columns + j ] = element;
    			temporary_matrix [ j * rows + i ] = conj ( element );
    		}
    	}

    	lpp :: heev ( & jobs , & uplo , & rows , temporary_matrix , & columns , eigenvalues_double , & info );

    	math_functions::complex_matrix diagonal_matrix ( rows , columns );
    	diagonal_matrix = zero;

    // 	Dout ( dc::notice, "Sqrt" );
    	for ( int i = 0 ; i < columns ; i ++ )
    	{
    // 		Dout ( dc::notice, "Eigenvalue: " << eigenvalues_double [ i ] );

    		if ( eigenvalues_double [ i ] >= 0 )
    		{
    			diagonal_matrix ( i , i ) = sqrt ( eigenvalues_double [ i ] );
    		}
    		else
    		{
    			diagonal_matrix ( i , i ) = 0;
    		}
    	}

    	math_functions::complex_matrix transformation ( rows , columns );
    	for ( int i = 0 ; i < columns ; i ++ )
    	{
    		for ( int j = 0 ; j < rows ; j ++ )
    		{
    			transformation ( i , j ) = temporary_matrix [ j * columns + i ];
    		}
    	}

    	long int pivots [ columns ];

    	lpp :: getrf ( & rows , & columns , temporary_matrix , & columns , pivots , & info );
    	lpp :: getri ( & rows , temporary_matrix , & columns , pivots , & info );

    	math_functions::complex_matrix inverse_transformation ( rows , columns );
    	for ( int i = 0 ; i < columns ; i ++ )
    	{
    		for ( int j = 0 ; j < rows ; j ++ )
    		{
    			inverse_transformation ( i , j ) = temporary_matrix [ j * columns + i ];
    		}
    	}

    	math_functions::complex_matrix semiresult ( rows , columns );
    	math_functions::complex_matrix result_as_matrix ( rows , columns );
    	boost::numeric::ublas::hermitian_adaptor<math_functions::complex_matrix, boost::numeric::ublas::upper> hermitian_adaptor ( result_as_matrix );

    	semiresult = prec_prod ( transformation , diagonal_matrix );
    	result_as_matrix = prec_prod ( semiresult , inverse_transformation );
    	result = hermitian_adaptor;

    	return 0;
    */
}


int Sqrt ( math_functions::complex_matrix & result , math_functions::complex_matrix & base )
{
// 	std :: complex<double> ( * function ) ( const std :: complex<double> & ) = & std::sqrt;
    return GeneralFunction ( result , base , & std::sqrt );
    /*
    	long info;
    	char jobs = 'V';
    	char uplo = 'U';
    	long rows = base . size1 ( );
    	long columns = base . size2 ( );
    	boost::numeric::ublas::zero_matrix<std :: complex<double> > zero ( rows , columns );

    	double * eigenvalues_double = new double [ columns ];
    	std :: complex<double> * temporary_matrix = new std :: complex<double> [ rows * columns ];

    	for ( int i = 0 ; i < columns ; i ++ )
    	{
    		temporary_matrix [ i * columns + i ] = base ( i , i );

    		for ( int j = i ; j < rows ; j ++ )
    		{
    			std :: complex<double> element = base ( i , j );

    			temporary_matrix [ i * columns + j ] = element;
    			temporary_matrix [ j * rows + i ] = conj ( element );
    		}
    	}

    	lpp :: heev ( & jobs , & uplo , & rows , temporary_matrix , & columns , eigenvalues_double , & info );

    	math_functions::complex_matrix diagonal_matrix ( rows , columns );
    	diagonal_matrix = zero;

    	for ( int i = 0 ; i < columns ; i ++ )
    	{
    		if ( eigenvalues_double [ i ] >= 0 )
    		{
    			diagonal_matrix ( i , i ) = sqrt ( eigenvalues_double [ i ] );
    		}
    		else
    		{
    			diagonal_matrix ( i , i ) = 0;
    		}
    	}

    	math_functions::complex_matrix transformation ( rows , columns );
    	for ( int i = 0 ; i < columns ; i ++ )
    	{
    		for ( int j = 0 ; j < rows ; j ++ )
    		{
    			transformation ( i , j ) = temporary_matrix [ j * columns + i ];
    		}
    	}

    	long int pivots [ columns ];

    	lpp :: getrf ( & rows , & columns , temporary_matrix , & columns , pivots , & info );
    	lpp :: getri ( & rows , temporary_matrix , & columns , pivots , & info );

    	math_functions::complex_matrix inverse_transformation ( rows , columns );
    	for ( int i = 0 ; i < columns ; i ++ )
    	{
    		for ( int j = 0 ; j < rows ; j ++ )
    		{
    			inverse_transformation ( i , j ) = temporary_matrix [ j * columns + i ];
    		}
    	}

    	math_functions::complex_matrix semiresult ( rows , columns );
    	math_functions::complex_matrix result_as_matrix ( rows , columns );
    	boost::numeric::ublas::hermitian_adaptor<math_functions::complex_matrix, boost::numeric::ublas::upper> hermitian_adaptor ( result_as_matrix );

    	semiresult = prec_prod ( transformation , diagonal_matrix );
    	result_as_matrix = prec_prod ( semiresult , inverse_transformation );
    	result = hermitian_adaptor;

    	return 0;
    */
}


double QuantumEntropy ( math_functions::complex_hermitian_matrix & input )
{
    math_functions::complex_hermitian_matrix logarithm ( input . size1 ( ) , input . size2 ( ) );

    if ( Logarithm ( logarithm , input ) )
    {
        return 0;
    }
    else
    {
        math_functions::complex_matrix before_trace ( input . size1 ( ) , input . size2 ( ) );
        boost::numeric::ublas::hermitian_adaptor<math_functions::complex_matrix, boost::numeric::ublas::upper> hermitian_adaptor ( before_trace );

        before_trace = prec_prod ( input , logarithm );

        return ( - TrOperator ( hermitian_adaptor ) );
    }
}


int Sqr ( math_functions::complex_hermitian_matrix & result , math_functions::complex_hermitian_matrix & base )
{
    result = prod ( base , base );

    return 0;
}


int Sqr ( math_functions::complex_matrix & result , math_functions::complex_matrix & base )
{
    result = prod ( base , base );

    return 0;
}


void Join ( math_functions::complex_matrix & result , math_functions::complex_matrix & base )
{
    if ( ( base.size1() != result.size2() ) && ( base.size2() != result.size1() ) )
    {
        std::cerr << "Join: Non-compatible matrices" << std::endl;
    }
    else
    {
        for ( unsigned int i = 0 ; i < base.size1() ; i ++ )
        {
            for ( unsigned int j = 0 ; j < base.size2() ; j ++ )
            {
                result ( j , i ) = conj ( base ( i , j ) );
            }
        }

    }
}


void HamiltonianContinuousQuantumRandomWalk ( math_functions::complex_hermitian_matrix & result , double phase )
{
    if ( result.size1() != result.size2() )
    {
        std::cerr << "HamiltonianContinuousQuantumRandomWalk: Illegal matrix" << std::endl;
    }
    else
    {
        boost::numeric::ublas::zero_matrix<std :: complex<double> > zero ( result.size1 ( ) , result.size2 ( ) );
        result = zero;
        std::complex<double> j ( 0 , 1 );
        std::complex<double> phase_shift = exp ( j * phase );
        std::complex<double> minus_phase_shift = exp ( - j * phase );

        for ( unsigned int i = 1 ; i < result. size1 () ; i ++)
        {
            result ( i , i -1 ) = 1;
            result ( i - 1 , i ) = 1;
        }
        result ( 0 , result.size1 () - 1 ) = phase_shift;
        result ( result.size1 () - 1 , 0 ) = minus_phase_shift;
    }
}

}
