/*! \file math_functions.cpp
\brief Mathematical functions definition file

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

#include <iostream>
#include "math_functions.h"
#include <math.h>
#include <stdio.h>

#include <gsl/gsl_sf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>

#define JAKOBI_METHOD_ITERATION_LIMIT 50
#define JAKOBI_METHOD_ABSOLUTE_NONDIAGONAL_ELEMENT_UPPER_LIMIT 0.00001

namespace math_functions
{

double absolut_value ( double a )
{
    double result = a;
    if ( result >= 0 )
    {
        return result;
    }
    else
    {
        return - result;
    }
}


int absolut_value ( int a )
{
    int result = a;
    if ( result >= 0 )
    {
        return result;
    }
    else
    {
        return - result;
    }
}

double factorial ( int i )
{
    return gsl_sf_fact ( i );
}


double pythagoras_hypotenuse ( double a , double b )
{
    return gsl_hypot ( a , b );
}


double conditional_sign ( double a , double b )
{
    return ( ( b ) >= 0.0 ? fabs ( a ) : - fabs ( a ) );
}


//#define INDEX_VECTOR(v,i,N) *(v+i)
double & INDEX_VECTOR ( double * v , int i , int N )
{
    int index = i;
    return ( * ( v + index ) );
}


//#define INDEX_MATRIX(A,i,j,N) *(A+j*N+i)
double & INDEX_MATRIX ( double * A , int i , int j , int N )
{
    int index = i * N + j;
    return ( * ( A + index ) );
}


std :: complex<double> & INDEX_MATRIX ( std :: complex<double> * A , int i , int j , int N )
{
    int index = i * N + j;
    std :: complex<double> & element = * ( A + index );
    return ( element );
}


int ROTATE ( double * a , int  i , int  j , int  k , int  l , const int N , double * s , double * tau )
{
    double g = INDEX_MATRIX ( a , i , j , N );
    double h = INDEX_MATRIX ( a , k , l , N);
    INDEX_MATRIX ( a , i , j , N ) = g - ( * s ) * ( h + g * ( * tau )  );
// 	INDEX_MATRIX ( a , j , i , N ) = INDEX_MATRIX ( a , i , j , N );
    INDEX_MATRIX ( a , k , l , N ) = h + ( * s ) * ( g - h * ( * tau ) );
// 	INDEX_MATRIX ( a , l , k , N ) = INDEX_MATRIX ( a , k , l , N );

    return 0;
}


int Jakobi_eigenvalue_method_with_eigenvectors ( const int N , double * A , double * eigenvalue , double * eigenvectors , int * number_of_iterations )
{
    int i;
    int j;
    int iq;
    int ip;

    double tresh;
    double theta;
    double tau;
    double t;
    double sm;
    double s;
    double h;
    double g;
    double c;
    double * b;
    double * z;

    b = ( double * ) malloc ( N * sizeof ( double ) );
    z = ( double * ) malloc ( N * sizeof ( double ) );

    for ( ip = 0 ; ip < N ; ip ++ )
    {
        for ( iq = 0 ; iq < N ; iq ++ )
        {
            INDEX_MATRIX ( eigenvectors , ip , iq , N ) = 0.0;
        }
        INDEX_MATRIX ( eigenvectors , ip , ip , N ) = 1.0;
    }

    for ( ip = 0 ; ip < N ; ip ++ )
    {
        b [ ip ] = INDEX_MATRIX ( A , ip , ip , N );
        INDEX_VECTOR ( eigenvalue , ip , N ) = INDEX_MATRIX ( A , ip , ip , N );
        z [ ip ] = 0.0;
    }

    * number_of_iterations = 0;

    for ( i = 0 ; i < JAKOBI_METHOD_ITERATION_LIMIT ; i ++ )
    {
        sm = 0.0;
        for ( ip = 0 ; ip <= N - 2 ; ip ++ )
        {
            for ( iq = ip + 1 ; iq <= N - 1 ; iq ++ )
            {
                sm += fabs ( INDEX_MATRIX ( A , ip , iq , N ) );
            }
        }

        if ( sm <= JAKOBI_METHOD_ABSOLUTE_NONDIAGONAL_ELEMENT_UPPER_LIMIT )
        {
            free ( b );
            free ( z );

            return 0;
        }

        if ( i < 4 )
        {
            tresh = 0.2 * sm / ( N * N );
        }
        else
        {
            tresh = 0.0;
        }

        for ( ip = 0 ; ip <= N - 2 ; ip ++ )
        {
            for ( iq = ip + 1 ; iq <= N - 1 ; iq ++ )
            {
                g = 100.0 * fabs ( INDEX_MATRIX ( A , ip , iq , N ) );

                if ( ( ( i > 4 ) && ( fabs ( eigenvalue [ ip ] + g ) == fabs ( eigenvalue [ ip ] ) ) ) && ( fabs ( eigenvalue [ iq ] + g ) == fabs ( eigenvalue [ iq ] ) ) )
                {
                    INDEX_MATRIX ( A , ip , iq , N ) = 0.0;
                }
                else
                {
                    if ( fabs ( INDEX_MATRIX ( A , ip , iq , N ) ) > tresh )
                    {
                        h = eigenvalue [ iq ] - eigenvalue [ ip ];
                        if ( ( fabs ( h ) + g ) == fabs ( h ) )
                        {
                            t = INDEX_MATRIX ( A , ip , iq , N ) / h;
                        }
                        else
                        {
                            theta = h / 2.0 / INDEX_MATRIX ( A , ip , iq , N );

                            if ( theta > 0.0 )
                            {
                                t = - theta + sqrt ( 1.0 + theta * theta );
                            }
                            else
                            {
                                t = - theta - sqrt ( 1.0 + theta * theta );
                            }
                        }
                        c = 1.0 / sqrt ( 1.0 + t * t );
                        s = t * c;
                        tau = ( 1.0 - c ) / s;
                        h = t * INDEX_MATRIX ( A , ip , iq , N );
                        z [ ip ] -= h;
                        z [ iq ] += h;
                        eigenvalue [ ip ] -= h;
                        eigenvalue [ iq ] += h;
                        INDEX_MATRIX ( A , ip , iq , N ) = 0.0;
                        INDEX_MATRIX ( A , iq , ip , N ) = 0.0;

                        for ( j = 0 ; j <= ip -1 ; j ++ )
                        {
                            ROTATE ( A , j , ip , j , iq , N , & s , & tau );
                        }
                        for ( j = ip + 1 ; j  <= iq - 1 ; j ++ )
                        {
                            ROTATE (  A , ip , j , j , iq , N , & s , & tau );
                        }
                        for ( j = iq + 1 ; j < N ; j ++ )
                        {
                            ROTATE ( A , ip , j , iq , j , N , & s , & tau );
                        }
                        for ( j = 0 ; j < N ; j ++ )
                        {
                            ROTATE ( eigenvectors , j , ip , j , iq , N , & s , & tau );
                        }
                        ( * number_of_iterations ) ++;
                    }
                }
            }
        }
        for ( ip = 0 ; ip < N ; ip ++ )
        {
            b [ ip ] += z [ ip ];
            eigenvalue [ ip ] = b [ ip ];
            z [ ip ] = 0.0;
        }
    }
    fprintf ( stdout , "Jakobi_eigenvalue_method: Error! \n" );
    return 1;
}


int Jakobi_eigenvalue_method ( const int N , double * A , double * eigenvalue , int * number_of_iterations )
{
    int i;
    int j;
    int iq;
    int ip;

    double tresh;
    double theta;
    double tau;
    double t;
    double sm;
    double s;
    double h;
    double g;
    double c;
    double * b;
    double * z;

    b = ( double * ) malloc ( N * sizeof ( double ) );
    z = ( double * ) malloc ( N * sizeof ( double ) );

    for ( ip = 0 ; ip < N ; ip ++ )
    {
        b [ ip ] = INDEX_MATRIX ( A , ip , ip , N );
        INDEX_VECTOR ( eigenvalue , ip , N ) = INDEX_MATRIX ( A , ip , ip , N );
        z [ ip ] = 0.0;
    }

    * number_of_iterations = 0;

    for ( i = 0 ; i < JAKOBI_METHOD_ITERATION_LIMIT ; i ++ )
    {
        sm = 0.0;
        for ( ip = 0 ; ip <= N - 2 ; ip ++ )
        {
            for ( iq = ip + 1 ; iq <= N - 1 ; iq ++ )
            {
                sm += fabs ( INDEX_MATRIX ( A , ip , iq , N ) );
            }
        }

        if ( sm <= JAKOBI_METHOD_ABSOLUTE_NONDIAGONAL_ELEMENT_UPPER_LIMIT )
        {
            free ( b );
            free ( z );

            return 0;
        }

        if ( i < 4 )
        {
            tresh = 0.2 * sm / ( N * N );
        }
        else
        {
            tresh = 0.0;
        }

        for ( ip = 0 ; ip <= N - 2 ; ip ++ )
        {
            for ( iq = ip + 1 ; iq <= N - 1 ; iq ++ )
            {
                g = 100.0 * fabs ( INDEX_MATRIX ( A , ip , iq , N ) );

                if ( ( i > 4 ) && ( fabs ( eigenvalue [ ip ] + g ) == fabs ( eigenvalue [ ip ] ) ) && ( fabs ( eigenvalue [ iq ] + g ) == fabs ( eigenvalue [ iq ] ) ) )
                {
                    INDEX_MATRIX ( A , ip , iq , N ) = 0.0;
                }
                else
                {
                    if ( fabs ( INDEX_MATRIX ( A , ip , iq , N ) ) > tresh )
                    {
                        h = eigenvalue [ iq ] - eigenvalue [ ip ];
                        if ( fabs ( h ) + g == fabs ( h ) )
                        {
                            t = INDEX_MATRIX ( A , ip , iq , N ) / h;
                        }
                        else
                        {
                            theta = h / 2.0 / INDEX_MATRIX ( A , ip , iq , N );

                            if ( theta > 0.0 )
                            {
                                t = - theta + sqrt ( 1.0 + theta * theta );
                            }
                            else
                            {
                                t = - theta - sqrt ( 1.0 + theta * theta );
                            }
                        }
                        c = 1.0 / sqrt ( 1.0 + t * t );
                        s = t * c;
                        tau = ( 1.0 - c ) / s;
                        h = t * INDEX_MATRIX ( A , ip , iq , N );
                        z [ ip ] -= h;
                        z [ iq ] += h;
                        eigenvalue [ ip ] -= h;
                        eigenvalue [ iq ] += h;
                        INDEX_MATRIX ( A , ip , iq , N ) = 0.0;
                        INDEX_MATRIX ( A , iq , ip , N ) = 0.0;

                        for ( j = 0 ; j <= ip -1 ; j ++ )
                        {
                            ROTATE ( A , j , ip , j , iq , N , & s , & tau );
                        }
                        for ( j = ip + 1 ; j  <= iq - 1 ; j ++ )
                        {
                            ROTATE (  A , ip , j , j , iq , N , & s , & tau );
                        }
                        for ( j = iq + 1 ; j < N ; j ++ )
                        {
                            ROTATE ( A , ip , j , iq , j , N , & s , & tau );
                        }
                        ( * number_of_iterations ) ++;
                    }
                }
            }
        }
        for ( ip = 0 ; ip < N ; ip ++ )
        {
            b [ ip ] += z [ ip ];
            eigenvalue [ ip ] = b [ ip ];
            z [ ip ] = 0.0;
        }
    }
    fprintf ( stdout , "Jakobi_eigenvalue_method: Error! \n" );
    return 1;
}


int Householder_reduction ( const int N , double * A , double * diagonal , double * extradiagonal )
{
    int l , k , j = 0 , i;
    double scale , hh , h , g , f;

    for ( i = N - 1 ; i >= 1 ; i -- )
    {
        l = i - 1;
        scale = 0.0;
        h = 0.0;

        if ( l > 0 )
        {
            for ( k = 0 ; k <= l ; k ++ )
            {
                scale += fabs ( INDEX_MATRIX ( A , i , k , N ) );
            }

            if ( scale == 0.0 )
            {
                extradiagonal [ i ] = INDEX_MATRIX ( A , i , l , N );
            }
            else
            {
                for ( k = 0 ; k <= l ; k ++ )
                {
                    INDEX_MATRIX ( A , i , k , N ) /= scale;
                    h += ( INDEX_MATRIX ( A , i , k , N ) * INDEX_MATRIX ( A , i , k , N ) );
                }
            }
            f = INDEX_MATRIX ( A , i , l , N );
            g = ( f >= 0.0 ? - sqrt ( h ) : sqrt ( h ) );
            extradiagonal [ i ] = scale * g;
            h -= f * g;
            INDEX_MATRIX ( A , i , l , N ) = f - g;
            f = 0.0;

            for ( j = 0 ; j <= l ; j ++ )
            {
                INDEX_MATRIX ( A , j , i , N ) = INDEX_MATRIX ( A , i , j , N ) / h;
                g = 0.0;

                for ( k = 0 ; k <= j ; k ++ )
                {
                    g += ( INDEX_MATRIX ( A , j , k , N ) * INDEX_MATRIX ( A , i , k , N ) );
                }
                for ( k = j +1 ; k <= l ; k ++ )
                {
                    g += ( INDEX_MATRIX ( A , k , j , N ) * INDEX_MATRIX ( A , i , k , N ) );
                }
                extradiagonal [ j ] = g / h;
                f += extradiagonal [ j ] * INDEX_MATRIX ( A , i , j , N );
            }
            hh = f / ( 2.0 * h);

            for ( j = 0 ; j <= l ; j ++ )
            {
                f = INDEX_MATRIX ( A , i , j , N );
                g = extradiagonal [ j ] - hh * f;
                extradiagonal [ j ] = g;

                for ( k = 0 ; k <= j ; k ++ )
                {
                    INDEX_MATRIX ( A , j , k , N ) -= ( f * extradiagonal [ k ] + g * INDEX_MATRIX ( A , i , k , N ) );
                }
            }
        }
        else
        {
            extradiagonal [ i ] = INDEX_MATRIX ( A , i , l , N );
        }
        diagonal [ i ] = h;
    }

    diagonal [ i ] = 0.0;
    extradiagonal [ i ] = 0.0;

    for ( i = 0 ; i <= N - 1 ; i ++ )
    {
        l = i - 1;
        if ( diagonal [ i ] != 0.0 )
        {
            for ( j = 0 ; j <= l ; j ++ )
            {
                g = 0.0;

                for ( k = 0 ; k <= l ; k ++ )
                {
                    g += ( INDEX_MATRIX ( A , i , k , N ) * INDEX_MATRIX ( A , k , j , N ) );
                }
                for ( k = 0 ; k <= l ; k ++ )
                {
                    INDEX_MATRIX ( A , k , j , N ) -= g * INDEX_MATRIX ( A , k , i , N );
                }
            }
        }
        diagonal [ i ] = INDEX_MATRIX ( A , i , j , N );
        INDEX_MATRIX ( A , i , i , N ) = 1.0;

        for ( j = 0 ; j <= l ; j ++ )
        {
            INDEX_MATRIX ( A , j , i , N ) = 0.0;
            INDEX_MATRIX ( A , i , j , N ) = 0.0;
        }
    }
    return 0;
}


int QL_algorithm_for_tridiagonal_with_eigenvectors ( int N , double * diagonal , double * extradiagonal , double * eigenvectors )
{
    int m , l , iter , i , k;
    double s , r , p , g , f , dd , c , b;

    for ( i = 1 ; i <= N - 1 ; i ++ )
    {
        extradiagonal [ i - 1 ] = extradiagonal [ i ];
    }
    extradiagonal [ N - 1 ] = 0.0;

    for ( l = 0 ; l <= N - 1 ; l ++ )
    {
        iter = 0;
        do
        {
            for ( m = l ; m <= N - 2 ; m ++ )
            {
                dd = fabs ( diagonal [ m ] ) + fabs ( diagonal [ m + 1 ] );
                if ( ( fabs ( extradiagonal [ m ] ) + dd ) == dd )
                {
                    break;
                }
            }
            if ( m != l )
            {
                if ( iter ++ == 30 )
                {
                    fprintf ( stdout , "QL_algorithm_for_tridiagonal_with_eigenstd :: vectors: Error! \n" );
                }
                g = ( diagonal [ l + 1 ] - diagonal [ l ] ) / ( 2.0 * extradiagonal [ l ] );
                r = pythagoras_hypotenuse ( g , 1.0 );
                g = diagonal [ m ] - diagonal [ l ] + extradiagonal [ l ] / ( g + conditional_sign ( r , g ) );
                s = 1.0;
                c = 1.0;
                p = 0.0;

                for ( i = m - 1 ; i >= l ; i -- )
                {
                    f = s * extradiagonal [ i ];
                    b = c * extradiagonal [ i ];
                    r = pythagoras_hypotenuse ( f , g );
                    extradiagonal [ i + 1 ] = r;
                    if ( r == 0.0 )
                    {
                        diagonal [ i + 1 ] -= p;
                        extradiagonal [ m ] = 0.0;
                        break;
                    }

                    s = f / r;
                    c = g / r;
                    g = diagonal [ i + 1 ] - p;
                    r = ( diagonal [ i ] - g ) * s + 2.0 * c * b;
                    p = s * r;
                    diagonal [ i + 1 ] = g + p;
                    g = c * r - b;

                    for ( k = 0 ; k <= N - 1 ; k ++ )
                    {
                        f = INDEX_MATRIX ( eigenvectors , k , i + 1 , N );
                        INDEX_MATRIX ( eigenvectors , k , i + 1 , N ) = s * INDEX_MATRIX ( eigenvectors , k , i , N ) + c * f;
                        INDEX_MATRIX ( eigenvectors , k , i , N ) = c * INDEX_MATRIX ( eigenvectors , k , i , N ) - s * f;
                    }
                }
                if ( ( r == 0.0 ) && ( i >= 0 ) )
                {
                    continue;
                }
                diagonal [ l ] -= p;
                extradiagonal [ l ] = g;
                extradiagonal [ m ] = 0.0;
            }
        }
        while ( m != l );
    }
    return 0;
}


int QL_algorithm_and_household_reduction_with_eigenvectors ( int N , double * A , double * eigenvalues , double * eigenvectors  )
{
    double * extradiagonal = ( double * ) malloc ( sizeof ( double ) * N );

    Householder_reduction ( N , A , eigenvalues , extradiagonal );

    QL_algorithm_for_tridiagonal_with_eigenvectors ( N , eigenvalues , extradiagonal , eigenvectors );

    free ( extradiagonal );

    return 0;
}


int RealMatrixMultiplication ( double * A , int columns_A , int rows_A , double * B , int columns_B , int rows_B , double * C , int columns_C , int rows_C )
{
    if ( ( columns_A != rows_B ) || ( columns_B != columns_C ) || ( rows_C != rows_A ) )
    {
        printf ( "RealMatrixMultiplication :Incorrect matrices." );
        return 1;
    }
    else
    {
        for ( int i = 0 ; i < columns_B ; i ++ )
        {
            for ( int j = 0 ; j < rows_A ; j ++ )
            {
                INDEX_MATRIX ( C , j , i , columns_B ) = 0.0;
                for ( int k = 0 ; k < columns_A ; k ++ )
                {
                    INDEX_MATRIX ( C , j , i , columns_C ) += ( INDEX_MATRIX ( A , j , k , columns_A ) * INDEX_MATRIX ( B , k , i , columns_B ) );
                }
            }
        }
        return 0;
    }
}


int ComplexMatrixMultiplication ( std :: complex<double> * A , int columns_A , int rows_A , std :: complex<double> * B , int columns_B , int rows_B , std :: complex<double> * C , int columns_C , int rows_C )
{
    if ( ( columns_A != rows_B ) || ( columns_B != columns_C ) || ( rows_C != rows_A ) )
    {
        printf ( "RealMatrixMultiplication :Incorrect matrices." );
        return 1;
    }
    else
    {
        for ( int i = 0 ; i < columns_B ; i ++ )
        {
            for ( int j = 0 ; j < rows_A ; j ++ )
            {
                std :: complex<double> & element = INDEX_MATRIX ( C , i , j , columns_C );
                element = 0.0;
                for ( int k = 0 ; k < columns_A ; k ++ )
                {
                    std :: complex<double> & element_1 = INDEX_MATRIX ( A , i , k , columns_A );
                    std :: complex<double> & element_2 = INDEX_MATRIX ( B , k , j , columns_B );
                    element += ( element_1 * element_2 );
                }
            }
        }
        return 0;
    }
}


int ComplexErrorPrintMatrix ( std :: complex<double> * A , int columns_A , int rows_A )
{
    for ( int count_i = 0 ; count_i < columns_A ; count_i ++ )
    {
        for ( int count_j = 0 ; count_j < rows_A ; count_j ++ )
        {
            int index = count_i * columns_A + count_j;
            std :: complex<double> & element = A [ index ];
            fprintf ( stderr , "%.2f +%.2f i" , real ( element ) , imag ( element ) );
            if ( count_j < columns_A - 1 )
            {
                fprintf ( stderr , " ," );
            }
        }
        fprintf ( stderr , "\n" );
    }
    return 0;
}


int ComplexMatrixTransposition ( std :: complex<double> * A , int columns_A , int rows_A , std :: complex<double> * B , int columns_B , int rows_B )
{
    if ( ( columns_A != rows_B ) || ( rows_A != columns_B ) )
    {
        printf ( "RealMatrixMultiplication :Incorrect matrices." );
        return 1;
    }
    else
    {
        for ( int i = 0 ; i < columns_A ; i ++ )
        {
            for ( int j = 0 ; j < rows_A ; j ++ )
            {
                INDEX_MATRIX ( B , j , i , columns_B ) = INDEX_MATRIX ( A , i , j , columns_A );
            }
        }
        return 0;
    }
}


int ComplexMatrixConjugation ( std :: complex<double> * A , int columns_A , int rows_A , std :: complex<double> * B , int columns_B , int rows_B )
{
    if ( ( columns_A != columns_B ) || ( rows_A != rows_B ) )
    {
        printf ( "RealMatrixMultiplication :Incorrect matrices." );
        return 1;
    }
    else
    {
        for ( int i = 0 ; i < columns_A ; i ++ )
        {
            for ( int j = 0 ; j < rows_A ; j ++ )
            {
                std :: complex<double> element = INDEX_MATRIX ( A , i , j , columns_A );
                INDEX_MATRIX ( B , i , j , columns_B ) = conj ( element );
            }
        }
        return 0;
    }
}


int ComplexMatrixHermitianConjugation ( std :: complex<double> * A , int columns_A , int rows_A , std :: complex<double> * B , int columns_B , int rows_B )
{
    if ( ( columns_A != rows_B ) || ( rows_A != columns_B ) )
    {
        printf ( "RealMatrixMultiplication :Incorrect matrices." );
        return 1;
    }
    else
    {
        for ( int i = 0 ; i < columns_A ; i ++ )
        {
            for ( int j = 0 ; j < rows_A ; j ++ )
            {
                std :: complex<double> element = INDEX_MATRIX ( A , j , i , columns_A );
                INDEX_MATRIX ( B , i , j , columns_B ) = conj ( element );
            }
        }
        return 0;
    }
}


long GetDigit ( long number , long decimals , long rank )
{
    long divider_1 = ( long ) pow ( decimals , rank );
    long divider_2 = ( long ) pow ( decimals , rank + 1 );

    long modulus_1 = number % divider_1;
    long modulus_2 = number % divider_2;

    return ( modulus_2 - modulus_1 ) / divider_1;
}


double sqr ( double x )
{
    return x * x;
}


double cube ( double x )
{
    return x * x * x;
}


double sqrt_inverse ( double x )
{
    return 1.0 / sqrt ( x );
}


double Signum ( double x )
{
	return GSL_SIGN (x);
}


double HeavysideFunction ( double x )
{
    if ( x > 0 )
    {
        return 1;
    }
    else
    {
        return 0;
    }
}


double Factorial ( int i )
{
    if ( i < 0 )
    {
        return 1.0;
    }
    else
    {
        return exp ( gsl_sf_lngamma ( i + 1 ) );
    }
}

double LogarithmusFactorial ( int i )
{
    if ( i < 0 )
    {
        return 0.0;
    }
    else
    {
        return gsl_sf_lngamma ( i + 1 );
    }
}


double & ScalarMultiplication ( boost::numeric::ublas::vector<double> & vector1 , boost::numeric::ublas::vector<double> & vector2 )
{
    double * result = new double;
    * result = 0;

    if ( vector1.size() != vector2.size() )
    {
        std::cout << "ScalarMultiplication: Incorrect sizes" << std :: endl;
    }
    else
    {
        for ( int i = 0 ; i < ( int ) vector1.size ( ) ; i ++ )
        {
            * result += vector1 [ i ] * vector2 [ i ];
        }
    }

    return * result;
}


std :: complex<double> & ScalarMultiplication ( boost::numeric::ublas::vector<std :: complex<double> > & vector1 , boost::numeric::ublas::vector<std :: complex<double> > & vector2 )
{
    std :: complex<double> * result = new std :: complex<double>;
    * result = ( std :: complex<double> ) 0;

    if ( vector1.size() != vector2.size() )
    {
        std::cout << "ScalarMultiplication: Incorrect sizes" << std :: endl;
    }
    else
    {
        for ( int i = 0 ; i < ( int ) vector1.size ( ) ; i ++ )
        {
            * result += vector1 [ i ] * vector2 [ i ];
        }
    }

    return * result;
}


double Combination ( int n , int m )
{
    return gsl_sf_choose ( n , m );
}


double LogarithmusCombination ( int n , int m )
{
    return gsl_sf_lnchoose ( n , m );
}


int minimum ( int a , int b )
{
    if ( a < b )
    {
        return a;
    }
    else
    {
        return b;
    }
}

double GeneralizedHarmonicNumber ( unsigned int N , double s )
{
	double sum = 0;
	for ( unsigned int i = 1 ; i <= N ; i ++ )
	{
		sum += 1.0 / pow ( i , s );
	}

	return sum;
}

}
