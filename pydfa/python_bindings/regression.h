/***************************************************************************
 *   Copyright (C) 2012 by Hynek Lavička   *
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

#include <list>
#include <vector>
#include <iostream>
#include <math.h>
#include <unordered_map>

#include <boost/thread/shared_mutex.hpp>

#include "math_functions.h"
#include "stl_wrapper_for_lapack.h"

#ifndef REGRESSION_H
#define REGRESSION_H

namespace statistical_functions
{

struct KeyHasherMatricesForRegression
{
	std::size_t operator()(const std::tuple< unsigned int , unsigned int> & key ) const
	{
		using std::size_t;
		using std::hash;

		return ( ( hash<unsigned int>() ( std::get<0> ( key ) ) ^ ( hash<unsigned int>() ( std::get<1> ( key ) ) << 1) ) >> 1 );
	}
};

double HarmonicFunction ( double x , unsigned int order , double domain );

// void PolynomialRegression ( std::list<double> & y_data , unsigned int order_of_method , std::list<double> & coefficients );
void PrepareMatricesForRegression ( unsigned int no_of_data  , unsigned int order_of_method , math_functions::real_matrix & transformation_y_data_cached  );

template<typename STORAGE , typename STORAGE_ITERATOR>
void PolynomialRegression ( STORAGE & y_data , unsigned int order_of_method , std::list<double> & coefficients , math_functions::real_matrix * transformation_y_data_cached = nullptr )
{
	if ( ( y_data.size() <= order_of_method ) )
	{
		std::cerr << "PolynomialRegression: Illegal data or illegal order of method" << std::endl;
	}
	else
	{
		unsigned int no_of_data = y_data.size();

		coefficients.clear();

		math_functions::real_column_vector Y_vector ( no_of_data );
		math_functions::real_column_vector coefficients_vector ( order_of_method + 1 );
		STORAGE_ITERATOR y_data_iter = y_data.begin();

		if ( transformation_y_data_cached == nullptr )
		{
			math_functions::real_matrix transformation_y_data ( order_of_method + 1 , no_of_data );
			math_functions::real_matrix X_T_X_matrix ( order_of_method  + 1, order_of_method  + 1 );
			math_functions::real_matrix X_T_X_matrix_invesion ( order_of_method + 1, order_of_method + 1 );
			boost::numeric::ublas::zero_matrix<double> zero1 ( order_of_method + 1 , order_of_method + 1 );
			boost::numeric::ublas::zero_matrix<double> zero2 ( order_of_method + 1 , no_of_data );
			math_functions::real_matrix X_matrix_T ( order_of_method  + 1, no_of_data );
			X_T_X_matrix = zero1;
			X_T_X_matrix_invesion = zero1;
			X_matrix_T = zero2;

			for ( unsigned int i = 0 ; i < no_of_data ; i ++ )
			{
				for ( unsigned int j = 0 ; j <= order_of_method ; j ++ )
				{
					double x = i;
					X_matrix_T ( j , i ) = pow ( x , ( double ) j );

					for ( unsigned int k = 0 ; k <= order_of_method ; k ++ )
					{
						X_T_X_matrix ( j , k ) += pow ( x , ( double ) j ) * pow ( x , ( double ) k );
					}
				}
				Y_vector ( i ) = ( * y_data_iter );
				++ y_data_iter;
			}

			math_tools::Inversion ( X_T_X_matrix_invesion , X_T_X_matrix );
			boost::numeric::ublas::axpy_prod ( X_T_X_matrix_invesion , X_matrix_T , transformation_y_data , true );
			boost::numeric::ublas::axpy_prod ( transformation_y_data , Y_vector, coefficients_vector , true);
		}
		else
		{
			for ( unsigned int i = 0 ; i < no_of_data ; i ++ )
			{
				Y_vector ( i ) = ( * y_data_iter );
				++ y_data_iter;
			}
			boost::numeric::ublas::axpy_prod ( * transformation_y_data_cached , Y_vector, coefficients_vector , true);
		}

		for ( unsigned int i = 0 ; i <= order_of_method ; i ++ )
		{
			coefficients.push_back ( coefficients_vector ( i ) );
		}
	}
}

template <class TYPENAME>
class default_NAN_float
{
	private:
		TYPENAME _value;

	public:
		default_NAN_float ()
		: _value ( std::numeric_limits<TYPENAME>::quiet_NaN() )
		{
		}

		TYPENAME & operator() ()
		{
			return _value;
		}
};


template <class TYPENAME>
TYPENAME GetPowerMapping ( std::unordered_map<std::tuple< unsigned int , unsigned int> , TYPENAME , KeyHasherMatricesForRegression> & power , unsigned int position , unsigned int exponent
#ifdef THREADING
, boost::shared_mutex & power_mutex 
#endif
)
{
	std::tuple< unsigned int , unsigned int> index = std::make_tuple ( position , exponent );
#ifdef THREADING
	power_mutex.lock_shared();
#endif
	typename std::unordered_map<std::tuple< unsigned int , unsigned int> , TYPENAME , KeyHasherMatricesForRegression>::iterator iterator_power;
	iterator_power = power.find ( index );

	if ( iterator_power == power.end() )
	{
#ifdef THREADING
		power_mutex.unlock_shared();
#endif
		TYPENAME result = math_functions::multiprecision::power ( static_cast<TYPENAME> ( position ) , static_cast<TYPENAME> ( exponent ) );
#ifdef THREADING
		power_mutex.lock();
#endif
		power [ index ] = result;
#ifdef THREADING
		power_mutex.unlock();
#endif
		return result;
	}
	else
	{
		TYPENAME result = ( * iterator_power ).second;
#ifdef THREADING
		power_mutex.unlock_shared();
#endif
		return result;
	}
}

/*
template <class MATRIX_TYPE>
MATRIX_TYPE GetMatrixInversionMapping ( std::unordered_map< std::tuple< unsigned int , unsigned int> , MATRIX_TYPE , KeyHasherMatricesForRegression> & matrix_inversions , unsigned int order , unsigned int no_of_data ,  MATRIX_TYPE & matrix_to_invert
#ifdef THREADING
, boost::shared_mutex & matrix_inversions_mutex 
#endif
)
{
	typename std::unordered_map< std::tuple< unsigned int , unsigned int> , MATRIX_TYPE , KeyHasherMatricesForRegression>::iterator iterator_matrix;

#ifdef THREADING
	matrix_inversions_mutex.lock_shared();
	std::tuple< unsigned int , unsigned int> index = std::make_tuple ( order , no_of_data );
	iterator_matrix = matrix_inversions.find ( index );
#endif

	if ( iterator_matrix == matrix_inversions.end() )
	{
#ifdef THREADING
		matrix_inversions_mutex.unlock_shared();
#endif
		MATRIX_TYPE inversion = math_tools::Inversion<MATRIX_TYPE> ( matrix_to_invert );
#ifdef THREADING
		matrix_inversions_mutex.lock();
#endif
		matrix_inversions [ index ] = inversion;
#ifdef THREADING
		matrix_inversions_mutex.unlock();
#endif
		return inversion;
	}
	else
	{
		MATRIX_TYPE inversion = ( * iterator_matrix ).second;
#ifdef THREADING
		matrix_inversions_mutex.unlock_shared();
#endif
		return inversion;
	}
}
*/

template <class MATRIX_TYPE>
void PrepareMatricesForRegression ( unsigned int no_of_data , unsigned int order_of_method , MATRIX_TYPE & transformation_y_data_cached  )
{
	typedef typename MATRIX_TYPE::value_type TYPE_NAME;

	static std::unordered_map<std::tuple< unsigned int , unsigned int> , TYPE_NAME , KeyHasherMatricesForRegression> powers;
#ifdef THREADING
	static boost::shared_mutex power_mutex;
#endif

	static std::unordered_map<std::tuple< unsigned int , unsigned int> , MATRIX_TYPE , KeyHasherMatricesForRegression> matrix_transforms;
#ifdef THREADING
	static boost::shared_mutex matrix_inversions_mutex;
#endif

	std::tuple< unsigned int , unsigned int> index = std::make_tuple ( order_of_method , no_of_data );
#ifdef THREADING
	matrix_inversions_mutex.lock_shared();
#endif
	typename std::unordered_map< std::tuple< unsigned int , unsigned int> , MATRIX_TYPE , KeyHasherMatricesForRegression>::iterator iterator_matrix;
	iterator_matrix = matrix_transforms.find ( index );

	if ( iterator_matrix == matrix_transforms.end() )
	{
#ifdef THREADING
		matrix_inversions_mutex.unlock_shared();
#endif

		transformation_y_data_cached.resize ( order_of_method + 1 , no_of_data );
		MATRIX_TYPE X_T_X ( order_of_method  + 1, order_of_method  + 1 );
		MATRIX_TYPE X_matrix_T ( order_of_method  + 1, no_of_data );

		{
			for ( unsigned int i = 0 ; i <= order_of_method ; ++ i )
			{
				for ( unsigned int j = 0 ; j <= order_of_method ; ++ j )
				{
					X_T_X ( i , j ) = static_cast< TYPE_NAME > ( 0 );
				}
			}
		}
#ifdef DEBUG
		std::cout << X_T_X << std::endl;
#endif

		for ( unsigned int i = 0 ; i < no_of_data ; ++ i )
		{
			for ( unsigned int j = 0 ; j <= order_of_method ; ++ j )
			{
				// TYPE_NAME x = i;
				TYPE_NAME power_1 = GetPowerMapping<TYPE_NAME> ( powers , i , j
#ifdef THREADING
				, power_mutex
#endif
				);

				X_matrix_T ( j , i ) = power_1;
	// 			pow ( x , ( double ) j );

				for ( unsigned int k = 0 ; k <= order_of_method ; k ++ )
				{
					TYPE_NAME power_2 = GetPowerMapping<TYPE_NAME> ( powers ,  i , k
#ifdef THREADING
					, power_mutex
#endif
					);

					TYPE_NAME result = power_1 * power_2;
	// 				pow ( x , ( double ) j ) * pow ( x , ( double ) k );
#ifdef DEBUG
	// 				std::cout << std::setprecision(std::numeric_limits<TYPE_NAME>::digits10 ) << "Power " << power_1 << " " << power_2 << " " << result << std::endl;
#endif
					X_T_X ( j , k ) += result;
				}
			}
		}

		MATRIX_TYPE inversion = math_tools::Inversion<MATRIX_TYPE> ( X_T_X );
		transformation_y_data_cached = boost::numeric::ublas::prod ( inversion , X_matrix_T );
	// 	boost::numeric::ublas::axpy_prod ( inversion , X_matrix_T , transformation_y_data_cached , true );

#ifdef THREADING
		matrix_inversions_mutex.lock();
#endif
		matrix_transforms [ index ] = transformation_y_data_cached;
#ifdef THREADING
		matrix_inversions_mutex.unlock();
#endif

#ifdef DEBUG
		std::cout << X_T_X << std::endl;
		std::cout << inversion << std::endl;
		std::cout << prod ( X_T_X , inversion ) << std::endl;
		std::cout << X_matrix_T << std::endl;
		std::cout << transformation_y_data_cached << std::endl;
#endif
	}
	else
	{
		transformation_y_data_cached = ( * iterator_matrix ).second;
#ifdef THREADING
		matrix_inversions_mutex.unlock_shared();
#endif
	}
}


template< typename STORAGE , typename INTERMEDIATE_TYPE , typename MATRIX_TYPE , typename VECTOR_TYPE >
void PolynomialRegressionTemplate ( STORAGE & y_data , unsigned int order_of_method , std::list<INTERMEDIATE_TYPE> & coefficients , MATRIX_TYPE * transformation_y_data_cached = nullptr )
{
	typedef typename STORAGE::iterator STORAGE_ITERATOR;
	typedef typename MATRIX_TYPE::value_type INTERNAL_TYPE;

	if ( ( y_data.size() <= order_of_method ) )
	{
		std::cerr << "PolynomialRegression: Illegal data or illegal order of method" << std::endl;
	}
	else
	{
		unsigned int no_of_data = y_data.size();

		coefficients.clear();

		VECTOR_TYPE Y_vector ( no_of_data );
		VECTOR_TYPE coefficients_vector ( order_of_method + 1 );
		STORAGE_ITERATOR y_data_iter = y_data.begin();

		MATRIX_TYPE transformation_y_data_calculated;
		if ( transformation_y_data_cached == nullptr )
		{
			PrepareMatricesForRegression<MATRIX_TYPE> ( no_of_data , order_of_method , transformation_y_data_calculated );
			transformation_y_data_cached = & transformation_y_data_calculated;
		}

		for ( unsigned int i = 0 ; i < no_of_data ; ++ i )
		{
			Y_vector ( i ) = static_cast<INTERNAL_TYPE> ( ( * y_data_iter ) );
			++ y_data_iter;
		}
// 		boost::numeric::ublas::axpy_prod ( * transformation_y_data_cached , Y_vector, coefficients_vector , true);
		coefficients_vector = boost::numeric::ublas::prod ( ( * transformation_y_data_cached ) , Y_vector );

		for ( unsigned int i = 0 ; i <= order_of_method ; ++ i )
		{
			coefficients.push_back ( static_cast <INTERMEDIATE_TYPE> ( coefficients_vector ( i ) ) );
		}
	}
}

// this function take one additional double variable that is for polynomial regression redundant
template<typename STORAGE , typename STORAGE_ITERATOR>
void PolynomialRegression ( STORAGE & y_data , unsigned int order_of_method , std::list<double> & coefficients , double domain , math_functions::real_matrix * transformation_y_data_cached = nullptr )
{
	PolynomialRegression<STORAGE,STORAGE_ITERATOR> ( y_data , order_of_method , coefficients , transformation_y_data_cached );
}


// void PolynomialRegression ( std::list<double> & x_data , std::list<double> & y_data , unsigned int order_of_method , std::list<double> & coefficients );

template<typename STORAGE , typename STORAGE_ITERATOR>
void PolynomialRegression ( STORAGE & x_data , STORAGE & y_data , unsigned int order_of_method , std::list<double> & coefficients )
{
	if ( ( x_data.size() != y_data.size() ) || ( x_data.size() <= order_of_method ) )
	{
		std::cerr << "PolynomialRegression: Illegal data or illegal order of method" << std::endl;
	}
	else
	{
		unsigned int no_of_data = x_data.size();

		coefficients.clear();

		boost::numeric::ublas::zero_matrix<double> zero1 ( order_of_method + 1 , order_of_method + 1 );
		boost::numeric::ublas::zero_matrix<double> zero2 ( order_of_method + 1 , no_of_data );
		math_functions::real_matrix X_matrix_T ( order_of_method + 1, no_of_data );
		math_functions::real_matrix X_T_X_matrix ( order_of_method + 1, order_of_method + 1 );
		math_functions::real_matrix X_T_X_matrix_invesion ( order_of_method + 1, order_of_method + 1 );
		math_functions::real_matrix transformation_y_data ( order_of_method + 1 , no_of_data );
		math_functions::real_column_vector Y_vector ( no_of_data );
		math_functions::real_column_vector coefficients_vector ( order_of_method + 1 );
		STORAGE_ITERATOR x_data_iter = x_data.begin();
		STORAGE_ITERATOR y_data_iter = y_data.begin();
		X_T_X_matrix = zero1;
		X_T_X_matrix_invesion = zero1;
		X_matrix_T = zero2;

		for ( unsigned int i = 0 ; i < no_of_data ; i ++ )
		{
			for ( unsigned int j = 0 ; j <= order_of_method ; j ++ )
			{
				double x = * x_data_iter;
				X_matrix_T ( j , i ) = pow ( x , ( double ) j );

				for ( unsigned int k = 0 ; k <= order_of_method ; k ++ )
				{
					X_T_X_matrix ( j , k ) += pow ( x , ( double ) j ) * pow ( x , ( double ) k );
				}
			}
			Y_vector ( i ) = ( * y_data_iter );

			++ x_data_iter;
			++ y_data_iter;
		}

		math_tools::Inversion ( X_T_X_matrix_invesion , X_T_X_matrix );
		boost::numeric::ublas::axpy_prod ( X_T_X_matrix_invesion , X_matrix_T , transformation_y_data , true );
		boost::numeric::ublas::axpy_prod ( transformation_y_data , Y_vector, coefficients_vector , true);

		for ( unsigned int i = 0 ; i <= order_of_method ; i ++ )
		{
			coefficients.push_back ( coefficients_vector ( i ) );
		}
	}
}

template<typename STORAGE , typename STORAGE_ITERATOR>
void PolynomialRegression ( STORAGE & x_data , STORAGE & y_data , unsigned int order_of_method , std::list<double> & coefficients , double domain )
{
	PolynomialRegression<STORAGE,STORAGE_ITERATOR> ( x_data , y_data , order_of_method , coefficients );
}

// void HarmonicRegression ( std::list<double> & y_data , unsigned int order_of_method , std::list<double> & coefficients , double domain );

// void HarmonicRegression ( std::list<double> & y_data , unsigned int order_of_method , std::list<double> & coefficients )
// {
// 	HarmonicRegression( y_data , order_of_method , coefficients , y_data.size () );
// }

template<typename STORAGE , typename STORAGE_ITERATOR>
void HarmonicRegression ( STORAGE & y_data , unsigned int order_of_method , std::list<double> & coefficients , double domain )
{
	if ( ( y_data.size() <= order_of_method ) )
	{
		std::cerr << "HarmonicRegresion: Illegal data or illegal order of method" << std::endl;
	}
	else
	{
		unsigned int no_of_data = y_data.size();

		coefficients.clear();

#ifdef DEBUG
		double minimum_x = 0;
		double maximum_x = no_of_data - 1;
#endif

		boost::numeric::ublas::zero_matrix<double> zero1 ( order_of_method + 1 , order_of_method + 1 );
		boost::numeric::ublas::zero_matrix<double> zero2 ( order_of_method + 1 , no_of_data );
		math_functions::real_matrix X_matrix_T ( order_of_method + 1, no_of_data );
		math_functions::real_matrix X_T_X_matrix ( order_of_method + 1, order_of_method + 1 );
		math_functions::real_matrix X_T_X_matrix_invesion ( order_of_method + 1, order_of_method + 1 );
		math_functions::real_matrix transformation_y_data ( order_of_method  + 1, no_of_data );
		math_functions::real_column_vector Y_vector ( no_of_data );
		math_functions::real_column_vector coefficients_vector ( order_of_method  + 1);
		STORAGE_ITERATOR y_data_iter = y_data.begin();
		X_T_X_matrix = zero1;
		X_T_X_matrix_invesion = zero1;
		X_matrix_T = zero2;

		for ( unsigned int i = 0 ; i < no_of_data ; i ++ )
		{
			for ( unsigned int j = 0 ; j <= order_of_method ; j ++ )
			{
				double x = i;
				X_matrix_T ( j , i ) = HarmonicFunction ( x , j , domain );

				for ( unsigned int k = 0 ; k <= order_of_method ; k ++ )
				{
					X_T_X_matrix ( j , k ) += HarmonicFunction ( x , j , domain ) * HarmonicFunction ( x , k , domain );
				}
			}
			Y_vector ( i ) = ( * y_data_iter );

			++ y_data_iter;
		}

// 		math_tools::Inversion ( X_T_X_matrix_invesion , X_T_X_matrix );
// 		boost::numeric::ublas::axpy_prod ( X_T_X_matrix_invesion , X_matrix_T , transformation_y_data , true );
		boost::numeric::ublas::axpy_prod ( X_matrix_T , Y_vector, coefficients_vector , true);

		for ( unsigned int i = 0 ; i <= order_of_method ; i ++ )
		{
			coefficients.push_back ( coefficients_vector ( i ) );
		}

#ifdef DEBUG
// 		std::cout << X_matrix_T << std::endl;
// 		std::cout << X_T_X_matrix << std::endl;
// 		std::cout << X_T_X_matrix_invesion << std::endl;
// 		std::cout << coefficients_vector << std::endl;
#endif
// 		coefficients.clear();
//
// 		for ( unsigned int j = 0 ; j <= order_of_method ; j ++ )
// 		{
// 			double res = 0;
// 			STORAGE_ITERATOR y_data_iter = y_data.begin();
// 			for ( unsigned int i = 0 ; i < no_of_data ; i ++ )
// 			{
// 				res += ( * y_data_iter ) * HarmonicFunction ( i , j , domain );
// 				y_data_iter ++;
// 			}
// 			coefficients.push_back ( res );
// #ifdef DEBUG
// 			std::cout << coefficients.back ( ) << " ";
// #endif
// 		}
// #ifdef DEBUG
// 		std::cout << std::endl;
// #endif
	}
}

// void HarmonicRegression ( std::list<double> & x_data , std::list<double> & y_data , unsigned int order_of_method , std::list<double> & coefficients , double domain );

// inline void HarmonicRegression ( std::list<double> & x_data , std::list<double> & y_data , unsigned int order_of_method , std::list<double> & coefficients )
// {
// 	HarmonicRegression ( x_data , y_data , order_of_method , coefficients , y_data.size () );
// }

template<typename STORAGE , typename STORAGE_ITERATOR>
void HarmonicRegression ( STORAGE & x_data , STORAGE & y_data , unsigned int order_of_method , std::list<double> & coefficients , double domain )
{
	if ( ( x_data.size() != y_data.size() ) || ( x_data.size() <= order_of_method ) )
	{
		std::cerr << "HarmonicRegresion: Illegal data or illegal order of method" << std::endl;
	}
	else
	{
		unsigned int no_of_data = x_data.size();

		coefficients.clear();

		double minimum_x = x_data.front();
		double maximum_x = x_data.back();

		boost::numeric::ublas::zero_matrix<double> zero1 ( order_of_method + 1 , order_of_method + 1 );
		boost::numeric::ublas::zero_matrix<double> zero2 ( order_of_method + 1 , no_of_data );
		math_functions::real_matrix X_matrix_T ( order_of_method + 1, no_of_data );
		math_functions::real_matrix X_T_X_matrix ( order_of_method + 1, order_of_method + 1 );
		math_functions::real_matrix X_T_X_matrix_invesion ( order_of_method + 1, order_of_method + 1 );
		math_functions::real_matrix transformation_y_data ( order_of_method + 1 , no_of_data );
		math_functions::real_column_vector Y_vector ( no_of_data );
		math_functions::real_column_vector coefficients_vector ( order_of_method + 1 );
		STORAGE_ITERATOR x_data_iter = x_data.begin();
		STORAGE_ITERATOR y_data_iter = y_data.begin();
		X_T_X_matrix = zero1;
		X_T_X_matrix_invesion = zero1;
		X_matrix_T = zero2;

		for ( unsigned int i = 0 ; i < no_of_data ; i ++ )
		{
			for ( unsigned int j = 0 ; j <= order_of_method ; j ++ )
			{
				double x = * x_data_iter;
				X_matrix_T ( j , i ) = HarmonicFunction ( x , j , domain );

				for ( unsigned int k = 0 ; k <= order_of_method ; k ++ )
				{
					X_T_X_matrix ( j , k ) += HarmonicFunction ( x , j , domain ) * HarmonicFunction ( x , k , domain );
				}
			}
			Y_vector ( i ) = ( * y_data_iter );

			++ x_data_iter;
			++ y_data_iter;
		}

		math_tools::Inversion ( X_T_X_matrix_invesion , X_T_X_matrix );
// 		boost::numeric::ublas::axpy_prod ( X_T_X_matrix_invesion , X_matrix_T , transformation_y_data , true );
		boost::numeric::ublas::axpy_prod ( X_matrix_T , Y_vector, coefficients_vector , true);

		for ( unsigned int i = 0 ; i <= order_of_method ; i ++ )
		{
			coefficients.push_back ( coefficients_vector ( i ) );
		}

#ifdef DEBUG
		std::cout << X_matrix_T << std::endl;
		std::cout << X_T_X_matrix << std::endl;
		std::cout << X_T_X_matrix_invesion << std::endl;
		std::cout << coefficients_vector << std::endl;
#endif
// 		coefficients.clear();
//
// 		for ( unsigned int j = 0 ; j <= order_of_method ; j ++ )
// 		{
// 			double res = 0;
// 			STORAGE_ITERATOR x_data_iter = x_data.begin();
// 			STORAGE_ITERATOR y_data_iter = y_data.begin();
// 			for ( unsigned int i = 0 ; i < no_of_data ; i ++ )
// 			{
// 				res += ( * y_data_iter ) * HarmonicFunction ( * x_data_iter , j , domain );
// 				x_data_iter ++;
// 				y_data_++ iter;
// 			}
// 			coefficients.push_back ( res );
// #ifdef DEBUG
// 			std::cout << coefficients.back ( ) << " ";
// #endif
// 		}
// #ifdef DEBUG
// 		std::cout << std::endl;
// #endif
	}
}

double Polynomial ( double x , std::list<double> & coefficients );

inline double Polynomial ( double x , std::list<double> & coefficients , double domain )
{
	return Polynomial ( x , coefficients );
}

template<class INTERMEDIATE_TYPE>
INTERMEDIATE_TYPE PolynomialTemplate ( INTERMEDIATE_TYPE x , std::list<INTERMEDIATE_TYPE> & coefficients , INTERMEDIATE_TYPE domain )
{
	INTERMEDIATE_TYPE return_value = 0;
	typename std::list<INTERMEDIATE_TYPE>::iterator iter_end = coefficients.end();
	INTERMEDIATE_TYPE n = 0;
	for ( typename std::list<INTERMEDIATE_TYPE>::iterator iter = coefficients.begin() ; iter != iter_end ; ++ iter )
	{
		return_value += ( * iter ) * math_functions::multiprecision::power<INTERMEDIATE_TYPE> ( x , n ) ;
		n += 1;
	}

	return return_value;
}

double HarmonicSerie ( double x , std::list<double> & coefficients , double domain );

}
#endif
