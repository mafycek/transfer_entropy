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

#include <unordered_map>
#include <list>
#include <tuple>
#include <limits>
#include <cmath>
#include <boost/thread/shared_mutex.hpp>

#include "regression.h"

#include <gsl/gsl_sf_legendre.h>

namespace statistical_functions
{

double HarmonicFunction ( double x , unsigned int order , double domain )
{
	if ( order == 0 )
	{
		return sqrt ( 1 / domain );
	}
	else
	{
		if ( order % 2 == 1 )
		{
			return sqrt ( 2 / domain ) * cos ( 2 * M_PI * x * order / domain );
		}
		else
		{
			return sqrt ( 2 / domain ) * sin ( 2 * M_PI * x * order / domain );
		}
	}
// 	return ( order != 0 ) ? sqrt ( 2 / domain ) * cos ( M_PI * x * order / domain ) : sqrt ( 1 / domain );
}

double LegendreFunction ( double x , unsigned int order , double domain )
{
	return gsl_sf_legendre_Pl ( order , 2 * x / domain - 1 );
// 	return ( order != 0 ) ? sqrt ( 2 / domain ) * cos ( M_PI * x * order / domain ) : sqrt ( 1 / domain );
}

double Polynomial ( double x , std::list<double> & coefficients )
{
	double return_value = 0;
	std::list<double>::iterator iter_end = coefficients.end();
	double n = 0;
	for ( std::list<double>::iterator iter = coefficients.begin() ; iter != iter_end ; ++ iter )
	{
		return_value += ( * iter ) * pow ( x , n ) ;
		n += 1;
	}

	return return_value;
}

double HarmonicSerie ( double x , std::list<double> & coefficients , double domain )
{
	double return_value = 0;
	std::list<double>::iterator iter_end = coefficients.end();
	double n = 0;
	for ( std::list<double>::iterator iter = coefficients.begin() ; iter != iter_end ; ++ iter )
	{
		return_value += ( * iter ) * HarmonicFunction ( x , n , domain );
		n += 1;
	}

	return return_value;
}

double LegendreSerie ( double x , std::list<double> & coefficients , double domain )
{
	double return_value = 0;
	std::list<double>::iterator iter_end = coefficients.end();
	double n = 0;
	for ( std::list<double>::iterator iter = coefficients.begin() ; iter != iter_end ; ++ iter )
	{
		return_value += ( * iter ) * LegendreFunction ( x , n , domain );
		n += 1;
	}

	return return_value;
}


class default_double
{
	private:
		double _value;

	public:
		default_double ()
		: _value ( std::numeric_limits<double>::quiet_NaN() )
		{
		}

		double & operator() ()
		{
			return _value;
		}
};


inline double GetPowerMapping ( std::unordered_map<std::tuple< unsigned int , unsigned int> , double , KeyHasherMatricesForRegression> & power , unsigned int position , unsigned int exponent
#ifdef THREADING
, boost::shared_mutex & power_mutex 
#endif
)
{
	std::tuple< unsigned int , unsigned int> index = std::make_tuple ( position , exponent );
#ifdef THREADING
	power_mutex.lock_shared();
#endif
	typename std::unordered_map<std::tuple< unsigned int , unsigned int> , double , KeyHasherMatricesForRegression>::iterator iterator_power;
	iterator_power = power.find ( index );

	if ( iterator_power == power.end() )
	{
#ifdef THREADING
		power_mutex.unlock_shared();
		power_mutex.lock();
#endif
		double result = pow ( position , exponent );
		power [ index ] = result;
#ifdef THREADING
		power_mutex.unlock();
#endif
		return result;
	}
	else
	{
		double result = ( * iterator_power ) . second;
#ifdef THREADING
		power_mutex.unlock_shared();
#endif
		return result;
	}
}

// TODO: prepare the function ready for other regression methods
void PrepareMatricesForRegression ( unsigned int no_of_data , unsigned int order_of_method , math_functions::real_matrix & transformation_y_data_cached  )
{
	static std::unordered_map<std::tuple< unsigned int , unsigned int> , double , KeyHasherMatricesForRegression> powers;
#ifdef THREADING
	static boost::shared_mutex power_mutex;
#endif

	transformation_y_data_cached.resize(order_of_method + 1 , no_of_data);
	math_functions::real_matrix X_T_X ( order_of_method  + 1, order_of_method  + 1 );
	math_functions::real_matrix X_T_X_inversion_cached ( order_of_method + 1, order_of_method + 1 );
	boost::numeric::ublas::zero_matrix<double> zero1 ( order_of_method + 1 , order_of_method + 1 );
	boost::numeric::ublas::zero_matrix<double> zero2 ( order_of_method + 1 , no_of_data );
	math_functions::real_matrix X_matrix_T ( order_of_method  + 1, no_of_data );
	X_T_X = zero1;
	X_matrix_T = zero2;

	for ( unsigned int i = 0 ; i < no_of_data ; i ++ )
	{
		for ( unsigned int j = 0 ; j <= order_of_method ; j ++ )
		{
			double power_1 = GetPowerMapping ( powers , i , j
#ifdef THREADING
			, power_mutex
#endif
			);

			X_matrix_T ( j , i ) = power_1;
// 			pow ( x , ( double ) j );

			for ( unsigned int k = 0 ; k <= order_of_method ; k ++ )
			{
				double power_2 = GetPowerMapping ( powers , i , k
#ifdef THREADING
				, power_mutex
#endif
				);

				double result = power_1 * power_2;
// 				pow ( x , ( double ) j ) * pow ( x , ( double ) k );
				X_T_X ( j , k ) += result;
			}
		}
	}

	math_tools::Inversion ( X_T_X_inversion_cached , X_T_X );
	boost::numeric::ublas::axpy_prod ( X_T_X_inversion_cached , X_matrix_T , transformation_y_data_cached , true );
#ifdef DEBUG
	std::cout << X_T_X << std::endl;
	std::cout << X_T_X_inversion_cached << std::endl;
	std::cout << X_matrix_T << std::endl;
	std::cout << transformation_y_data_cached << std::endl;
#endif
}

}
