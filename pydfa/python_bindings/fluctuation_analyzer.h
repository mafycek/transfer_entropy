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

#ifndef FLUCTUATION_ANALYZER
#define FLUCTUATION_ANALYZER

#include <list>
#include <vector>
#include <string>
#include <fstream>
#include <tuple>
#include <map>
#include <cmath> 

#include <boost/concept_check.hpp>
#include <boost/thread/pthread/thread_data.hpp>
#include <boost/thread/detail/thread.hpp>

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/hypot.hpp>

#include <boost/numeric/ublas/hermitian.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "math_functions.h"
#include "regression.h"
#include "debug_support.h"

namespace statistical_functions
{

void enqueue_thread_for_cleanup ( );

void IntegrationOfDataset ( std::list<double> & data );

class SegmentWorker
{
protected:
	std::list<double> & _data;
	std::vector<double> & _q_th_order_fluctuation;
	
	std::ofstream & _output_stream;

	unsigned int _segment_size;
	bool _moving_segment_adjustment;
	unsigned int _moving_segment_shift;
	bool _sample_variance;

	SegmentWorker ( std::list<double> & data , std::vector<double> & q_th_order_fluctuation , std::ofstream & output_stream , unsigned int segment_size , bool moving_segment_adjustment = false , unsigned int moving_segment_shift = 2 , bool sample_variance = false )
	: _data ( data ) , _q_th_order_fluctuation ( q_th_order_fluctuation ) , _output_stream ( output_stream ) , _segment_size ( segment_size ) , _moving_segment_adjustment ( moving_segment_adjustment ) , _moving_segment_shift ( moving_segment_shift ) , _sample_variance ( sample_variance )
	{
	}

};


void DetrendedFluctuationAnalysisExponentiallyScaledDomain ( std::list<double> & data , std::vector<unsigned int> & segment_sizes , std::map<std::tuple<double, unsigned int, unsigned int> , double > & DFAvariance , std::vector<double> & q_th_order_fluctuation , unsigned int start_index , unsigned int end_index , unsigned segments , const std::vector<long> & order_of_method , std::string output_analysis , bool moving_segment_adjustment , unsigned int moving_segment_shift , bool sample_variance = false , bool double_integration = false
#ifdef THREADING
, unsigned int number_threads = 2
#ifdef MICROTHREADING
, unsigned int microthreads = 2
#endif
#endif
, void ( * fluctuation_analyzer ) ( std::vector<double> &  , unsigned int , std::list<double> &  , double , math_functions::real_matrix * ) = &PolynomialRegression<std::vector<double>,std::vector<double>::iterator> , double ( * polynomial ) ( double  , std::list<double> & , double  ) = &Polynomial );

void DetrendedFluctuationAnalysisExponentiallyScaledDomainSampleVarianceTest ( std::list<double> & data , std::list<unsigned int> & segment_sizes , std::map<std::tuple<double, unsigned int, unsigned int> , double > & DFAvariance , std::vector<double> & q_th_order_fluctuation , unsigned int start_index , unsigned int end_index , unsigned segments , unsigned int order_of_method , std::string output_analysis , bool moving_segment_adjustment , unsigned int moving_segment_shift , bool sample_variance = false , bool double_integration = false
#ifdef THREADING
, unsigned int number_threads = 2
#ifdef MICROTHREADING
, unsigned int microthreads = 2
#endif
#endif
, void ( * fluctuation_analyzer ) ( std::vector<double> &  , unsigned int , std::list<double> &  , double , math_functions::real_matrix * ) = &PolynomialRegression<std::vector<double>,std::vector<double>::iterator> , double ( * polynomial ) ( double  , std::list<double> & , double  ) = &Polynomial );


void DetrendedFluctuationAnalysis ( std::list<double> & data , std::vector<unsigned int> & segment_sizes , std::map<std::tuple<double, unsigned int, unsigned int> , double > & DFAvariance , std::vector<double> & q_th_order_fluctuation , const std::vector<long> & order_of_method , std::string output_analysis , bool moving_segment_adjustment , unsigned int moving_segment_shift , bool sample_variance , bool double_integration
#ifdef THREADING
, unsigned int number_threads = 2
#ifdef MICROTHREADING
, unsigned int microthreads = 2
#endif
#endif
, void ( * fluctuation_analyzer ) ( std::vector<double> &  , unsigned int , std::list<double> &  , double , math_functions::real_matrix * ) = &PolynomialRegression<std::vector<double>,std::vector<double>::iterator> , double ( * polynomial ) ( double  , std::list<double> & , double  ) = &Polynomial );

void DetrendedCrossCorrelationAnalysis ( std::list<double> & data1 , std::list<double> & data2 , std::vector<unsigned int> & segment_sizes , std::map<std::tuple<double, unsigned int, unsigned int> , std::array<double,4> > & DCCA , std::vector<double> & q_th_order_fluctuation , std::vector<unsigned int> & order_of_method , std::string output_analysis , bool moving_segment_adjustment , unsigned int moving_segment_shift , bool sample_variance , bool double_integration
#ifdef THREADING
, unsigned int number_threads = 2
#endif
, void ( * fluctuation_analyzer ) ( std::vector<double> &  , unsigned int , std::list<double> &  , double , math_functions::real_matrix * ) = &PolynomialRegression<std::vector<double>,std::vector<double>::iterator> , double ( * polynomial ) ( double  , std::list<double> & , double  ) = &Polynomial );

void DetrendedFluctuationAnalysisLinearlyScaledDomain ( std::list<double> & data , std::vector<unsigned int> & segment_sizes , std::map<std::tuple<double, unsigned int, unsigned int> , double > & DFAvariance , std::vector<double> & q_th_order_fluctuation , unsigned int start_index , unsigned int end_index , unsigned segment_step = 1 , const std::vector<long> & order_of_method = std::vector<long>().operator=(std::vector<long>()) , std::string output_analysis = std::string ( "" ) , bool moving_segment_adjustment = false , unsigned int moving_segment_shift = 2 , bool sample_variance = false , bool double_integration = false
#ifdef THREADING
, unsigned int number_threads = 2
#ifdef MICROTHREADING
, unsigned int microthreads = 2
#endif
#endif
, void ( * fluctuation_analyzer ) ( std::vector<double> &  , unsigned int , std::list<double> & , double , math_functions::real_matrix * transformation_y_data_cached ) = &PolynomialRegression<std::vector<double>,std::vector<double>::iterator> , double ( * polynomial ) ( double  , std::list<double> & , double  ) = &Polynomial
);


void DetrendedHarmonicFluctuationAnalysis ( std::list<double> & data , std::list<double> & DFAvariance , double q_th_order_fluctuation = 1.0 , unsigned int order_of_method = 1 , unsigned segment_step = 1 , unsigned int min_number_of_segments = 2 , std::string output_analysis = std::string ( "" ) );

void Autocorrelation ( std::vector<double> data , std::list<double> & autocorrelation , unsigned segment_step = 1 , unsigned int min_number_of_segments = 2 );

void AutocorrelationCyclicData ( std::vector<double> data , std::list<double> & autocorrelation , unsigned segment_step = 1 , unsigned int  min_number_of_segments = 2 );

void StructureFunction ( std::vector<double> data , std::list<double> & timelagfunction , double q_th_order = 1 , unsigned segment_step = 1 , unsigned int min_number_of_segments = 2 );

void StructureFunctionCyclicData ( std::vector<double> data , std::list<double> & timelagfunction , double q_th_order = 1 , unsigned segment_step = 1 , unsigned int min_number_of_segments = 2 );

struct extremum
{
	int _position;
	double _value;

	extremum ( int position , double value )
	: _position( position ) , _value ( value )
	{
	}
};

bool IsLocalMaximum ( double y1 , double y2 , double y3 );

bool IsLocalMinimum ( double y1 , double y2 , double y3 );

template <typename STORAGE_ITERATOR>
bool IsLocalMaximum ( STORAGE_ITERATOR y1 , STORAGE_ITERATOR y2 , STORAGE_ITERATOR y3 )
{
	if ( ( ( ( * y2 ) > ( * y1 ) ) || ( ( ( * y2 ) == ( * y1 ) ) && ( ( * y2 ) > ( * ( -- y1 ) ) ) ) ) && ( ( ( * y2 ) > ( * y3 ) ) || ( ( ( * y2 ) == ( * y3 ) ) && ( ( * y2 ) > ( * ( ++ y3 ) ) ) ) ) )
	{
		return true;
	}
	else
	{
		return false;
	}
}

template <typename STORAGE_ITERATOR>
bool IsLocalMinimum ( STORAGE_ITERATOR y1 , STORAGE_ITERATOR y2 , STORAGE_ITERATOR y3 )
{
	if ( ( ( ( * y2 ) < ( * y1 ) ) || ( ( ( * y2 ) == ( * y1 ) ) && ( ( * y2 ) < ( * ( -- y1 ) ) ) ) ) && ( ( ( * y2 ) < ( * y3 ) ) || ( ( ( * y2 ) == ( * y3 ) ) && ( ( * y2 ) < ( * ( ++ y3 ) ) ) ) ) )
	{
		return true;
	}
	else
	{
		return false;
	}
}


template<typename STORAGE , typename STORAGE_ITERATOR>
bool SiftingOfDataset ( STORAGE & data , STORAGE_ITERATOR data_begin , STORAGE_ITERATOR data_end , STORAGE & sifted_dataset , STORAGE & fluctuation_part_of_dataset , STORAGE_ITERATOR & fluctuation_part_of_dataset_begin , STORAGE_ITERATOR & fluctuation_part_of_dataset_end , bool extended_algorithm = false , bool extended_algorithm_shift_limits = false , bool accept_sifting = false , double level_of_difference = 1e-6 )
{
// searching for maxima and minima of the dataset
	std::list<extremum *> maxima;
	std::list<extremum *> minima;
	STORAGE_ITERATOR fluctuation_part_segment_begin;
	STORAGE_ITERATOR fluctuation_part_segment_end;
	int fluctuation_part_segment_begin_coordinate = 0;
	int fluctuation_part_segment_end_coordinate = 0;

// possible enhancement regarding maximum and minimum from part ahead of the segment
	if ( ( extended_algorithm == false ) || ( data_begin == data . begin ( ) ) )
	{
		minima . push_back ( new extremum ( 0 , ( * data_begin ) ) );
#ifdef DEBUG
		std::cout << "Additional minimum: " << 0 << " " << ( * data_begin ) << std::endl;
#endif
		maxima . push_back ( new extremum ( 0 , ( * data_begin ) ) );
#ifdef DEBUG
		std::cout << "Additional maximum: " << 0 << " " << ( * data_begin ) << std::endl;
#endif
		fluctuation_part_segment_begin = data_begin;
	}
	else
	{
		STORAGE_ITERATOR before_iter = data_begin;
		before_iter --;
		STORAGE_ITERATOR iter = data_begin;
		STORAGE_ITERATOR after_iter = data_begin;
		++ after_iter;
		int datapoint = 0;
		bool minimum_found = false;
		bool maximum_found = false;

		for ( ; ( before_iter != data.begin() ) && ( ( minimum_found == false ) || ( maximum_found == false ) ) ; )
		{
			if ( IsLocalMinimum ( ( before_iter ) , ( iter ) , ( after_iter ) ) )
			{
				minima . push_back ( new extremum ( datapoint , ( * iter ) ) );
				minimum_found = true;
#ifdef DEBUG
				std::cout << "Additional minimum: " << datapoint  << " " << ( * iter ) << std::endl;
#endif
				if ( maximum_found == true )
				{
					maxima . push_front ( new extremum ( datapoint , ( * iter ) ) );
#ifdef DEBUG
					std::cout << "Additional maximum: " << datapoint  << " " << ( * iter ) << std::endl;
#endif
				}
			}

			if ( IsLocalMaximum ( ( before_iter ) , ( iter ) , ( after_iter ) ) )
			{
				maxima . push_back ( new extremum ( datapoint , ( * iter ) ) );
				maximum_found = true;
#ifdef DEBUG
				std::cout << "Additional maximum: " << datapoint  << " " << ( * iter ) << std::endl;
#endif
				if ( minimum_found == true )
				{
					minima . push_front ( new extremum ( datapoint , ( * iter ) ) );
#ifdef DEBUG
					std::cout << "Additional minimum: " << datapoint  << " " << ( * iter ) << std::endl;
#endif
				}
			}

			after_iter = iter;
			iter = before_iter;
			before_iter --;
			datapoint --;
		}

		if ( before_iter == data.begin() )
		{
			minima . push_back ( new extremum ( 0 , ( * data_begin ) ) );
#ifdef DEBUG
			std::cout << "Minimum: " << 0 << " " << ( * data_begin ) << std::endl;
#endif
		}

		if ( before_iter == data.begin() )
		{
			maxima . push_back ( new extremum ( 0 , ( * data_begin ) ) );
#ifdef DEBUG
			std::cout << "Maximum: " << 0 << " " << ( * data_begin ) << std::endl;
#endif
		}

		if ( extended_algorithm_shift_limits )
		{
			fluctuation_part_segment_begin = after_iter;
			fluctuation_part_segment_begin_coordinate = datapoint + 1;
		}
	}

	unsigned int datapoint = 1;
	STORAGE_ITERATOR end_iter = data_end;
	STORAGE_ITERATOR iter = data_begin;
	++ iter;
	STORAGE_ITERATOR before_iter = data_begin;
	STORAGE_ITERATOR after_iter = iter;
	++ after_iter;

// searching for local extrema in the interval
	for ( ; after_iter != end_iter ; )
	{
		double actual_datapoint = ( * iter );

		if ( IsLocalMinimum ( ( before_iter ) , ( iter ) , ( after_iter ) ) )
		{
			minima . push_back ( new extremum ( datapoint , ( * iter ) ) );
#ifdef DEBUG
			std::cout << "Minimum: " << datapoint  << " " << ( * iter ) << std::endl;
#endif
		}
		else
		{
			if ( IsLocalMaximum ( ( before_iter ) , ( iter ) , ( after_iter ) ) )
			{
				maxima . push_back ( new extremum ( datapoint , ( * iter ) ) );
#ifdef DEBUG
				std::cout << "Maximum: " << datapoint << " " << ( * iter ) << std::endl;
#endif
			}
		}

		++ datapoint;
		before_iter = iter;
		iter = after_iter;
		++ after_iter;
	}

	unsigned int size_of_segment = datapoint + 1;

// looking for extrema in after the segment
	if ( ( extended_algorithm == false ) || ( data_end == data . end ( ) ) )
	{
		STORAGE_ITERATOR temp_iter = data_end;
		temp_iter --;
		double actual_datapoint = ( * temp_iter );
		minima . push_back ( new extremum ( datapoint , actual_datapoint ) );
#ifdef DEBUG
		std::cout << "Additional minimum: " << datapoint  << " " << actual_datapoint << std::endl;
#endif
		maxima . push_back ( new extremum ( datapoint , actual_datapoint ) );
#ifdef DEBUG
		std::cout << "Additional maximum: " << datapoint << " " << actual_datapoint << std::endl;
#endif
		fluctuation_part_segment_end = after_iter;
		fluctuation_part_segment_end_coordinate = datapoint;
	}
	else
	{
		int temp_datapoint = datapoint;
		iter = data_end;
		iter --;
		before_iter = data_end;
		before_iter --;
		before_iter --;
		after_iter = data_end;
		bool minimum_found = false;
		bool maximum_found = false;

		for ( ; ( after_iter != data.end() ) && ( ( minimum_found == false ) || ( maximum_found == false ) ); )
		{
			if ( IsLocalMinimum ( ( before_iter ) , ( iter ) , ( after_iter ) ) )
			{
				minima . push_back ( new extremum ( temp_datapoint , ( * iter ) ) );
				minimum_found = true;
#ifdef DEBUG
				std::cout << "Additional minimum: " << temp_datapoint  << " " << ( * iter ) << std::endl;
#endif
				if ( maximum_found == true )
				{
					maxima . push_back ( new extremum ( temp_datapoint , ( * iter ) ) );
#ifdef DEBUG
					std::cout << "Additional maximum: " << temp_datapoint  << " " << ( * iter ) << std::endl;
#endif
				}
			}

			if ( IsLocalMaximum ( ( before_iter ) , ( iter ) , ( after_iter ) ) )
			{
				maxima . push_back ( new extremum ( temp_datapoint , ( * iter ) ) );
				maximum_found = true;
#ifdef DEBUG
				std::cout << "Additional maximum: " << temp_datapoint  << " " << ( * iter ) << std::endl;
#endif
				if ( minimum_found == true )
				{
					minima . push_back ( new extremum ( temp_datapoint , ( * iter ) ) );
#ifdef DEBUG
					std::cout << "Additional minimum: " << temp_datapoint  << " " << ( * iter ) << std::endl;
#endif
				}
			}

			before_iter = iter;
			iter = after_iter;
			++ after_iter;
			++ temp_datapoint;
		}

		if ( after_iter == data.end() )
		{
			minima . push_back ( new extremum ( temp_datapoint , ( * iter ) ) );
#ifdef DEBUG
			std::cout << "Additional minimum: " << temp_datapoint  << " " << ( * iter ) << std::endl;
#endif
		}

		if ( after_iter == data.end() )
		{
			maxima . push_back ( new extremum ( temp_datapoint , ( * iter ) ) );
#ifdef DEBUG
			std::cout << "Additional maximum: " << temp_datapoint  << " " << ( * iter ) << std::endl;
#endif
		}

		if ( extended_algorithm_shift_limits )
		{
			fluctuation_part_segment_end = iter;
			fluctuation_part_segment_end_coordinate = temp_datapoint;
		}
	}

	if ( extended_algorithm_shift_limits )
	{
		fluctuation_part_of_dataset.clear ();
// copy data in neighborhood of segment
		fluctuation_part_of_dataset = STORAGE ( fluctuation_part_segment_begin , fluctuation_part_segment_end );

// setting up limits of segment within its neighborhood
		int count = fluctuation_part_segment_begin_coordinate;
		for ( auto iter = fluctuation_part_of_dataset.begin() ; iter != fluctuation_part_of_dataset.end() ; ++ iter , count ++ )
		{
			if ( count == 0 )
			{
				fluctuation_part_of_dataset_begin = iter;
			}

			if ( count == size_of_segment - 1 )
			{
				fluctuation_part_of_dataset_end = iter;
			}

			sifted_dataset . push_back ( 0 );
		}
	}

	int number_of_extrema = 0;

// interpolation of maxima and minima
// initialization of the library
	gsl_interp_accel * acc_maxima = gsl_interp_accel_alloc ();
	gsl_interp_accel * acc_minima = gsl_interp_accel_alloc ();
	const gsl_interp_type * type_of_interpolation = gsl_interp_cspline;

// condition that prevents GSL crash due to small number of maxima or minima
	if ( ( maxima.size() <= 2 ) || ( minima.size() <= 2 ) )
	{
		std::cerr << "Unable to execute sifting. Not enough datapoints. Actually we have " << maxima.size() << " maxima and " << minima.size() << " minima." << std::endl;

#ifdef DEBUG
		STORAGE_ITERATOR data_iter = data_begin;
		for ( unsigned int xi = 0 ; xi <= datapoint ; xi += 1 , data_iter ++ )
		{
			std::cout << xi << " datapoint: " << ( * data_iter ) << std::endl;
		}
#endif
	}

	gsl_spline * spline_maxima = gsl_spline_alloc ( type_of_interpolation , maxima.size() );
	gsl_spline * spline_minima = gsl_spline_alloc ( type_of_interpolation , minima.size() );

// preparation of data
	double * positions_maxima = ( double * ) malloc ( maxima.size() * sizeof ( double ) );
	double * dataset_maxima = ( double * ) malloc ( maxima.size() * sizeof ( double ) );

	std::list<extremum *>::iterator iter2 = maxima.begin();
	for ( unsigned int i = 0 ; iter2 != maxima . end () ; i ++ )
	{
		positions_maxima [ i ] = ( * iter2 ) ->_position;
		dataset_maxima [ i ] = ( * iter2 ) ->_value;

		if ( ( ( positions_maxima [ i ] > 0 ) || ( ( positions_maxima [ i ] == 0 ) && ( i > 0 ) ) ) && ( positions_maxima [ i ] < datapoint ) )
		{
			number_of_extrema ++;
		}
		++ iter2;
	}

	double * positions_minima = ( double * ) malloc ( minima.size() * sizeof ( double ) );
	double * dataset_minima = ( double * ) malloc ( minima.size() * sizeof ( double ) );

	std::list<extremum *>::iterator iter3 = minima . begin ( );
	for ( unsigned int i = 0 ; i < minima . size() ; ++ i )
	{
		positions_minima [ i ] = ( * iter3 ) ->_position;
		dataset_minima [ i ] = ( * iter3 ) ->_value;

		if ( ( ( positions_minima [ i ] > 0 ) || ( ( positions_minima [ i ] == 0 ) && ( i > 0 ) ) ) && ( positions_minima [ i ] < datapoint ) )
		{
			number_of_extrema ++;
		}
		++ iter3;
	}

// interpolation
	gsl_spline_init ( spline_maxima , positions_maxima , dataset_maxima , maxima . size ( ) );
	gsl_spline_init ( spline_minima , positions_minima , dataset_minima , minima . size ( ) );

// filtration
#ifdef DEBUG
	STORAGE_ITERATOR data_iter;
	STORAGE_ITERATOR sift_iter = sifted_dataset.begin ();
	int xi = fluctuation_part_segment_begin_coordinate;
	for ( data_iter = fluctuation_part_segment_begin ; data_iter != fluctuation_part_segment_end ; data_iter ++  )
	{
// xi < fluctuation_part_segment_end_coordinate
		std::cout << xi << " datapoint: " << ( * data_iter ) << " max: " << gsl_spline_eval (spline_maxima, xi, acc_maxima) << " min: " << gsl_spline_eval (spline_minima, xi, acc_minima) << " mean: " << ( gsl_spline_eval (spline_maxima, xi, acc_maxima) + gsl_spline_eval (spline_minima, xi, acc_minima) ) / 2.0  << " difference: " << ( * data_iter ) - ( gsl_spline_eval (spline_maxima, xi, acc_maxima) + gsl_spline_eval (spline_minima, xi, acc_minima) ) / 2.0 << " sifted function: " << ( * sift_iter ) << std::endl;
		xi ++;
		++ sift_iter;
	}
#endif

// test
	int crossings = 0;
	double maxima_mean = 0;
	double minima_mean = 0;
	STORAGE_ITERATOR iter4 = data_begin;
	bool above = ( ( gsl_spline_eval ( spline_maxima, 0, acc_maxima ) + gsl_spline_eval ( spline_minima, 0, acc_minima ) ) / 2.0 > ( * iter4 ) ) ? false : true;
	iter4 ++;
	for ( unsigned int xi = 1 ; xi <= datapoint ; xi += 1 )
	{
		double maxima_total = gsl_spline_eval ( spline_maxima, xi, acc_maxima );
		double minima_total = gsl_spline_eval ( spline_minima, xi, acc_minima );
		double mean = ( maxima_total + minima_total ) / 2.0;
		double actual_data = ( * iter4 );

		maxima_mean += maxima_total - mean;
		minima_mean += minima_total - mean;
		if ( ( mean > actual_data ) && ( above == true ) )
		{
			crossings ++;
			above = false;
		}
		else
		{
			if ( ( mean < actual_data ) && ( above == false ) )
			{
				crossings ++;
				above = true;
			}
		}
		iter4 ++;
	}

#ifdef DEBUG
	std::cout << "We encouted " << crossings << " crossings and " << number_of_extrema << " local extrema " << std::endl;
	std::cout << "Maxima " << maxima_mean << " minima " << minima_mean << " difference " << maxima_mean + minima_mean << std::endl;
#endif
	bool return_value;

// decision to accept sifting
	if ( accept_sifting || ( ( number_of_extrema + 1 >= crossings ) && ( crossings >= number_of_extrema - 1 ) && ( fabs ( maxima_mean + minima_mean ) <= level_of_difference ) ) )
	{
#ifdef DEBUG
		std::cout << "Sifting accepted" << std::endl;
#endif
		STORAGE_ITERATOR dataset_iter = fluctuation_part_segment_begin;
		STORAGE_ITERATOR sift_end = sifted_dataset.end();
		STORAGE_ITERATOR sift_iter = sifted_dataset.begin();
		STORAGE_ITERATOR fluctuation_iter = fluctuation_part_of_dataset.begin();
		STORAGE_ITERATOR fluctuation_end = fluctuation_part_of_dataset.end();
		int position = fluctuation_part_segment_begin_coordinate;

		for ( ; dataset_iter != fluctuation_part_segment_end ; )
		{
			double approximation = ( gsl_spline_eval (spline_maxima, position, acc_maxima) + gsl_spline_eval (spline_minima, position, acc_minima) ) / 2.0;
			double & sift_data = ( * sift_iter );
			double & fluc_data = ( * fluctuation_iter );
			double & original_data = ( * dataset_iter );

			if ( sift_iter != sift_end )
			{
				( * sift_iter ) += approximation;
			}

			if ( fluctuation_iter != fluctuation_end )
			{
				( * fluctuation_iter ) = ( * dataset_iter ) - approximation;
			}

			++ position;
			++ sift_iter;
			++ fluctuation_iter;
			++ dataset_iter;
		}
		return_value = ( ( number_of_extrema + 1 >= crossings ) && ( crossings >= number_of_extrema - 1 ) && ( fabs ( maxima_mean + minima_mean ) <= level_of_difference ) );
	}
	else
	{
#ifdef DEBUG
		std::cout << "Sifting is not accepted" << std::endl;
#endif
		STORAGE_ITERATOR dataset_iter = data_begin;
		STORAGE_ITERATOR sift_end = sifted_dataset.end();
		STORAGE_ITERATOR sift_iter = sifted_dataset.begin();
		STORAGE_ITERATOR fluctuation_iter = fluctuation_part_of_dataset.begin();
		STORAGE_ITERATOR fluctuation_end = fluctuation_part_of_dataset.end();
		int position = 0;
		for ( ; dataset_iter != data_end ; )
		{
			double & sift_data = ( * sift_iter );
			double & fluc_data = ( * fluctuation_iter );
			double & original_data = ( * dataset_iter );
			if ( sift_iter != sift_end)
			{
				( * sift_iter ) += 0;
			}
			if ( fluctuation_iter != fluctuation_end )
			{
				( * fluctuation_iter ) = ( * dataset_iter );
			}

			++ position;
			++ sift_iter;
			++ fluctuation_iter;
			++ dataset_iter;
		}

		return_value = false;
	}

// finishing the library
	free ( positions_maxima );
	free ( dataset_maxima );
	gsl_spline_free ( spline_maxima );
	gsl_spline_free ( spline_minima );
	gsl_interp_accel_free ( acc_maxima);
	gsl_interp_accel_free ( acc_minima );

// clean up of the memory
	std::list<extremum *>::iterator end_iter2 = maxima.end();
	for ( std::list<extremum *>::iterator iter = maxima.begin() ; iter != end_iter2 ; ++ iter )
	{
		extremum * actual_data = ( * iter );
		delete actual_data;
	}

	std::list<extremum *>::iterator end_iter3 = minima.end();
	for ( std::list<extremum *>::iterator iter = minima.begin() ; iter != end_iter3 ; ++ iter )
	{
		extremum * actual_data = ( * iter );
		delete actual_data;
	}

	minima.clear();
	maxima.clear();

	return return_value;
}

void EmpiricalModeDecomposionDetrendedFluctuationAnalysisExponentiallyScaledDomain ( std::list<double> & data , std::list<unsigned int> & segment_sizes , std::vector<std::list<double> > & DFAvariance , std::vector<double> & q_th_order_fluctuation , unsigned int start_index , unsigned int end_index , unsigned segments , std::string output_analysis , bool moving_segment_adjustment , unsigned int moving_segment_shift , bool sample_variance , bool double_integration , int minimal_sifting
#ifdef THREADING
, unsigned int number_threads
#ifdef MICROTHREADING
, unsigned int microthreads
#endif
#endif
, bool ( * sifting_method )  ( std::list<double> & , std::list<double>::iterator , std::list<double>::iterator , std::list<double> & , std::list<double> & , std::list<double>::iterator &, std::list<double>::iterator &, bool , bool , bool , double ) = &SiftingOfDataset<std::list<double>,std::list<double>::iterator> );

void EmpiricalModeDecomposionDetrendedFluctuationAnalysis ( std::list<double> & data , std::list<unsigned int> & segment_sizes , std::vector<std::list<double> > & DFAvariance , std::vector<double> & q_th_order_fluctuation , unsigned int start_index , unsigned int end_index , unsigned segment_step , std::string output_analysis , bool moving_segment_adjustment , unsigned int moving_segment_shift , bool sample_variance , bool double_integration , int minimal_sifting
#ifdef THREADING
, unsigned int number_threads
#ifdef MICROTHREADING
, unsigned int microthreads
#endif
#endif
, bool ( * sifting_method )  ( std::list<double> & , std::list<double>::iterator , std::list<double>::iterator , std::list<double> & , std::list<double> & , std::list<double>::iterator & , std::list<double>::iterator & , bool , bool , bool , double ) );

}

#endif
