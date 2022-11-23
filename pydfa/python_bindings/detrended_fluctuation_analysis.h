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
#include <string>
#include <fstream>
#include <tuple>
#include <map>
#include <cmath> 
#include <algorithm>
#include <chrono>
#include <iostream>
#include <functional>
#include <exception>
#include <stdexcept>
#include <locale>

#include <boost/concept_check.hpp>
#include <boost/thread/pthread/thread_data.hpp>
#include <boost/thread/detail/thread.hpp>

#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/concept_check.hpp>

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

#include <boost/date_time/posix_time/posix_time.hpp>

#include "debug_support.h"
#include "regression.h"

#ifndef DETRENDED_FLUCTUATION_ANALYZER
#define DETRENDED_FLUCTUATION_ANALYZER

namespace statistical_functions
{

namespace detrended_fluctuation_analysis
{

void enqueue_thread_for_cleanup_template ( );

class DFASynchronization
{
	private:
// variables for conditional variables
		boost::mutex _add_jobs_mutex;

		boost::mutex _text_output_mutex;

		boost::condition_variable _add_jobs_CV;

		std::list< boost::thread::id > _finished_jobs;

	public:
		static DFASynchronization & GetDFASynchronization ()
		{
			static DFASynchronization _sychronization;

			return ( _sychronization );
		}

		boost::mutex & GetAddJobsMutex ()
		{
			return _add_jobs_mutex;
		}
		
		boost::condition_variable & GetAddJobsCV()
		{
			return _add_jobs_CV;
		}

		std::list< boost::thread::id > & GetFinishedJobs ()
		{
			return _finished_jobs;
		}

		boost::mutex & GetTextOutputMutex ()
		{
			return _text_output_mutex;
		}

	private:

		DFASynchronization ()
		{
		}

		~DFASynchronization ()
		{
		}

		DFASynchronization ( const DFASynchronization & ) = delete;
		DFASynchronization ( DFASynchronization && ) = delete;
		DFASynchronization & operator= ( const DFASynchronization & );
		DFASynchronization & operator= ( DFASynchronization && );
};


template<typename DATATYPE , template <typename,typename,typename> class MATRIX = boost::numeric::ublas::matrix >
class DFATransformationMatrices
{
	private:
		std::map < unsigned int , MATRIX<DATATYPE, boost::numeric::ublas::row_major , boost::numeric::ublas::unbounded_array<DATATYPE> > > _matrix_storage;

		boost::mutex _access_mutex;

	public:
		boost::mutex & GetAccessMutex ()
		{
			return _access_mutex;
		}

		std::map < unsigned int , MATRIX<DATATYPE, boost::numeric::ublas::row_major , boost::numeric::ublas::unbounded_array<DATATYPE> > > & GetMatrixStorage ()
		{
			return _matrix_storage;
		}

		static DFATransformationMatrices & GetDFATransformationMatrices ()
		{
			static DFATransformationMatrices _transformation_matrices;

			return _transformation_matrices;
		}

	private:

		DFATransformationMatrices ()
		{
		}

		~DFATransformationMatrices ()
		{
		}

		DFATransformationMatrices ( const DFATransformationMatrices & ) = delete;
		DFATransformationMatrices ( DFATransformationMatrices && ) = delete;
		DFATransformationMatrices & operator= ( const DFATransformationMatrices & );
		DFATransformationMatrices & operator= ( DFATransformationMatrices && );
};

template <class INPUT_TYPE>
class SegmentWorkerTemplate
{
protected:

	std::list<INPUT_TYPE> & _data;

	std::vector<double> & _q_th_order_fluctuation;

	std::ofstream & _output_stream;

	unsigned int _segment_size;

	bool _moving_segment_adjustment;

	unsigned int _moving_segment_shift;

	bool _sample_variance;

public:
	SegmentWorkerTemplate ( std::list<INPUT_TYPE> & data , std::vector<double> & q_th_order_fluctuation , std::ofstream & output_stream , unsigned int segment_size , bool moving_segment_adjustment = false , unsigned int moving_segment_shift = 2 , bool sample_variance = false )
	: _data ( data ) , _q_th_order_fluctuation ( q_th_order_fluctuation ) , _output_stream ( output_stream ) , _segment_size ( segment_size ) , _moving_segment_adjustment ( moving_segment_adjustment ) , _moving_segment_shift ( moving_segment_shift ) , _sample_variance ( sample_variance )
	{
	}

	virtual ~SegmentWorkerTemplate ()
	{
	}

	std::list<INPUT_TYPE> & GetData ()
	{
		return _data;
	}

	std::vector<double> & GetQthOrderFluctuation ()
	{
		return _q_th_order_fluctuation;
	}

	unsigned int GetSegmentSize ()
	{
		return _segment_size;
	}

	bool GetMovingSegmentAdjustment ()
	{
		return _moving_segment_adjustment;
	}

	unsigned int GetMovingSegmentShift ()
	{
		return _moving_segment_shift;
	}

	bool GetSampleVariance ()
	{
		return _sample_variance;
	}

	std::ofstream & GetOutputStream ()
	{
		return _output_stream;
	}

};


template <class INPUT_TYPE , class INTERNAL_TYPE , class INTERMEDIATE_TYPE , class OUTPUT_TYPE >
class SegmentDFAWorkerTemplate : public SegmentWorkerTemplate<INPUT_TYPE>
#ifdef QTHREADS
, public QThread
#endif
{
public:
	typedef void ( * FLUCTUATION_ANALYSER_TYPE ) ( std::vector<INTERMEDIATE_TYPE> & y_data , unsigned int order_of_method , std::list<INTERMEDIATE_TYPE> & coefficients , boost::numeric::ublas::matrix< INTERNAL_TYPE > * ); 

	typedef INTERMEDIATE_TYPE ( * POLYNOMIAL_TYPE ) ( INTERMEDIATE_TYPE x , std::list<INTERMEDIATE_TYPE> & coefficients , INTERMEDIATE_TYPE domain );

protected:

	std::shared_ptr<std::vector<OUTPUT_TYPE *> > _DFAvariance_result;

	unsigned int _order_of_method;

#ifdef MICROTHREADING
	int _microthreads;
#endif

	FLUCTUATION_ANALYSER_TYPE _fluctuation_analyzer;

	POLYNOMIAL_TYPE _polynomial;

	bool _old_algorithm;

public:

	SegmentDFAWorkerTemplate ( std::list<INPUT_TYPE> & data , std::vector<double> & q_th_order_fluctuation , std::shared_ptr< std::vector<OUTPUT_TYPE *> > DFAvariance_result , std::ofstream & output_stream , unsigned int order_of_method , unsigned int segment_size , bool moving_segment_adjustment = false , unsigned int moving_segment_shift = 2 , bool sample_variance = false , FLUCTUATION_ANALYSER_TYPE fluctuation_analyzer = &statistical_functions::PolynomialRegressionTemplate < std::vector<INTERMEDIATE_TYPE> , INTERMEDIATE_TYPE ,  boost::numeric::ublas::matrix< INTERNAL_TYPE > , boost::numeric::ublas::vector<INTERNAL_TYPE> > , POLYNOMIAL_TYPE polynomial  = &statistical_functions::PolynomialTemplate<INTERMEDIATE_TYPE>
#ifdef MICROTHREADING
, int microthreads
#endif
, bool old_algorithm = false )
	: SegmentWorkerTemplate<INPUT_TYPE> ( data , q_th_order_fluctuation , output_stream ,  segment_size , moving_segment_adjustment ,  moving_segment_shift , sample_variance ) , _DFAvariance_result ( DFAvariance_result ) , _order_of_method ( order_of_method ) , _fluctuation_analyzer ( fluctuation_analyzer )  , _polynomial ( polynomial )
#ifdef MICROTHREADING
, _microthreads ( microthreads )
#endif
, _old_algorithm ( old_algorithm )
	{
	}

	~SegmentDFAWorkerTemplate ()
	{
		_DFAvariance_result.reset();
	}

	unsigned int GetOrderOfMethod ()
	{
		return _order_of_method;
	}

	FLUCTUATION_ANALYSER_TYPE & GetFluctuationAnalyser ()
	{
		return _fluctuation_analyzer;
	}

	POLYNOMIAL_TYPE & GetPolynomial ()
	{
		return _polynomial;
	}


#ifdef QTHREADS
	void run()
	{
		(* this) ( );
	}
#endif

	void operator () ()
	{
// 		std::locale locales ("C");
// 		std::locale::global(locales);
// 		std::cout.imbue(locales);
// 		std::cerr.imbue(locales);

		if ( _old_algorithm == true )
		{

		}
		else
		{
// new algorithm
// register a cleanup function
			boost::this_thread::at_thread_exit ( statistical_functions::detrended_fluctuation_analysis::enqueue_thread_for_cleanup_template );

// new algorithm calculation that does not make list of all segments
			unsigned int moving_segment_shift = SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize ();
#ifdef DEBUG
			unsigned int counter = 0;
#endif

			if ( SegmentWorkerTemplate<INPUT_TYPE>::GetMovingSegmentAdjustment () == true )
			{
				if ( SegmentWorkerTemplate<INPUT_TYPE>::GetMovingSegmentShift () <= SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () )
				{
					moving_segment_shift = SegmentWorkerTemplate<INPUT_TYPE>::GetMovingSegmentShift ();
				}
				else
				{
					moving_segment_shift = SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize ();
				}
			}

	// initial setup of positions of q-th order fluctiation function
			std::vector<INTERMEDIATE_TYPE> fluctuation_function ( SegmentWorkerTemplate<INPUT_TYPE>::GetQthOrderFluctuation (). size ( ) );
			for ( unsigned int i = 0 ; i < SegmentWorkerTemplate<INPUT_TYPE>::GetQthOrderFluctuation () . size ( ) ; i ++ )
			{
				fluctuation_function [ i ] = static_cast <INTERMEDIATE_TYPE> ( 0 );
			}

#ifdef DEBUG
			{
				boost::posix_time::ptime time_of_log ( boost::posix_time::second_clock::local_time() );
				std::cerr << time_of_log << " Thread start with segment size " << SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () << std::endl;
			}
#endif

// preparation of cached matrices for regression procedure
			boost::numeric::ublas::matrix< INTERNAL_TYPE > transformation_y_data_cached;
			PrepareMatricesForRegression<boost::numeric::ublas::matrix< INTERNAL_TYPE > > ( SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () , GetOrderOfMethod () , transformation_y_data_cached );

// preparation of segments starting from first data to last
			unsigned int count_segments = 0;
			typename std::list<INPUT_TYPE>::iterator end_iter = SegmentWorkerTemplate<INPUT_TYPE>::GetData().end();
			for ( typename std::list<INPUT_TYPE>::iterator iter_data = SegmentWorkerTemplate<INPUT_TYPE>::GetData().begin() ; iter_data != end_iter ; )
			{
				std::vector<INTERMEDIATE_TYPE> new_segment ( SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () );
				unsigned int count_elements_of_segment = 0;

// copy data to new destination to process
				for ( typename std::list<INPUT_TYPE>::iterator iter_data_segment = iter_data ;  ( count_elements_of_segment < SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () ) && (  iter_data_segment != end_iter ) ; ++ iter_data_segment , ++ count_elements_of_segment )
				{
					INPUT_TYPE actual_data = ( * iter_data_segment );
					new_segment [ count_elements_of_segment ] = static_cast<INTERMEDIATE_TYPE> ( actual_data );
				}

				if ( count_elements_of_segment /*new_segment.size() */ == SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () )
				{
#ifdef DEBUG
					{
						std::stringstream data_for_debug;
						for ( unsigned int i = 0 ; i < SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () ; ++ i )
						{
							data_for_debug << new_segment [ i ] << " ";
						}
						std::string new_segment_data_for_debug = data_for_debug.str();
						boost::posix_time::ptime time_of_log ( boost::posix_time::second_clock::local_time() );
						math_model::debug::DebugLog ( 6 ) << time_of_log << " segment: " << new_segment_data_for_debug << std::endl;
					}
#endif

					std::list<INTERMEDIATE_TYPE> coefficients ( GetOrderOfMethod () );
					GetFluctuationAnalyser () ( new_segment, GetOrderOfMethod () , coefficients , & transformation_y_data_cached );
#ifdef DEBUG
					{
						std::stringstream data_for_debug2;
						for ( typename std::list<INTERMEDIATE_TYPE>::iterator iter_coefficients = coefficients.begin() ; iter_coefficients != coefficients.end() ; ++ iter_coefficients )
						{
							data_for_debug2 << ( * iter_coefficients ) << " ";
						}
						std::string new_coefficient_data = data_for_debug2.str();
						boost::posix_time::ptime time_of_log2 ( boost::posix_time::second_clock::local_time() );
						math_model::debug::DebugLog ( 6 ) << time_of_log2 << " coefficients: " << new_coefficient_data << std::endl;
					}
#endif

// calculation of regression
					std::vector<INTERMEDIATE_TYPE> differences ( new_segment.size() );
					unsigned int counter_of_position_in_segment = 0;
					typename std::vector<INTERMEDIATE_TYPE>::iterator end_segment_iter = new_segment.end();
					for ( typename std::vector<INTERMEDIATE_TYPE>::iterator iter = new_segment.begin() ; iter != end_segment_iter ; ++ iter )
					{
						INTERMEDIATE_TYPE actual_value = ( * iter );
						INTERMEDIATE_TYPE approximation = GetPolynomial () ( counter_of_position_in_segment , coefficients , SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () );

						differences [ counter_of_position_in_segment ] = actual_value - approximation;
// 						( * iter ) = actual_value - approximation;
#ifdef DEBUG
						{
							boost::posix_time::ptime time_of_log ( boost::posix_time::second_clock::local_time() );
							math_model::debug::DebugLog ( 6 ) << time_of_log << " " << SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () << " " << count_segments << " " << counter_of_position_in_segment << " " << actual_value << " " << approximation << " " << actual_value - approximation << std::endl;
						}
#endif
						if ( SegmentWorkerTemplate<INPUT_TYPE>::GetOutputStream ().is_open ( ) )
						{
							SegmentWorkerTemplate<INPUT_TYPE>::GetOutputStream () << SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () << " " << count_segments << " " << counter_of_position_in_segment << " " << actual_value << " " << approximation << " " << actual_value - approximation << std::endl;
						}
						++ counter_of_position_in_segment;
					}
					++ count_segments;

	#ifdef DEBUG
					{
						boost::posix_time::ptime time_of_log = boost::posix_time::second_clock::local_time();
						std::cerr << time_of_log << " Regression in segments finished with segment size " << SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () << std::endl;
					}
	#endif

// calculation of fluctuation function
					INTERMEDIATE_TYPE segment_variance = 0;
					{
						typename std::vector<INTERMEDIATE_TYPE>::iterator end_segment_iter = differences.end();
						for ( typename std::vector<INTERMEDIATE_TYPE>::iterator iter = differences.begin() ; iter != end_segment_iter ; ++ iter )
						{
							INTERMEDIATE_TYPE actual_value = ( * iter );

							segment_variance += ( actual_value * actual_value );
						}
					}

// normalization to get sampled variance or variance
					if ( SegmentWorkerTemplate<INPUT_TYPE>::GetSampleVariance () == false )
					{
						segment_variance /= SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize ();
					}
					else
					{
						segment_variance /= ( SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () - 1 );
					}

					for ( unsigned int i = 0 ; i < SegmentWorkerTemplate<INPUT_TYPE>::GetQthOrderFluctuation() . size ( ) ; i ++ )
					{
						if ( SegmentWorkerTemplate<INPUT_TYPE>::GetQthOrderFluctuation() [ i ] != 0 )
						{
							fluctuation_function [ i ] += math_functions::multiprecision::power<INTERMEDIATE_TYPE> ( segment_variance , static_cast<INTERMEDIATE_TYPE> ( SegmentWorkerTemplate<INPUT_TYPE>::GetQthOrderFluctuation() [ i ] / 2.0 ) );
// 							fluctuation_function [ i ] += boost::math::powm1 ( segment_variance , SegmentWorkerTemplate<INPUT_TYPE>::GetQthOrderFluctuation() [ i ] / 2.0 ) + 1;
						}
						else
						{
							fluctuation_function [ i ] += math_functions::multiprecision::logarithm<INTERMEDIATE_TYPE> ( math_functions::multiprecision::absolute_value<INTERMEDIATE_TYPE> ( segment_variance ) );
// 							fluctuation_function [ i ] += boost::math::log1p ( fabs ( segment_variance - 1.0 ) );
						}
					}
				}

// shifting iterator of start of new segment
				for ( unsigned int i = 0 ; ( ( i < moving_segment_shift ) && ( iter_data != end_iter ) ) ;  iter_data ++, i ++ );
			}
// averaging data
			for ( unsigned int i = 0 ; i < SegmentWorkerTemplate<INPUT_TYPE>::GetQthOrderFluctuation() . size ( ) ; i ++ )
			{
				if ( SegmentWorkerTemplate<INPUT_TYPE>::GetQthOrderFluctuation() [ i ] != 0.0 )
				{
					INTERMEDIATE_TYPE data = fluctuation_function [ i ] / count_segments;
					* ( ( * _DFAvariance_result ) [ i ] ) = static_cast<OUTPUT_TYPE> ( math_functions::multiprecision::power<INTERMEDIATE_TYPE> ( data , static_cast<INTERMEDIATE_TYPE> ( 1.0 / SegmentWorkerTemplate<INPUT_TYPE>::GetQthOrderFluctuation() [ i ] ) ) );
// 					* ( ( * _DFAvariance_result ) [ i ] ) = boost::math::powm1 ( data , 1.0 / _q_th_order_fluctuation [ i ] ) + 1;
				}
				else
				{
					INTERMEDIATE_TYPE data = fluctuation_function [ i ] / ( 2 * count_segments );
					* ( ( * _DFAvariance_result ) [ i ] ) = static_cast<OUTPUT_TYPE> ( math_functions::multiprecision::exponential<INTERMEDIATE_TYPE> ( data ) );
// 					* ( ( * _DFAvariance_result ) [ i ] ) = boost::math::expm1 ( data ) + 1;
				}

				if ( math_functions::multiprecision::isnan<OUTPUT_TYPE> ( * ( ( * _DFAvariance_result ) [ i ] ) ) )
				{
					std::cerr << "Problem with result. Segment: " << SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () << " q: " << SegmentWorkerTemplate<INPUT_TYPE>::GetQthOrderFluctuation() [ i ] << std::endl;
				}
			}
		}

#ifdef DEBUG
		{
			boost::posix_time::ptime time_of_log ( boost::posix_time::second_clock::local_time() );
			std::cerr << time_of_log << " Thread finished with segment size " << SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () << std::endl;
		}
#endif

		{
			boost::posix_time::ptime actual_time ( boost::posix_time::microsec_clock::local_time() );
			std::cout << actual_time << " Finalising of the segment size " << SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () << " of order " << GetOrderOfMethod () << std::endl;
		}
	}
};


template <class INPUT_TYPE>
void IntegrationOfDatasetTemplate ( std::list<INPUT_TYPE> & data )
{
// determination  of mean
	INPUT_TYPE complete_sum = 0;
	for ( typename std::list<INPUT_TYPE>::iterator iter = data.begin() ; iter != data.end() ; ++ iter )
	{
		INPUT_TYPE & value = ( * iter );
		complete_sum += value;
	}
	INPUT_TYPE mean = complete_sum / data.size();

// integration of dataset
	INPUT_TYPE sum = 0;
	for ( typename std::list<INPUT_TYPE>::iterator iter = data.begin() ; iter != data.end() ; ++ iter )
	{
		INPUT_TYPE & value = ( * iter );
		sum += ( value - mean );
		value = sum;
	}
}

/// Segment polynomial MF-DCCA worker
template< class INPUT_TYPE , class INTERNAL_TYPE , class INTERMEDIATE_TYPE , class OUTPUT_TYPE>
class SegmentDCCAWorkerTemplate : public SegmentWorkerTemplate < INPUT_TYPE >
#ifdef QTHREADS
, public QThread
#endif
{
public:
	typedef void ( * FLUCTUATION_ANALYSER_TYPE ) ( std::vector<INTERMEDIATE_TYPE> & y_data , unsigned int order_of_method , std::list<INTERMEDIATE_TYPE> & coefficients , boost::numeric::ublas::matrix< INTERNAL_TYPE > * ); 

	typedef INTERMEDIATE_TYPE ( * POLYNOMIAL_TYPE ) ( INTERMEDIATE_TYPE x , std::list<INTERMEDIATE_TYPE> & coefficients , INTERMEDIATE_TYPE domain );

protected:
	std::shared_ptr<std::vector< std::array < double * , 4 > > > _DCCA_result;

	std::list<INPUT_TYPE> & _data2;

	unsigned int _order_of_method;

	FLUCTUATION_ANALYSER_TYPE _fluctuation_analyzer;

	POLYNOMIAL_TYPE _polynomial;

// 	void ( * _fluctuation_analyzer ) ( std::vector<double> & y_data , unsigned int order_of_method , std::list<double> & coefficients , double , math_functions::real_matrix * );

// 	double ( * _polynomial ) ( double x , std::list<double> & coefficients , double domain );

public:
	SegmentDCCAWorkerTemplate ( std::list<INPUT_TYPE> & data1 , std::list<INPUT_TYPE> & data2 , std::vector<double> & q_th_order_fluctuation , std::shared_ptr<std::vector< std::array < double * , 4 > > > DCCA_result , std::ofstream & output_stream , unsigned int order_of_method , unsigned int segment_size , bool moving_segment_adjustment = false , unsigned int moving_segment_shift = 2 , bool sample_variance = false , FLUCTUATION_ANALYSER_TYPE fluctuation_analyzer = &statistical_functions::PolynomialRegressionTemplate < std::vector<INTERMEDIATE_TYPE> , INTERMEDIATE_TYPE ,  boost::numeric::ublas::matrix< INTERNAL_TYPE > , boost::numeric::ublas::vector<INTERNAL_TYPE> > , POLYNOMIAL_TYPE polynomial  = &statistical_functions::PolynomialTemplate<INTERMEDIATE_TYPE> )
	: SegmentWorkerTemplate < INPUT_TYPE > ( data1 , q_th_order_fluctuation , output_stream ,  segment_size , moving_segment_adjustment ,  moving_segment_shift , sample_variance ) , _DCCA_result ( DCCA_result ) , _data2 ( data2 ) , _order_of_method ( order_of_method ) , _fluctuation_analyzer ( fluctuation_analyzer ) , _polynomial ( polynomial )
	{
		if ( SegmentWorkerTemplate < INPUT_TYPE >::GetData().size() != GetData2().size() )
		{
			std::cout << "SegmentDCCAWorker::SegmentDCCAWorker: Datasets of interest does not have the same size." << std::endl;
			exit ( EXIT_FAILURE );
		}
	}

	unsigned int GetOrderOfMethod ()
	{
		return _order_of_method;
	}

	FLUCTUATION_ANALYSER_TYPE & GetFluctuationAnalyser ()
	{
		return _fluctuation_analyzer;
	}

	POLYNOMIAL_TYPE & GetPolynomial ()
	{
		return _polynomial;
	}

	std::list<INPUT_TYPE> & GetData2 ()
	{
		return _data2;
	}

#ifdef QTHREADS
	void run()
	{
		(* this) ( );
	}
#endif

	void operator () ()
	{
// register a cleanup function
		boost::this_thread::at_thread_exit ( statistical_functions::detrended_fluctuation_analysis::enqueue_thread_for_cleanup_template );

// new algorithm calculation that does not make list of all segments
		unsigned int moving_segment_shift = SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize ();
#ifdef DEBUG
		unsigned int counter = 0;
#endif

		if ( SegmentWorkerTemplate<INPUT_TYPE>::GetMovingSegmentAdjustment() == true )
		{
			if ( SegmentWorkerTemplate<INPUT_TYPE>::GetMovingSegmentShift () <= SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () )
			{
				moving_segment_shift = SegmentWorkerTemplate<INPUT_TYPE>::GetMovingSegmentShift ();
			}
			else
			{
				moving_segment_shift = SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize ();
			}
		}

// initial setup of positions of q-th order fluctiation function
		std::vector<INTERMEDIATE_TYPE> fluctuation_function ( SegmentWorkerTemplate<INPUT_TYPE>::GetQthOrderFluctuation() . size ( ) );
		std::vector<INTERMEDIATE_TYPE> fluctuation_function2 ( SegmentWorkerTemplate<INPUT_TYPE>::GetQthOrderFluctuation() . size ( ) );
		std::vector<INTERMEDIATE_TYPE> crosscorrelation_function ( SegmentWorkerTemplate<INPUT_TYPE>::GetQthOrderFluctuation() . size ( ) );
		std::vector<INTERMEDIATE_TYPE> abs_crosscorrelation_function ( SegmentWorkerTemplate<INPUT_TYPE>::GetQthOrderFluctuation() . size ( ) );

		for ( unsigned int i = 0 ; i < SegmentWorkerTemplate<INPUT_TYPE>::GetQthOrderFluctuation() . size ( ) ; i ++ )
		{
			fluctuation_function [ i ] = static_cast <INTERMEDIATE_TYPE> ( 0 );
			fluctuation_function2 [ i ] = static_cast <INTERMEDIATE_TYPE> ( 0 );
			crosscorrelation_function [ i ] = static_cast <INTERMEDIATE_TYPE> ( 0 );
			abs_crosscorrelation_function [ i ] = static_cast <INTERMEDIATE_TYPE> ( 0 );
		}

#ifdef DEBUG
		{
			boost::mutex::scoped_lock lock( DFASynchronization::GetDFASynchronization ().GetTextOutputMutex () );
			boost::posix_time::ptime time_of_log ( boost::posix_time::second_clock::local_time() );
			std::cerr << time_of_log << " Thread start with segment size " << SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () << std::endl;
		}
#endif

// preparation of cached matrices for regression procedure
		boost::numeric::ublas::matrix< INTERNAL_TYPE > transformation_y_data_cached;
		PrepareMatricesForRegression<boost::numeric::ublas::matrix< INTERNAL_TYPE > > ( SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize ()  , GetOrderOfMethod () , transformation_y_data_cached );

// preparation of segments starting from first data to last
		unsigned int count_segments = 0;
		typename std::list<INPUT_TYPE>::iterator end_iter = SegmentWorkerTemplate<INPUT_TYPE>::GetData().end();
		typename std::list<INPUT_TYPE>::iterator iter_data2 = GetData2().begin();
		for ( typename std::list<INPUT_TYPE>::iterator iter_data = SegmentWorkerTemplate<INPUT_TYPE>::GetData().begin() ; iter_data != end_iter ; )
		{
			std::vector<INTERMEDIATE_TYPE> new_segment ( SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () );
			std::vector<INTERMEDIATE_TYPE> new_segment2 ( SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () );
			unsigned int count_elements_of_segment = 0;

			typename std::list<INPUT_TYPE>::iterator iter_data_segment2 = iter_data2;
			for ( typename std::list<INPUT_TYPE>::iterator iter_data_segment = iter_data ;  ( count_elements_of_segment < SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () ) && (  iter_data_segment != end_iter ) ; iter_data_segment ++ , count_elements_of_segment ++ , iter_data_segment2 ++ )
			{
				INPUT_TYPE actual_data = ( * iter_data_segment );
				INPUT_TYPE actual_data2 = ( * iter_data_segment2 );
				new_segment [ count_elements_of_segment ] = actual_data;
				new_segment2 [ count_elements_of_segment ] = actual_data2;
			}

			if ( count_elements_of_segment == SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () )
			{
#ifdef DEBUG
				{
					boost::mutex::scoped_lock lock( DFASynchronization::GetDFASynchronization ().GetTextOutputMutex () );
					std::stringstream data_for_debug;
					for ( unsigned int i = 0 ; i < SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () ; i ++ )
					{
						data_for_debug << new_segment [ i ] << " ";
					}
					std::string new_segment_data_for_debug = data_for_debug.str();
					boost::posix_time::ptime time_of_log ( boost::posix_time::second_clock::local_time() );
					math_model::debug::DebugLog ( 6 ) << time_of_log << " segment: " << new_segment_data_for_debug << std::endl;
				}
#endif

				std::list<INTERMEDIATE_TYPE> coefficients ( GetOrderOfMethod () );
				std::list<INTERMEDIATE_TYPE> coefficients2 ( GetOrderOfMethod () );
				GetFluctuationAnalyser () ( new_segment, GetOrderOfMethod () , coefficients , & transformation_y_data_cached );
				GetFluctuationAnalyser () ( new_segment2, GetOrderOfMethod () , coefficients2 , & transformation_y_data_cached );
#ifdef DEBUG
				{
					boost::mutex::scoped_lock lock( DFASynchronization::GetDFASynchronization ().GetTextOutputMutex () );
					std::stringstream data_for_debug2;
					for ( typename std::list<INTERMEDIATE_TYPE>::iterator iter_coefficients = coefficients.begin() ; iter_coefficients != coefficients.end() ; iter_coefficients ++ )
					{
						data_for_debug2 << ( * iter_coefficients ) << " ";
					}
					std::string new_coefficient_data = data_for_debug2.str();
					boost::posix_time::ptime time_of_log2 ( boost::posix_time::second_clock::local_time() );
					math_model::debug::DebugLog ( 6 ) << time_of_log2 << " coefficients: " << new_coefficient_data << std::endl;
				}
#endif

// calculation of regression
				unsigned int counter_of_position_in_segment = 0;
				typename std::vector<INTERMEDIATE_TYPE>::iterator end_segment_iter = new_segment.end();
				typename std::vector<INTERMEDIATE_TYPE>::iterator iter2 = new_segment2.begin();
				for ( typename std::vector<INTERMEDIATE_TYPE>::iterator iter = new_segment.begin() ; iter != end_segment_iter ; ++ iter , ++ iter2 )
				{
					INTERMEDIATE_TYPE actual_value = ( * iter );
					INTERMEDIATE_TYPE actual_value2 = ( * iter2 );

					INTERMEDIATE_TYPE approximation = GetPolynomial () ( counter_of_position_in_segment , coefficients , SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () );
					INTERMEDIATE_TYPE approximation2 = GetPolynomial () ( counter_of_position_in_segment , coefficients2 , SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () );

					( * iter ) = actual_value - approximation;
					( * iter2 ) = actual_value2 - approximation2;
#ifdef DEBUG
					{
						boost::mutex::scoped_lock lock( DFASynchronization::GetDFASynchronization ().GetTextOutputMutex () );
						boost::posix_time::ptime time_of_log ( boost::posix_time::second_clock::local_time() );
						math_model::debug::DebugLog ( 6 ) << time_of_log << " " << SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () << " " << count_segments << " " << counter_of_position_in_segment << " " << actual_value << " " << approximation << " " << actual_value - approximation << std::endl;
					}
#endif
					if ( SegmentWorkerTemplate<INPUT_TYPE>::GetOutputStream ().is_open ( ) )
					{
						SegmentWorkerTemplate<INPUT_TYPE>::GetOutputStream () << SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () << " " << count_segments << " " << counter_of_position_in_segment << " " << actual_value << " " << approximation << " " << actual_value - approximation << std::endl;
					}
					counter_of_position_in_segment ++;
				}
				count_segments ++;

#ifdef DEBUG
				{
					boost::mutex::scoped_lock lock( DFASynchronization::GetDFASynchronization ().GetTextOutputMutex () );
					boost::posix_time::ptime time_of_log = boost::posix_time::second_clock::local_time();
					std::cerr << time_of_log << " Regression in segments finished with segment size " << SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () << std::endl;
				}
#endif

// calculation of fluctuation function
				INTERMEDIATE_TYPE segment_variance = 0;
				INTERMEDIATE_TYPE segment_variance2 = 0;
				INTERMEDIATE_TYPE segment_covariance = 0;
				INTERMEDIATE_TYPE segment_abs_covariance = 0;
				INTERMEDIATE_TYPE segment_correlation_coefficient = 0;

				iter2 = new_segment2.begin();
				for ( auto iter = new_segment.begin() ; iter != end_segment_iter ; ++ iter , ++ iter2 )
				{
					INTERMEDIATE_TYPE actual_value = ( * iter );
					INTERMEDIATE_TYPE actual_value2 = ( * iter2 );

					segment_variance += ( actual_value * actual_value );
					segment_variance2 += ( actual_value2 * actual_value2 );
					segment_covariance += ( actual_value * actual_value2 );
					segment_abs_covariance += ( math_functions::multiprecision::absolute_value<INTERMEDIATE_TYPE> ( actual_value ) * math_functions::multiprecision::absolute_value<INTERMEDIATE_TYPE>  ( actual_value2 ) );
				}

// normalization to get sampled variance or variance
				if ( SegmentWorkerTemplate<INPUT_TYPE>::GetSampleVariance() == false )
				{
					segment_variance /= SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize ();
					segment_variance2 /= SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize ();
					segment_covariance /= SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize ();
					segment_abs_covariance /= SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize ();
				}
				else
				{
					segment_variance /= ( SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () - 1 );
					segment_variance2 /= ( SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () - 1 );
					segment_covariance /= ( SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () - 1 );
					segment_abs_covariance /= ( SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () - 1 );
				}

				for ( unsigned int i = 0 ; i < SegmentWorkerTemplate<INPUT_TYPE>::GetQthOrderFluctuation() . size ( ) ; i ++ )
				{
					if ( SegmentWorkerTemplate<INPUT_TYPE>::GetQthOrderFluctuation() [ i ] != 0 )
					{
						INTERMEDIATE_TYPE exponent_of_power = SegmentWorkerTemplate<INPUT_TYPE>::GetQthOrderFluctuation() [ i ] / 2;
						fluctuation_function [ i ] += math_functions::multiprecision::power<INTERMEDIATE_TYPE> ( segment_variance , exponent_of_power );
						fluctuation_function2 [ i ] += math_functions::multiprecision::power<INTERMEDIATE_TYPE> ( segment_variance2 , exponent_of_power );
						crosscorrelation_function [ i ] += math_functions::multiprecision::power<INTERMEDIATE_TYPE> ( math_functions::multiprecision::absolute_value<INTERMEDIATE_TYPE> ( segment_covariance ) , exponent_of_power );
						abs_crosscorrelation_function [ i ] += math_functions::multiprecision::power<INTERMEDIATE_TYPE> ( ( segment_abs_covariance ) , exponent_of_power );
					}
					else
					{
						fluctuation_function [ i ] += math_functions::multiprecision::logarithm<INTERMEDIATE_TYPE> ( segment_variance );
						fluctuation_function2 [ i ] += math_functions::multiprecision::logarithm<INTERMEDIATE_TYPE> ( segment_variance2 );
						crosscorrelation_function [ i ] += math_functions::multiprecision::logarithm<INTERMEDIATE_TYPE> ( math_functions::multiprecision::absolute_value<INTERMEDIATE_TYPE> ( segment_covariance ) );
						abs_crosscorrelation_function [ i ] += math_functions::multiprecision::logarithm<INTERMEDIATE_TYPE> ( ( segment_abs_covariance ) );
					}
				}
			}

// shifting iterator of start of new segment
			for ( unsigned int i = 0 ; ( ( i < moving_segment_shift ) && ( iter_data != end_iter ) ) ;  iter_data ++ , iter_data2 ++ , i ++ );
		}

// averaging data
		for ( unsigned int i = 0 ; i < SegmentWorkerTemplate<INPUT_TYPE>::GetQthOrderFluctuation() . size ( ) ; i ++ )
		{
			if ( ( ( * _DCCA_result ) [ i ] [ 0 ] == nullptr ) || ( ( * _DCCA_result ) [ i ] [ 1 ] == nullptr ) || ( ( * _DCCA_result ) [ i ] [ 2 ] == nullptr ) || ( ( * _DCCA_result ) [ i ] [ 3 ] == nullptr ) )
			{
				std::cerr << "SegmentDCCAWorker::operator: Error: Nullptr encoutered" << std::endl;
				std::cerr << "Parameters are: Order: " << _order_of_method << " , Segment size: " << SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () << " , i: " << i << " , affiliated q: " << SegmentWorkerTemplate<INPUT_TYPE>::GetQthOrderFluctuation() [ i ] << " , Pointers are: " << ( * _DCCA_result ) [ i ] [ 0 ] << " , " << ( * _DCCA_result ) [ i ] [ 1 ] << " , " << ( * _DCCA_result ) [ i ] [ 2 ] << " , " << ( * _DCCA_result ) [ i ] [ 3 ] << std::endl;
				exit ( EXIT_FAILURE );
			}

			if ( SegmentWorkerTemplate<INPUT_TYPE>::GetQthOrderFluctuation() [ i ] != 0 )
			{
				INTERMEDIATE_TYPE data1 = fluctuation_function [ i ] / count_segments;
				INTERMEDIATE_TYPE data2 = fluctuation_function2 [ i ] / count_segments;
				INTERMEDIATE_TYPE data3 = crosscorrelation_function [ i ] / count_segments;
				INTERMEDIATE_TYPE data4 = abs_crosscorrelation_function [ i ] / count_segments;

				INTERMEDIATE_TYPE exponent_of_power = 1 / SegmentWorkerTemplate<INPUT_TYPE>::GetQthOrderFluctuation() [ i ];
				( * ( * _DCCA_result ) [ i ] [ 0 ] ) = static_cast<OUTPUT_TYPE> ( math_functions::multiprecision::power<INTERMEDIATE_TYPE> ( data1 , exponent_of_power ) );
				( * ( * _DCCA_result ) [ i ] [ 1 ] ) = static_cast<OUTPUT_TYPE> ( math_functions::multiprecision::power<INTERMEDIATE_TYPE> ( data2 , exponent_of_power ) );
				( * ( * _DCCA_result ) [ i ] [ 2 ] ) = static_cast<OUTPUT_TYPE> ( math_functions::multiprecision::power<INTERMEDIATE_TYPE> ( math_functions::multiprecision::absolute_value<INTERMEDIATE_TYPE> ( data3 ) , exponent_of_power ) );
				( * ( * _DCCA_result ) [ i ] [ 3 ] ) = static_cast<OUTPUT_TYPE> ( math_functions::multiprecision::power<INTERMEDIATE_TYPE> ( math_functions::multiprecision::absolute_value<INTERMEDIATE_TYPE> ( data4 ) , exponent_of_power ) );
			}
			else
			{
				INTERMEDIATE_TYPE data1 = fluctuation_function [ i ] / ( 2 * count_segments );
				INTERMEDIATE_TYPE data2 = fluctuation_function2 [ i ] / ( 2 * count_segments );
				INTERMEDIATE_TYPE data3 = crosscorrelation_function [ i ] / ( 2 * count_segments );
				INTERMEDIATE_TYPE data4 = abs_crosscorrelation_function [ i ] / ( 2 * count_segments );

				( * ( * _DCCA_result ) [ i ] [ 0 ] ) = static_cast<OUTPUT_TYPE> ( math_functions::multiprecision::exponential<INTERMEDIATE_TYPE> ( data1 ) );
				( * ( * _DCCA_result ) [ i ] [ 1 ] ) = static_cast<OUTPUT_TYPE> ( math_functions::multiprecision::exponential<INTERMEDIATE_TYPE> ( data2 ) );
				( * ( * _DCCA_result ) [ i ] [ 2 ] ) = static_cast<OUTPUT_TYPE> ( math_functions::multiprecision::exponential<INTERMEDIATE_TYPE> ( data3 ) );
				( * ( * _DCCA_result ) [ i ] [ 3 ] ) = static_cast<OUTPUT_TYPE> ( math_functions::multiprecision::exponential<INTERMEDIATE_TYPE> ( data4 ) );
			}

			if ( ( math_functions::multiprecision::isnan<OUTPUT_TYPE> ( * ( * _DCCA_result ) [ i ] [ 0 ] ) ) || ( math_functions::multiprecision::isnan<OUTPUT_TYPE> ( * ( * _DCCA_result ) [ i ] [ 1 ] ) ) || ( math_functions::multiprecision::isnan<OUTPUT_TYPE> ( * ( * _DCCA_result ) [ i ] [ 2 ] ) ) || ( math_functions::multiprecision::isnan<OUTPUT_TYPE> ( * ( * _DCCA_result ) [ i ] [ 3 ] ) ) )
			{
				std::cerr << "Problem with result. Segment: " << SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () << " q: " << SegmentWorkerTemplate<INPUT_TYPE>::GetQthOrderFluctuation() [ i ] << std::endl;
				std::cerr << ( * ( * _DCCA_result ) [ i ] [ 0 ] ) << " " << ( * ( * _DCCA_result ) [ i ] [ 1 ] ) << " " << ( * ( * _DCCA_result ) [ i ] [ 2 ] ) << " " << ( * ( * _DCCA_result ) [ i ] [ 3 ] ) << std::endl;
			}
		}

#ifdef DEBUG
		{
			boost::mutex::scoped_lock lock( DFASynchronization::GetDFASynchronization ().GetTextOutputMutex () );
			boost::posix_time::ptime time_of_log = boost::posix_time::second_clock::local_time();
			std::cerr << time_of_log << " Thread finished with segment size " << SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize() << std::endl;
		}
#endif

		{
			boost::posix_time::ptime actual_time ( boost::posix_time::microsec_clock::local_time() );
			std::cout << actual_time << " Finalising of the segment size " << SegmentWorkerTemplate<INPUT_TYPE>::GetSegmentSize () << " of order " << GetOrderOfMethod () << std::endl;
		}
	}
};


template< class INPUT_TYPE , class INTERNAL_TYPE , class INTERMEDIATE_TYPE , class OUTPUT_TYPE>
void DetrendedFluctuationAnalysisTemplate ( std::list<INPUT_TYPE> & data , std::vector<unsigned int> & segment_sizes , std::map<std::tuple<double, unsigned int, unsigned int> , OUTPUT_TYPE > & DFA , std::vector<double> & q_th_order_fluctuation , const std::vector<long> & order_of_method , std::string output_analysis , bool moving_segment_adjustment , unsigned int moving_segment_shift , bool sample_variance , bool double_integration
#ifdef THREADING
, unsigned int number_threads
#ifdef MICROTHREADING
, unsigned int microthreads
#endif
#endif
, typename SegmentDFAWorkerTemplate <INPUT_TYPE, INTERNAL_TYPE, INTERMEDIATE_TYPE, OUTPUT_TYPE>::FLUCTUATION_ANALYSER_TYPE fluctuation_analyzer = &statistical_functions::PolynomialRegressionTemplate < std::vector<INTERMEDIATE_TYPE> , INTERMEDIATE_TYPE ,  boost::numeric::ublas::matrix< INTERNAL_TYPE > , boost::numeric::ublas::vector<INTERNAL_TYPE> > , typename SegmentDFAWorkerTemplate <INPUT_TYPE, INTERNAL_TYPE, INTERMEDIATE_TYPE, OUTPUT_TYPE>::POLYNOMIAL_TYPE  polynomial = &statistical_functions::PolynomialTemplate<INTERMEDIATE_TYPE> )
{
#ifdef DEBUG
	math_model::debug::OpenDebugStream ();
#endif

	std::locale locales ("C");
	std::locale::global(locales);
	std::cout.imbue(locales);
	std::cerr.imbue(locales);

// open output file
	std::ofstream output_analysis_file;
	if ( ! output_analysis.empty() )
	{
		output_analysis_file.open ( output_analysis.c_str() );
	}

	std::list<INPUT_TYPE> integrated_dataset ( data );
	IntegrationOfDatasetTemplate<INPUT_TYPE> ( integrated_dataset );

	if ( double_integration == true )
	{
		IntegrationOfDatasetTemplate<INPUT_TYPE> ( integrated_dataset );
	}

	std::vector<unsigned int>::iterator actual_segment_iter = segment_sizes.begin ();
	std::vector<long>::const_iterator order_of_method_iter = order_of_method.cbegin ();

#ifdef THREADING
#ifdef QTHREADS
	std::list<SegmentDFAWorker *> DFA_threads;
#else
	boost::thread_group * DFA_threads = new boost::thread_group;
	std::list<SegmentDFAWorkerTemplate<INPUT_TYPE , INTERNAL_TYPE , INTERMEDIATE_TYPE , OUTPUT_TYPE> *> workers_for_cleanup;
#endif

	{
		boost::posix_time::ptime actual_time ( boost::posix_time::microsec_clock::local_time() );
		std::cout << actual_time << " Start of the method " << std::endl;
	}

// loop over segment sizes that are going to be iterated
	while ( ( order_of_method_iter != order_of_method.cend () ) )
	{
		{
			boost::posix_time::ptime actual_time ( boost::posix_time::microsec_clock::local_time() );
			std::cout << actual_time << " Start to submit job of the segment size " << ( * actual_segment_iter ) << " of order " << ( * order_of_method_iter ) << std::endl;
		}

		const unsigned int number_of_q_elements = q_th_order_fluctuation.size();
		std::shared_ptr< std::vector< OUTPUT_TYPE * > > DFA_variance_results ( new std::vector< OUTPUT_TYPE * > ( number_of_q_elements ) );

		for ( unsigned int i = 0 ; i < number_of_q_elements ; i ++ )
		{
			DFA [ std::make_tuple ( q_th_order_fluctuation [ i ] , ( * order_of_method_iter ) , ( * actual_segment_iter ) ) ] = static_cast<typename std::map<std::tuple<double, unsigned int, unsigned int> , OUTPUT_TYPE >::mapped_type> ( 0 );

// allocation of data for results
			OUTPUT_TYPE & results = DFA [ std::make_tuple ( q_th_order_fluctuation [ i ] , ( * order_of_method_iter ) , ( * actual_segment_iter ) ) ];

			( * DFA_variance_results ) [ i ] = & ( results );

			if ( ( * DFA_variance_results ) [ i ] == nullptr )
			{
				std::cerr << "Error: We have got nullptr." << std::endl;
				exit (1);
			}
		}

		// sanity check to prevent problems with too large (not enought statistics) and too small segments (detrending failure)
		if ( ( data.size() <= ( * actual_segment_iter ) * ( 4 ) ) && ( ( * order_of_method_iter ) + 1 > ( * actual_segment_iter ) ) )
		{
			std::cerr << "DetrendedFluctuationAnalysis: End index of DFA to large. Higher values are irrelevant due to bad statistics" << std::endl;
			continue;
		}

		SegmentDFAWorkerTemplate<INPUT_TYPE , INTERNAL_TYPE , INTERMEDIATE_TYPE , OUTPUT_TYPE> * thread_worker = new SegmentDFAWorkerTemplate<INPUT_TYPE , INTERNAL_TYPE , INTERMEDIATE_TYPE , OUTPUT_TYPE> ( integrated_dataset , q_th_order_fluctuation , DFA_variance_results , output_analysis_file , ( * order_of_method_iter ) , ( * actual_segment_iter ) , moving_segment_adjustment , moving_segment_shift , sample_variance
#ifdef MICROTHREADING
 , microthreads
#endif
, fluctuation_analyzer , polynomial
		);
		workers_for_cleanup . push_back ( thread_worker );

		{
			size_t number_of_threads;
			size_t number_finished_jobs;

			{
				boost::mutex::scoped_lock lock ( DFASynchronization::GetDFASynchronization().GetAddJobsMutex () );
				boost::thread * new_thread = new boost::thread ( * thread_worker );

				DFA_threads -> add_thread ( new_thread );

				number_of_threads = DFA_threads->size ( );
				number_finished_jobs = DFASynchronization::GetDFASynchronization().GetFinishedJobs().size ( );

				if ( ( number_of_threads - number_finished_jobs ) >= number_threads )
				{
					// wait for a thread to wake up
					DFASynchronization::GetDFASynchronization().GetAddJobsCV().wait ( lock );
				}
			}

// increment job position
			++ actual_segment_iter;

			if ( ( actual_segment_iter == segment_sizes.end () ) && ( order_of_method_iter != order_of_method.end () ) )
			{
				boost::posix_time::ptime actual_time ( boost::posix_time::microsec_clock::local_time() );
				std::cout << actual_time << " We finished to enqueue jobs of the order " << ( * order_of_method_iter ) << std::endl;
				actual_segment_iter = segment_sizes.begin ();
				++ order_of_method_iter;
			}
		}
	}

	DFA_threads -> join_all();
	delete DFA_threads;
	DFASynchronization::GetDFASynchronization().GetFinishedJobs().clear ();
	for ( auto item = workers_for_cleanup.begin() ; item != workers_for_cleanup.end() ; ++ item )
	{
		delete ( * item );
	}

#endif
	if ( output_analysis_file.is_open() )
	{
		output_analysis_file.close ();
	}
}


template< class INPUT_TYPE , class INTERNAL_TYPE , class INTERMEDIATE_TYPE , class OUTPUT_TYPE>
void DetrendedCrossCorrelationAnalysisTemplate ( std::list<INPUT_TYPE> & data1 , std::list<INPUT_TYPE> & data2 , std::vector<unsigned int> & segment_sizes , std::map<std::tuple<double, unsigned int, unsigned int> , std::array<OUTPUT_TYPE,4> > & DCCA , std::vector<double> & q_th_order_fluctuation , std::vector<unsigned int> & order_of_method , std::string output_analysis , bool moving_segment_adjustment , unsigned int moving_segment_shift , bool sample_variance , bool double_integration
#ifdef THREADING
, unsigned int number_threads
#endif
, typename SegmentDFAWorkerTemplate <INPUT_TYPE, INTERNAL_TYPE, INTERMEDIATE_TYPE, OUTPUT_TYPE>::FLUCTUATION_ANALYSER_TYPE fluctuation_analyzer = &statistical_functions::PolynomialRegressionTemplate < std::vector<INTERMEDIATE_TYPE> , INTERMEDIATE_TYPE ,  boost::numeric::ublas::matrix< INTERNAL_TYPE > , boost::numeric::ublas::vector<INTERNAL_TYPE> > , typename SegmentDFAWorkerTemplate <INPUT_TYPE, INTERNAL_TYPE, INTERMEDIATE_TYPE, OUTPUT_TYPE>::POLYNOMIAL_TYPE  polynomial = &statistical_functions::PolynomialTemplate<INTERMEDIATE_TYPE> )
{
#ifdef DEBUG
	math_model::debug::OpenDebugStream ();
#endif

	std::locale::global(std::locale(""));
	std::cout.imbue(std::locale());
	std::cerr.imbue(std::locale());

// iterators of setup values of DCCA
	if ( ( segment_sizes . empty () ) || ( q_th_order_fluctuation . empty () ) || ( order_of_method .empty () ) )
	{
		std::cout << "No jobs to calculate: segment empty." << std::endl;
		exit ( EXIT_FAILURE );
	}

// open output file
	std::ofstream output_analysis_file;
	if ( ! output_analysis.empty() )
	{
		output_analysis_file.open ( output_analysis.c_str() );
	}

	// integration of the datasets
	// dataset 1
	std::list<INPUT_TYPE> integrated_dataset_1 ( data1 );
	IntegrationOfDatasetTemplate<INPUT_TYPE> ( integrated_dataset_1 );

	if ( double_integration == true )
	{
		IntegrationOfDatasetTemplate<INPUT_TYPE> ( integrated_dataset_1 );
	}

	// dataset 2
	std::list<INPUT_TYPE> integrated_dataset_2 ( data2 );
	IntegrationOfDatasetTemplate<INPUT_TYPE> ( integrated_dataset_2 );

	if ( double_integration == true )
	{
		IntegrationOfDatasetTemplate<INPUT_TYPE> ( integrated_dataset_2 );
	}

	std::vector<unsigned int>::iterator actual_segment_iter = segment_sizes.begin ();
	std::vector<unsigned int>::iterator order_of_method_iter = order_of_method.begin ();

#ifdef THREADING
#ifdef QTHREADS
	std::list<SegmentDCCAWorker *> DCCA_threads;
#else
	boost::thread_group * DCCA_threads = new boost::thread_group;
#endif

// loop over segment sizes that are going to be iterated
	while ( ( order_of_method_iter != order_of_method.end () ) )
	{
		{
			boost::posix_time::ptime actual_time ( boost::posix_time::microsec_clock::local_time() );
			std::cout << actual_time << " Start to submit job of the segment size " << ( * actual_segment_iter ) << " of order " << ( * order_of_method_iter ) << std::endl;
		}

		const unsigned int number_of_q_elements = q_th_order_fluctuation.size();
		std::shared_ptr<std::vector< std::array < OUTPUT_TYPE * , 4 > > > DCCA_variance_results ( new std::vector< std::array < OUTPUT_TYPE * , 4 > > ( number_of_q_elements ) );

		for ( unsigned int i = 0 ; i < number_of_q_elements ; i ++ )
		{
// loop over all fluctuations asked, F^{a}_q (s) , F^{b}_q (s) , XF^{a,b}_q (s) , CF^{a,b}_q (s)

			std::array < double , 4 > zero;
			zero.fill ( static_cast<typename std::map<std::tuple<double, unsigned int, unsigned int> , OUTPUT_TYPE >::mapped_type> ( 0 ) );
			DCCA [ std::make_tuple ( q_th_order_fluctuation [ i ] , ( * order_of_method_iter ) , ( * actual_segment_iter ) ) ] = zero;

// allocation of data for results
			std::array < OUTPUT_TYPE , 4 > & results = DCCA [ std::make_tuple ( q_th_order_fluctuation [ i ] , ( * order_of_method_iter ) , ( * actual_segment_iter ) ) ];

			for ( unsigned int j = 0 ; j < 4 ; j ++ )
			{
				( * DCCA_variance_results ) [ i ] [ j ] = & ( results [ j ] );

				if ( ( * DCCA_variance_results ) [ i ] [ j ] == nullptr )
				{
					std::cerr << "Error: We have got nullptr." << std::endl;
					exit (1);
				}
			}
		}

		if ( ( data1.size() <= ( * actual_segment_iter ) * ( 4 ) ) && ( ( * order_of_method_iter ) + 1 > ( * actual_segment_iter ) ) )
		{
			std::cerr << "DetrendedCrossCorrelationAnalysis: End index of DFA to large. Higher values are irrelevant due to bad statistics" << std::endl;
			continue;
		}

		SegmentDCCAWorkerTemplate<INPUT_TYPE , INTERNAL_TYPE , INTERMEDIATE_TYPE , OUTPUT_TYPE> * thread_worker = new SegmentDCCAWorkerTemplate<INPUT_TYPE , INTERNAL_TYPE , INTERMEDIATE_TYPE , OUTPUT_TYPE> ( integrated_dataset_1 , integrated_dataset_2 , q_th_order_fluctuation , DCCA_variance_results , output_analysis_file , ( * order_of_method_iter ) , ( * actual_segment_iter ) , moving_segment_adjustment , moving_segment_shift , sample_variance
#ifdef MICROTHREADING
 , microthreads
#endif
, fluctuation_analyzer , polynomial
		);

		{
			size_t number_of_threads;
			size_t number_finished_jobs;

			{
				boost::mutex::scoped_lock lock ( DFASynchronization::GetDFASynchronization().GetAddJobsMutex () );
				boost::thread * new_thread = new boost::thread ( * thread_worker );

				DCCA_threads -> add_thread ( new_thread );

				number_of_threads = DCCA_threads->size ( );
				number_finished_jobs = DFASynchronization::GetDFASynchronization().GetFinishedJobs().size ( );

				if ( ( number_of_threads - number_finished_jobs ) >= number_threads )
				{
					// wait for a thread to wake up
					DFASynchronization::GetDFASynchronization().GetAddJobsCV().wait ( lock );
				}
			}

// increment job position
			++ actual_segment_iter;

			if ( ( actual_segment_iter == segment_sizes.end () ) && ( order_of_method_iter != order_of_method.end () ) )
			{
				boost::posix_time::ptime actual_time ( boost::posix_time::microsec_clock::local_time() );
				std::cout << actual_time << " We finished to enqueue jobs of the order " << ( * order_of_method_iter ) << std::endl;
				actual_segment_iter = segment_sizes.begin ();
				++ order_of_method_iter;
			}
		}
	}

	DCCA_threads -> join_all();
	delete DCCA_threads;
	DFASynchronization::GetDFASynchronization().GetFinishedJobs().clear ();

#endif
	if ( output_analysis_file.is_open() )
	{
		output_analysis_file.close ();
	}
}

}

}

#endif
