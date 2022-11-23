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

#include <iostream>
#include <algorithm> 
#include <array>
#include <assert.h>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <mutex>
#include <thread>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#ifdef QTHREADS
#include <qt4/QtCore/QThread>
#endif

#include <boost/math/special_functions/sign.hpp>

#ifdef THREADING
#include <boost/thread.hpp>
#include <boost/date_time.hpp>
#include <boost/thread/shared_mutex.hpp>
#include <boost/thread/condition_variable.hpp>
#endif

#include <boost/math/special_functions.hpp>

#include "fluctuation_analyzer.h"
#include "regression.h"
#include "debug_support.h"

boost::mutex work_list_iterator_mutex;
std::list< std::list<double> * >::iterator work_list_iterator;
std::list< std::list<double> * >::iterator end_work_iterator;
boost::mutex text_output_mutex;
static bool finish_threads = false;
boost::mutex finish_thread_mutex;
boost::mutex map_write_mutex;

// variables for conditional variables
boost::mutex add_jobs_mutex;
boost::condition_variable add_jobs_CV;
std::list< boost::thread::id > finished_jobs;


namespace statistical_functions
{

/// Regression worker for microthreading
class PolynomialRegressionWorker
#ifdef QTHREAD
:public QThread
#endif
{
protected:
	int _id;
	std::list< std::list<double> * > & _list;
	std::ofstream & _output_stream;
	int _order_of_method;
	int _segment_size;

public:
	PolynomialRegressionWorker ( int id , std::list< std::list<double> * > & list , std::ofstream & output_stream , int order_of_method , int segment_size )
	: _id ( id ) , _list ( list ) , _output_stream ( output_stream ) , _order_of_method ( order_of_method ) , _segment_size ( segment_size )
	{
	}

#ifdef QTHREADS
	void run()
	{
		(* this) ( );
	}
#endif

	void operator () ()
	{
#ifdef DEBUG
		std::cout << "Thread id: " << _id << " just started" << std::endl;
#endif
		do
		{
#ifdef DEBUG
		std::cout << "Thread id: " << _id << " the main loop started" << std::endl;
#endif
			work_list_iterator_mutex.lock();
			if ( work_list_iterator != end_work_iterator )
			{
				std::list<double> * actual_segment = ( * work_list_iterator );
				work_list_iterator ++;
				work_list_iterator_mutex.unlock();
#ifdef DEBUG
				std::cout << "Thread id: " << _id << " work list iterator pass job here" << std::endl;
#endif

				std::list<double> coefficients;
				PolynomialRegression<std::list<double>,std::list<double>::iterator> ( * actual_segment, _order_of_method , coefficients );
#ifdef DEBUG
				std::cout << "Thread id: " << _id << " regression finished" << std::endl;
#endif

				unsigned int count_i = 0;
				for ( std::list<double>::iterator iter = actual_segment ->begin() ; iter != actual_segment->end() ; ++ iter )
				{
					double actual_value = ( * iter );
					double approximation = Polynomial ( count_i , coefficients );
					( * iter ) = actual_value - approximation;

					if ( _output_stream.is_open ( ) )
					{
						text_output_mutex.lock();
						_output_stream << _segment_size << " " << count_i << " " << actual_value << " " << approximation << " " << actual_value - approximation << std::endl;
						text_output_mutex.unlock();
					}
					count_i ++;
				}
#ifdef DEBUG
				std::cout << "Thread id: " << _id << " calculation differences" << std::endl;
#endif
			}
			else
			{
#ifdef DEBUG
				std::cout << "Thread id: " << _id << " test to finish the thread" << std::endl;
#endif
				work_list_iterator_mutex.unlock();
				finish_thread_mutex.lock();
				if ( finish_threads == true )
				{
					finish_thread_mutex.unlock();
#ifdef DEBUG
				std::cout << "Thread id: " << _id << " going to finish" << std::endl;
#endif
					return;
				}
				finish_thread_mutex.unlock();

#ifdef DEBUG
				std::cout << "Thread id: " << _id << " going to sleep" << std::endl;
#endif
				boost::posix_time::seconds workTime ( 5 );
				boost::this_thread::sleep(workTime);
			}
		}
		while ( true );
	}
};


/// Segment sifting MF-EMD-DFA 
class SegmentEMD_DFAWorker : public SegmentWorker
#ifdef QTHREADS
, public QThread
#endif
{
protected:
	std::vector<double *> _DFAvariance_result;
	bool ( * _sifting_method ) ( std::list<double> & , std::list<double>::iterator , std::list<double>::iterator , std::list<double> & , std::list<double> & , std::list<double>::iterator &, std::list<double>::iterator &, bool , bool , bool , double );

	int _minimal_number_of_sifting;

	int _maximal_number_of_sifting;
// TODO: finish the class declaration as well as definition

public:
	SegmentEMD_DFAWorker ( std::list<double> & data , std::vector<double> & q_th_order_fluctuation , std::vector<double *> DFAvariance_result , std::ofstream & output_stream , unsigned int segment_size , bool moving_segment_adjustment = false , unsigned int moving_segment_shift = 2 , bool sample_variance = false , int minimal_sifting = 0 , bool ( * sifting_method ) ( std::list<double> & , std::list<double>::iterator , std::list<double>::iterator , std::list<double> & , std::list<double> & , std::list<double>::iterator &, std::list<double>::iterator &, bool , bool , bool , double ) = &SiftingOfDataset<std::list<double>, std::list<double>::iterator>
#ifdef MICROTHREADING
, int microthreads
#endif
, int maximal_number_of_sifting = 1
)
	: SegmentWorker ( data , q_th_order_fluctuation , output_stream , segment_size , moving_segment_adjustment ,  moving_segment_shift , sample_variance ) , _DFAvariance_result ( DFAvariance_result ) , _sifting_method ( sifting_method ) , _minimal_number_of_sifting ( minimal_sifting ) , _maximal_number_of_sifting ( maximal_number_of_sifting )
// 	, _maximal_number_of_sifting ( floor( log2 ( segment_size ) ) + 1 )
#ifdef MICROTHREADING
, _microthreads ( microthreads )
#endif
	{
	}

	void operator () ()
	{
// new algorithm calculation that does not make list of all segments
		unsigned int moving_segment_shift = _segment_size;

		if ( _moving_segment_adjustment == true )
		{
			if ( _moving_segment_shift <= _segment_size )
			{
				moving_segment_shift = _moving_segment_shift;
			}
			else
			{
				moving_segment_shift = _segment_size;
			}
		}

	// initial setup of positions of q-th order fluctiation function
		std::vector<double> fluctuation_function ( _q_th_order_fluctuation . size ( ) );
		for ( unsigned int i = 0 ; i < _q_th_order_fluctuation . size ( ) ; i ++ )
		{
			fluctuation_function [ i ] = 0;
		}

#ifdef DEBUG
		{
			boost::posix_time::ptime time_of_log ( boost::posix_time::second_clock::local_time() );
			std::cerr << time_of_log << " Thread start with segment size " << _segment_size << std::endl;
		}
#endif

// preparation of segments starting from first data to the last
		unsigned int count_segments = 0;

// walking through dataset
		std::list<double>::iterator end_iter = _data.end();
		for ( std::list<double>::iterator iter_data = _data.begin() ; iter_data != end_iter ; )
		{
			std::list<double> new_segment;
			std::list<double> sift;
			unsigned int count_elements_of_segment = 0;

// filling data of investigated segment
			std::list<double>::iterator end_iter2 = iter_data; 
			for ( auto iter_data_segment = iter_data ; ( count_elements_of_segment < _segment_size ) && (  iter_data_segment != end_iter ) ; iter_data_segment ++ , count_elements_of_segment ++ )
			{
				double actual_data = ( * iter_data_segment );
				new_segment.push_back ( actual_data );
				++ end_iter2;
			}

			if ( count_elements_of_segment == _segment_size )
			{
#ifdef DEBUG
				std::stringstream data_for_debug;
				for ( std::list<double>::iterator iter = new_segment.begin() ; iter != new_segment.end() ; ++ iter )
				{
					data_for_debug << ( * iter ) << " ";
				}
				std::string new_segment_data_for_debug = data_for_debug.str();
				{
					boost::posix_time::ptime time_of_log ( boost::posix_time::second_clock::local_time() );
					math_model::debug::DebugLog ( 6 ) << time_of_log << " segment: " << new_segment_data_for_debug << std::endl;
				}
#endif
				bool decision;
				int round = 1;
				std::list<double>::iterator new_segment_begin;
				std::list<double>::iterator new_segment_end;

// initial sifting of dataset
				decision = _sifting_method ( _data , iter_data , end_iter2 , sift , new_segment , new_segment_begin , new_segment_end , true , true , false , 1e-6 );
#ifdef DEBUG
				int counter = 0;
				std::cerr << "Starting at value " << ( * new_segment_begin ) << " ending at " << ( * new_segment_end ) << std::endl;
				std::list<double>::iterator sift_iter = sift.begin ();
				for ( std::list<double>::iterator iter = new_segment.begin() ; iter != new_segment.end() ; ++ iter )
				{
					std::cerr << round << " " << counter << " " << ( * iter ) << " " << ( * sift_iter ) << std::endl;
					++ counter;
					++ sift_iter;
				}
#endif

				while ( ( ( _minimal_number_of_sifting <= round ) || ( decision ) ) && ( round < _maximal_number_of_sifting ) )
				{
					round ++;
					decision = _sifting_method ( new_segment , new_segment.begin() , new_segment.end() , sift , new_segment , new_segment_begin , new_segment_end , true , false , false , 1e-6 );

#ifdef DEBUG
					int counter = 0;
					std::cerr << "Starting at value " << ( * new_segment_begin ) << " ending at " << ( * new_segment_end ) << std::endl;
					std::list<double>::iterator sift_iter = sift.begin ();
					for ( std::list<double>::iterator iter = new_segment.begin() ; iter != new_segment.end() ; ++ iter )
					{
					std::cerr << round << " " << counter << " " << ( * iter ) << " " << ( * sift_iter ) << std::endl;
					++ counter;
					++ sift_iter;
					}
#endif
				}
				count_segments ++;

#ifdef DEBUG
				{
					boost::posix_time::ptime time_of_log = boost::posix_time::second_clock::local_time();
					std::cerr << time_of_log << " Empirical mode decomposition of segments finished with segment size " << _segment_size << std::endl;
				}
#endif

// calculation of fluctuation function
				double segment_variance = 0;
				auto end_segment_iter = new_segment.end();
				for ( auto iter = new_segment_begin ; iter != end_segment_iter ; ++ iter )
				{
					double actual_value = ( * iter );

					segment_variance += ( actual_value * actual_value );
				}
				if ( _sample_variance == false )
				{
					segment_variance /= _segment_size;
				}
				else
				{
					segment_variance /= ( _segment_size - 1 );
				}

				for ( unsigned int i = 0 ; i < _q_th_order_fluctuation . size ( ) ; i ++ )
				{
					if ( _q_th_order_fluctuation [ i ] != 0.0 )
					{
						fluctuation_function [ i ] += pow ( segment_variance , _q_th_order_fluctuation [ i ] / 2.0 );
					}
					else
					{
						fluctuation_function [ i ] += log ( segment_variance );
					}
				}
#ifdef DEBUG
				std::cout << "Segment size, rounds, segment variance: " << _segment_size << " " << round << " " << segment_variance << std::endl;
#endif

#ifdef DEBUG
				{
					boost::posix_time::ptime time_of_log ( boost::posix_time::second_clock::local_time() );
					math_model::debug::DebugLog ( 6 ) << time_of_log << " segment size: " << _segment_size << " round: " << round << " segment variance: " << segment_variance << std::endl;
				}
#endif 
			}

// shifting iterator of start of new segment
			for ( unsigned int i = 0 ; ( ( i < moving_segment_shift ) && ( iter_data != end_iter ) ) ;  iter_data ++, i ++ );
		}
// averaging data
		for ( unsigned int i = 0 ; i < _q_th_order_fluctuation . size ( ) ; i ++ )
		{
			if ( _q_th_order_fluctuation [ i ] != 0.0 )
			{
				double data = fluctuation_function [ i ] / count_segments;
				* ( _DFAvariance_result [ i ] ) = pow ( data , 1.0 / _q_th_order_fluctuation [ i ]  );
			}
			else
			{
				double data = fluctuation_function [ i ] / ( 2 * count_segments );
				* ( _DFAvariance_result [ i ] ) = exp ( data );
			}
		}

#ifdef DEBUG
		{
			boost::posix_time::ptime time_of_log ( boost::posix_time::second_clock::local_time() );
			std::cerr << time_of_log << " Thread finished with segment size " << _segment_size << std::endl;
		}
#endif

		std::cerr << "Segment size " << _segment_size << " finished." << std::endl;
	}

};

void enqueue_thread_for_cleanup ( )
{
	boost::mutex::scoped_lock lock(add_jobs_mutex);

	finished_jobs.push_back ( boost::this_thread::get_id() );
// 		std::cout << "Huston, we have a problem with thread" << std::endl;
	add_jobs_CV.notify_one ();
}

/// Segment polynomial MF-DFA worker 
class SegmentDFAWorker : public SegmentWorker
#ifdef QTHREADS
, public QThread
#endif
{
protected:
	std::shared_ptr<std::vector<double *> > _DFAvariance_result;
	unsigned int _order_of_method;
#ifdef MICROTHREADING
	int _microthreads;
#endif
	void ( * _fluctuation_analyzer ) ( std::vector<double> & y_data , unsigned int order_of_method , std::list<double> & coefficients , double , math_functions::real_matrix * );

	double ( * _polynomial ) ( double x , std::list<double> & coefficients , double domain );
	bool _old_algorithm;

public:
	SegmentDFAWorker ( std::list<double> & data , std::vector<double> & q_th_order_fluctuation , std::shared_ptr< std::vector<double *> > DFAvariance_result , std::ofstream & output_stream , unsigned int order_of_method , unsigned int segment_size , bool moving_segment_adjustment = false , unsigned int moving_segment_shift = 2 , bool sample_variance = false , void ( * fluctuation_analyzer ) ( std::vector<double> &  , unsigned int , std::list<double> & , double,  math_functions::real_matrix * ) = &statistical_functions::PolynomialRegression<std::vector<double>,std::vector<double>::iterator> , double ( * polynomial ) ( double  , std::list<double> & , double  ) = &statistical_functions::Polynomial
#ifdef MICROTHREADING
, int microthreads
#endif
, bool old_algorithm = false )
	: SegmentWorker ( data , q_th_order_fluctuation , output_stream ,  segment_size , moving_segment_adjustment ,  moving_segment_shift , sample_variance ) , _DFAvariance_result ( DFAvariance_result ) , _order_of_method ( order_of_method ) , _fluctuation_analyzer ( fluctuation_analyzer )  , _polynomial ( polynomial )
#ifdef MICROTHREADING
, _microthreads ( microthreads )
#endif
, _old_algorithm ( old_algorithm )
	{
	}

#ifdef QTHREADS
	void run()
	{
		(* this) ( );
	}
#endif

	void operator () ()
	{
#ifndef MICROTHREADING
		if ( _old_algorithm == true )
		{
			std::list< std::vector<double> * > segments;
			std::vector<double> * new_segment = new std::vector<double> ( _segment_size );

			unsigned int moving_segment_shift = _segment_size;
			unsigned int counter = 0;
			unsigned int total_counter = 0;

			if ( _moving_segment_adjustment == true )
			{
				if ( _moving_segment_shift <= _segment_size )
				{
					moving_segment_shift = _moving_segment_shift;
				}
				else
				{
					moving_segment_shift = _segment_size;
				}
			}

	#ifdef DEBUG
			{
				boost::posix_time::ptime time_of_log ( boost::posix_time::second_clock::local_time() );
				std::cerr << time_of_log << " Thread start" << std::endl;
			}
	#endif

	// preparation of cached matrices for regression procedure
			math_functions::real_matrix transformation_y_data_cached;
			PrepareMatricesForRegression ( _segment_size  , _order_of_method , transformation_y_data_cached );

			// preparation of segments starting from first data to last
			for ( std::list<double>::iterator iter_data = _data.begin() ; iter_data != _data.end() ; iter_data ++ )
			{
				double actual_data = ( * iter_data );
				total_counter++;

				if ( counter < _segment_size )
				{
					( * new_segment ) [ counter ] = actual_data;
					counter ++;
				}
				else
				{
					segments . push_back( new_segment );

					if ( _data.size() - total_counter >= moving_segment_shift /*_segment_size*/ )
					{
						counter = 0;
						std::vector<double> * old_segment = new_segment;
						new_segment = new std::vector<double> ( _segment_size );

						for ( unsigned int i = 0 ; i < _segment_size - moving_segment_shift ; i ++ )
						{
							( * new_segment ) [ i ] = ( * old_segment ) [ i + moving_segment_shift ];
						}

						counter = _segment_size - moving_segment_shift;
						( * new_segment ) [ counter ] = actual_data;
						counter ++;
					}
					else
					{
						break;
					}
				}
			}

			if ( _moving_segment_adjustment == false )
			{
			// preparation of segments starting from last data to first
				counter = 0;
				total_counter = 0;
				new_segment = new std::vector<double> ( _segment_size );
				for ( std::list<double>::reverse_iterator iter_data = _data.rbegin() ; iter_data != _data.rend() ; iter_data ++ )
				{
					double actual_data = ( * iter_data );
					total_counter++;

					if ( counter < _segment_size )
					{
						( * new_segment ) [ _segment_size - 1 - counter] = actual_data;
						counter ++;
					}
					else
					{
						segments . push_back( new_segment );

						if ( _data . size() - total_counter >= moving_segment_shift /*_segment_size*/ )
						{
							counter = 0;
							std::vector<double> * old_segment = new_segment;
							new_segment = new std::vector<double> ( _segment_size );

							for ( unsigned int i = 0 ; i < _segment_size - moving_segment_shift ; i ++ )
							{
								( * new_segment ) [ i + moving_segment_shift ] = ( * old_segment ) [ i ];
							}

							counter = _segment_size - moving_segment_shift;
							( * new_segment ) [ _segment_size - 1 - counter] = actual_data;
							counter ++;
						}
						else
						{
							break;
						}
					}
				}
			}

	#ifdef DEBUG
			{
				boost::posix_time::ptime time_of_log ( boost::posix_time::second_clock::local_time() );
				std::cerr << time_of_log << " Segments ready" << std::endl;
			}
	#endif

			total_counter = 0;
			for ( std::list< std::vector<double> * >::iterator segment_iter = segments.begin() ; segment_iter != segments.end() ; ++ segment_iter )
			{
				std::vector<double> * actual_segment = ( * segment_iter );
				std::list<double> coefficients;

				_fluctuation_analyzer ( * actual_segment, _order_of_method , coefficients , _segment_size , & transformation_y_data_cached );
	// 			PolynomialRegression<std::vector<double>,std::vector<double>::iterator> ( * actual_segment, _order_of_method , coefficients );

				unsigned int count_i = 0;
				for ( std::vector<double>::iterator iter = actual_segment ->begin() ; iter != actual_segment->end() ; ++ iter )
				{
					double actual_value = ( * iter );

					double approximation = _polynomial ( count_i , coefficients , _segment_size );
	// 				double approximation = Polynomial ( count_i , coefficients );

					( * iter ) = actual_value - approximation;
	#ifdef DEBUG
					{
						boost::posix_time::ptime time_of_log ( boost::posix_time::second_clock::local_time() );
						math_model::debug::DebugLog ( 6 ) << time_of_log << " " << _segment_size << " " << total_counter << " " << count_i << " " << actual_value << " " << approximation << " " << actual_value - approximation << std::endl;
					}
	#endif
					if ( _output_stream.is_open ( ) )
					{
						_output_stream << _segment_size << " " << total_counter << " " << count_i << " " << actual_value << " " << approximation << " " << actual_value - approximation << std::endl;
					}
					count_i ++;
					total_counter ++;
				}
			}

	#ifdef DEBUG
			{
				boost::posix_time::ptime time_of_log ( boost::posix_time::second_clock::local_time() );
				std::cerr << time_of_log << " Regression in segments finished" << std::endl;
			}
	#endif

	// initial setup of positions of q-th order fluctiation function
			std::vector<double> fluctuation_function ( _q_th_order_fluctuation . size ( ) );
			for ( unsigned int i = 0 ; i < _q_th_order_fluctuation . size ( ) ; i ++ )
			{
				fluctuation_function [ i ] = 0;
			}

	// calculation of fluctuation function from segments
			for ( std::list< std::vector<double> * >::iterator segment_iter = segments.begin() ; segment_iter != segments.end() ; ++ segment_iter )
			{
				std::vector<double> * actual_segment = ( * segment_iter );

				double segment_variance = 0;
				for ( std::vector<double>::iterator iter = actual_segment ->begin() ; iter != actual_segment->end() ; ++ iter )
				{
					double actual_value = ( * iter );

					segment_variance += ( actual_value * actual_value );
				}
				if ( _sample_variance == false )
				{
					segment_variance /= _segment_size;
				}
				else
				{
					segment_variance /= ( _segment_size - 1 );
				}

				for ( unsigned int i = 0 ; i < _q_th_order_fluctuation . size ( ) ; i ++ )
				{
					if ( _q_th_order_fluctuation [ i ] != 0.0 )
					{
						fluctuation_function [ i ] += pow ( segment_variance , _q_th_order_fluctuation [ i ] / 2.0 );
					}
					else
					{
						fluctuation_function [ i ] += log ( segment_variance );
					}
				}
			}

			for ( unsigned int i = 0 ; i < _q_th_order_fluctuation . size ( ) ; i ++ )
			{
				if ( _q_th_order_fluctuation [ i ] != 0.0 )
				{
					double number_of_segments = segments.size();
					fluctuation_function [ i ] /= number_of_segments;
					double data = pow ( fluctuation_function [ i ] , 1.0 / _q_th_order_fluctuation [ i ] );

					* ( ( * _DFAvariance_result ) [ i ] ) = data;
	#ifdef DEBUG
					math_model::debug::DebugLog ( 6 ) << "Thread with segment size " << _segment_size << " calculated DFA" << _order_of_method << ": " << data << std::endl;
	#endif
	// 			std::cerr << "Thread with segment size " << _segment_size << " calculated DFA" << _order_of_method << ": " << data << std::endl;
				}
				else
				{
					double number_of_segments = segments.size();
					fluctuation_function [ i ] /= 2 * number_of_segments;
					double data = exp ( fluctuation_function [ i ] );

					* ( ( * _DFAvariance_result ) [ i ] ) = data;
	#ifdef DEBUG
					math_model::debug::DebugLog ( 6 ) << "Thread with segment size " << _segment_size << " calculated DFA" << _order_of_method << ": " << data << std::endl;
	#endif
				}
			}

	#ifdef DEBUG
			{
				boost::posix_time::ptime time_of_log ( boost::posix_time::second_clock::local_time() );
				std::cerr << time_of_log << " Calculation finished" << std::endl;
			}
	#endif

	// clean up structures prepared
			for ( std::list< std::vector<double> * >::iterator segment_iter = segments.begin() ; segment_iter != segments.end() ; ++ segment_iter )
			{
				std::vector<double> * actual_list = ( * segment_iter );
				delete actual_list;
			}
			segments.clear();
		}
		else
		{
// register a cleanup function
			boost::this_thread::at_thread_exit ( enqueue_thread_for_cleanup );

// new algorithm calculation that does not make list of all segments
			unsigned int moving_segment_shift = _segment_size;

			if ( _moving_segment_adjustment == true )
			{
				if ( _moving_segment_shift <= _segment_size )
				{
					moving_segment_shift = _moving_segment_shift;
				}
				else
				{
					moving_segment_shift = _segment_size;
				}
			}

	// initial setup of positions of q-th order fluctiation function
			std::vector<double> fluctuation_function ( _q_th_order_fluctuation . size ( ) );
			for ( unsigned int i = 0 ; i < _q_th_order_fluctuation . size ( ) ; i ++ )
			{
				fluctuation_function [ i ] = 0;
			}

#ifdef DEBUG
			{
				boost::posix_time::ptime time_of_log ( boost::posix_time::second_clock::local_time() );
				std::cerr << time_of_log << " Thread start with segment size " << _segment_size << std::endl;
			}
#endif

// preparation of cached matrices for regression procedure
			math_functions::real_matrix transformation_y_data_cached;
			PrepareMatricesForRegression ( _segment_size  , _order_of_method , transformation_y_data_cached );

// preparation of segments starting from first data to last
			unsigned int count_segments = 0;
			std::list<double>::iterator end_iter = _data.end();
			for ( std::list<double>::iterator iter_data = _data.begin() ; iter_data != end_iter ; )
			{
				std::vector<double> new_segment ( _segment_size );
				unsigned int count_elements_of_segment = 0;

				for (  std::list<double>::iterator iter_data_segment = iter_data ;  ( count_elements_of_segment < _segment_size ) && (  iter_data_segment != end_iter ) ; iter_data_segment ++ , count_elements_of_segment ++ )
				{
					double actual_data = ( * iter_data_segment );
					new_segment [ count_elements_of_segment ] = actual_data;
				}

				if ( count_elements_of_segment /*new_segment.size() */ == _segment_size )
				{
#ifdef DEBUG
					std::stringstream data_for_debug;
					for ( unsigned int i = 0 ; i < _segment_size ; i ++ )
					{
						data_for_debug << new_segment [ i ] << " ";
					}
					std::string new_segment_data_for_debug = data_for_debug.str();
					boost::posix_time::ptime time_of_log ( boost::posix_time::second_clock::local_time() );
					math_model::debug::DebugLog ( 6 ) << time_of_log << " segment: " << new_segment_data_for_debug << std::endl;
#endif

					std::list<double> coefficients ( _order_of_method );
					_fluctuation_analyzer ( new_segment, _order_of_method , coefficients , _segment_size , & transformation_y_data_cached );
#ifdef DEBUG
					std::stringstream data_for_debug2;
					for ( std::list<double>::iterator iter_coefficients = coefficients.begin() ; iter_coefficients != coefficients.end() ; iter_coefficients ++ )
					{
						data_for_debug2 << ( * iter_coefficients ) << " ";
					}
					std::string new_coefficient_data = data_for_debug2.str();
					boost::posix_time::ptime time_of_log2 ( boost::posix_time::second_clock::local_time() );
					math_model::debug::DebugLog ( 6 ) << time_of_log2 << " coefficients: " << new_coefficient_data << std::endl;
#endif

// calculation of regression
					unsigned int counter_of_position_in_segment = 0;
					std::vector<double>::iterator end_segment_iter = new_segment.end();
					for ( std::vector<double>::iterator iter = new_segment.begin() ; iter != end_segment_iter ; ++ iter )
					{
						double actual_value = ( * iter );
						double approximation = _polynomial ( counter_of_position_in_segment , coefficients , _segment_size );

						( * iter ) = actual_value - approximation;
#ifdef DEBUG
						boost::posix_time::ptime time_of_log ( boost::posix_time::second_clock::local_time() );
						math_model::debug::DebugLog ( 6 ) << time_of_log << " " << _segment_size << " " << count_segments << " " << counter_of_position_in_segment << " " << actual_value << " " << approximation << " " << actual_value - approximation << std::endl;
#endif
						if ( _output_stream.is_open ( ) )
						{
							_output_stream << _segment_size << " " << count_segments << " " << counter_of_position_in_segment << " " << actual_value << " " << approximation << " " << actual_value - approximation << std::endl;
						}
						counter_of_position_in_segment ++;
					}
					count_segments ++;

	#ifdef DEBUG
					{
						
						time_of_log = boost::posix_time::second_clock::local_time();
						std::cerr << time_of_log << " Regression in segments finished with segment size " << _segment_size << std::endl;
					}
	#endif

// calculation of fluctuation function
					double segment_variance = 0;
					for ( std::vector<double>::iterator iter = new_segment.begin() ; iter != end_segment_iter ; ++ iter )
					{
						double actual_value = ( * iter );

						segment_variance += ( actual_value * actual_value );
					}

// normalization to get sampled variance or variance
					if ( _sample_variance == false )
					{
						segment_variance /= _segment_size;
					}
					else
					{
						segment_variance /= ( _segment_size - 1 );
					}

					for ( unsigned int i = 0 ; i < _q_th_order_fluctuation . size ( ) ; i ++ )
					{
						if ( _q_th_order_fluctuation [ i ] != 0.0 )
						{
							fluctuation_function [ i ] += pow ( segment_variance , _q_th_order_fluctuation [ i ] / 2.0 );
// 							fluctuation_function [ i ] += boost::math::powm1 ( segment_variance , _q_th_order_fluctuation [ i ] / 2.0 ) + 1;
						}
						else
						{
							fluctuation_function [ i ] += log ( fabs ( segment_variance ) );
// 							fluctuation_function [ i ] += boost::math::log1p ( fabs ( segment_variance - 1.0 ) );
						}
					}
				}

// shifting iterator of start of new segment
				for ( unsigned int i = 0 ; ( ( i < moving_segment_shift ) && ( iter_data != end_iter ) ) ;  iter_data ++, i ++ );
			}
// averaging data
			for ( unsigned int i = 0 ; i < _q_th_order_fluctuation . size ( ) ; i ++ )
			{
				if ( _q_th_order_fluctuation [ i ] != 0.0 )
				{
					double data = fluctuation_function [ i ] / count_segments;
					* ( ( * _DFAvariance_result ) [ i ] ) = pow ( data , 1.0 / _q_th_order_fluctuation [ i ] );
// 					* ( ( * _DFAvariance_result ) [ i ] ) = boost::math::powm1 ( data , 1.0 / _q_th_order_fluctuation [ i ] ) + 1;
				}
				else
				{
					double data = fluctuation_function [ i ] / ( 2 * count_segments );
					* ( ( * _DFAvariance_result ) [ i ] ) = exp ( data );
// 					* ( ( * _DFAvariance_result ) [ i ] ) = boost::math::expm1 ( data ) + 1;
				}

				if ( std::isnan ( * ( ( * _DFAvariance_result ) [ i ] ) ) )
				{
					std::cerr << "Problem with result. Segment: " << _segment_size << " q: " << _q_th_order_fluctuation [ i ] << std::endl;
				}

			}
		}

#ifdef DEBUG
		boost::posix_time::ptime time_of_log ( boost::posix_time::second_clock::local_time() );
		std::cerr << time_of_log << " Thread finished with segment size " << _segment_size << std::endl;
#endif

#else
// algorithm for microthreading

		finish_thread_mutex.lock();
		finish_threads = false;
		finish_thread_mutex.unlock();

		work_list_iterator_mutex.lock();
		work_list_iterator = segments.begin ( );
		end_work_iterator = segments.end ( );

		std::list<double> * new_segment = new std::list<double>;
		work_list_iterator_mutex.unlock();

		unsigned int counter = 0;
		unsigned int total_counter = 0;

		boost::thread_group threads;
		for ( int i = 0 ; i < _microthreads ; i ++ )
		{
			PolynomialRegressionWorker * new_thread = new PolynomialRegressionWorker ( i , segments , _output_stream , _order_of_method , _segment_size );
			threads.add_thread ( new boost::thread ( * new_thread ) );
		}

		// preparation of segments starting from first data to last
		for ( std::list<double>::iterator iter_data = _data.begin() ; iter_data != _data.end() ; iter_data ++ )
		{
			double actual_data = ( * iter_data );
			total_counter++;

			if ( counter < _segment_size )
			{
				new_segment->push_back( actual_data );
				counter ++;
			}
			else
			{
				counter = 0;
				if ( _data.size() - total_counter >= _segment_size )
				{

					work_list_iterator_mutex.lock();

					segments . push_back( new_segment );

					if ( work_list_iterator == end_work_iterator )
					{
						work_list_iterator = segments.begin ( );
					}
					work_list_iterator_mutex.unlock();

					new_segment = new std::list<double>;
					new_segment->push_back( actual_data );
					counter ++;
				}
				else
				{
					break;
				}
			}
		}

		// preparation of segments starting from last data to first
		counter = 0;
		total_counter = 0;
		new_segment = new std::list<double>;
		for ( std::list<double>::reverse_iterator iter_data = _data.rbegin() ; iter_data != _data.rend() ; iter_data ++ )
		{
			double actual_data = ( * iter_data );
			total_counter++;

			if ( counter < _segment_size )
			{
				new_segment->push_front( actual_data );
				counter ++;
			}
			else
			{
				counter = 0;
				if ( _data.size() - total_counter >= _segment_size )
				{

					work_list_iterator_mutex.lock();

					segments . push_back( new_segment );

					if ( work_list_iterator == end_work_iterator )
					{
						work_list_iterator = segments.begin ( );
					}
					work_list_iterator_mutex.unlock();

					new_segment = new std::list<double>;
					new_segment->push_back( actual_data );
					counter ++;
				}
				else
				{
					break;
				}
			}
		}

		finish_thread_mutex.lock();
		finish_threads = true;
		finish_thread_mutex.unlock();

		threads.join_all();

		double fluctuation_function = 0;
		for ( std::list< std::list<double> * >::iterator segment_iter = segments.begin() ; segment_iter != segments.end() ; segment_++ iter )
		{
			std::list<double> * actual_segment = ( * segment_iter );

			unsigned int count_i = 0;
			double segment_variance = 0;
			for ( std::list<double>::iterator iter = actual_segment ->begin() ; iter != actual_segment->end() ; ++ iter )
			{
				double actual_value = ( * iter );

				segment_variance += ( actual_value * actual_value );
			}
			if ( _sample_variance == false )
			{
				segment_variance /= _segment_size;
			}
			else
			{
				segment_variance /= ( _segment_size - 1 );
			}

			if ( _q_th_order_fluctuation != 0.0 )
			{
				fluctuation_function += pow ( segment_variance , _q_th_order_fluctuation / 2.0 );
			}
			else
			{
				fluctuation_function += log ( segment_variance );
			}
		}
		if ( _q_th_order_fluctuation != 0.0 )
		{
			fluctuation_function /= segments.size();
			double data = pow ( fluctuation_function , 1.0 / _q_th_order_fluctuation );

			_DFAvariance_result = data;
		}
		else
		{
			fluctuation_function /= 2 * segments.size();
			double data = exp ( fluctuation_function );

			_DFAvariance_result = data;
		}

// clean up structures prepared
		for ( std::list< std::list<double> * >::iterator segment_iter = segments.begin() ; segment_iter != segments.end() ; segment_++ iter )
		{
			std::list<double> * actual_list = ( * segment_iter );
			delete actual_list;
		}
		segments.clear();
#endif

// 		std::cerr << "Segment size " << _segment_size << " finished." << std::endl;
	}
};


/// Segment polynomial MF-DCCA worker 
class SegmentDCCAWorker : public SegmentWorker
#ifdef QTHREADS
, public QThread
#endif
{
protected:
	std::shared_ptr<std::vector< std::array < double * , 4 > > > _DCCA_result;
	std::list<double> & _data2;
	unsigned int _order_of_method;
	void ( * _fluctuation_analyzer ) ( std::vector<double> & y_data , unsigned int order_of_method , std::list<double> & coefficients , double , math_functions::real_matrix * );

	double ( * _polynomial ) ( double x , std::list<double> & coefficients , double domain );

public:
	SegmentDCCAWorker ( std::list<double> & data1 , std::list<double> & data2 , std::vector<double> & q_th_order_fluctuation , std::shared_ptr<std::vector< std::array < double * , 4 > > > DCCA_result , std::ofstream & output_stream , unsigned int order_of_method , unsigned int segment_size , bool moving_segment_adjustment = false , unsigned int moving_segment_shift = 2 , bool sample_variance = false , void ( * fluctuation_analyzer ) ( std::vector<double> &  , unsigned int , std::list<double> & , double,  math_functions::real_matrix * ) = &statistical_functions::PolynomialRegression<std::vector<double>,std::vector<double>::iterator> , double ( * polynomial ) ( double  , std::list<double> & , double  ) = &statistical_functions::Polynomial )
	: SegmentWorker ( data1 , q_th_order_fluctuation , output_stream ,  segment_size , moving_segment_adjustment ,  moving_segment_shift , sample_variance ) , _DCCA_result ( DCCA_result ) , _data2 ( data2 ) , _order_of_method ( order_of_method ) , _fluctuation_analyzer ( fluctuation_analyzer ) , _polynomial ( polynomial )
	{
		if ( _data.size() != _data2.size() )
		{
			std::cout << "SegmentDCCAWorker::SegmentDCCAWorker: Datasets of interest does not have the same size." << std::endl;
			exit ( EXIT_FAILURE );
		}
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
		boost::this_thread::at_thread_exit ( enqueue_thread_for_cleanup );

// new algorithm calculation that does not make list of all segments
		unsigned int moving_segment_shift = _segment_size;

		if ( _moving_segment_adjustment == true )
		{
			if ( _moving_segment_shift <= _segment_size )
			{
				moving_segment_shift = _moving_segment_shift;
			}
			else
			{
				moving_segment_shift = _segment_size;
			}
		}

// initial setup of positions of q-th order fluctiation function
		std::vector<double> fluctuation_function ( _q_th_order_fluctuation . size ( ) );
		std::vector<double> fluctuation_function2 ( _q_th_order_fluctuation . size ( ) );
		std::vector<double> crosscorrelation_function ( _q_th_order_fluctuation . size ( ) );
		std::vector<double> abs_crosscorrelation_function ( _q_th_order_fluctuation . size ( ) );
		for ( unsigned int i = 0 ; i < _q_th_order_fluctuation . size ( ) ; i ++ )
		{
			fluctuation_function [ i ] = 0;
			fluctuation_function2 [ i ] = 0;
			crosscorrelation_function [ i ] = 0;
			abs_crosscorrelation_function [ i ] = 0;
		}

#ifdef DEBUG
		{
			boost::mutex::scoped_lock lock( text_output_mutex );
			boost::posix_time::ptime time_of_log ( boost::posix_time::second_clock::local_time() );
			std::cerr << time_of_log << " Thread start with segment size " << _segment_size << std::endl;
		}
#endif

// preparation of cached matrices for regression procedure
		math_functions::real_matrix transformation_y_data_cached;
		PrepareMatricesForRegression ( _segment_size  , _order_of_method , transformation_y_data_cached );

// preparation of segments starting from first data to last
		unsigned int count_segments = 0;
		std::list<double>::iterator end_iter = _data.end();
		std::list<double>::iterator iter_data2 = _data2.begin();
		for ( std::list<double>::iterator iter_data = _data.begin() ; iter_data != end_iter ; )
		{
			std::vector<double> new_segment ( _segment_size );
			std::vector<double> new_segment2 ( _segment_size );
			unsigned int count_elements_of_segment = 0;

			std::list<double>::iterator iter_data_segment2 = iter_data2;
			for (  std::list<double>::iterator iter_data_segment = iter_data ;  ( count_elements_of_segment < _segment_size ) && (  iter_data_segment != end_iter ) ; iter_data_segment ++ , count_elements_of_segment ++ , iter_data_segment2 ++ )
			{
				double actual_data = ( * iter_data_segment );
				double actual_data2 = ( * iter_data_segment2 );
				new_segment [ count_elements_of_segment ] = actual_data;
				new_segment2 [ count_elements_of_segment ] = actual_data2;
			}

			if ( count_elements_of_segment == _segment_size )
			{
#ifdef DEBUG
				{
					boost::mutex::scoped_lock lock( text_output_mutex );
					std::stringstream data_for_debug;
					for ( unsigned int i = 0 ; i < _segment_size ; i ++ )
					{
						data_for_debug << new_segment [ i ] << " ";
					}
					std::string new_segment_data_for_debug = data_for_debug.str();
					boost::posix_time::ptime time_of_log ( boost::posix_time::second_clock::local_time() );
					math_model::debug::DebugLog ( 6 ) << time_of_log << " segment: " << new_segment_data_for_debug << std::endl;
				}
#endif

				std::list<double> coefficients ( _order_of_method );
				std::list<double> coefficients2 ( _order_of_method );
				_fluctuation_analyzer ( new_segment, _order_of_method , coefficients , _segment_size , & transformation_y_data_cached );
				_fluctuation_analyzer ( new_segment2, _order_of_method , coefficients2 , _segment_size , & transformation_y_data_cached );
#ifdef DEBUG
				{
					boost::mutex::scoped_lock lock( text_output_mutex );
					std::stringstream data_for_debug2;
					for ( std::list<double>::iterator iter_coefficients = coefficients.begin() ; iter_coefficients != coefficients.end() ; iter_coefficients ++ )
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
				std::vector<double>::iterator end_segment_iter = new_segment.end();
				std::vector<double>::iterator iter2 = new_segment2.begin();
				for ( std::vector<double>::iterator iter = new_segment.begin() ; iter != end_segment_iter ; ++ iter , ++ iter2 )
				{
					double actual_value = ( * iter );
					double actual_value2 = ( * iter2 );
					double approximation = _polynomial ( counter_of_position_in_segment , coefficients , _segment_size );
					double approximation2 = _polynomial ( counter_of_position_in_segment , coefficients2 , _segment_size );

					( * iter ) = actual_value - approximation;
					( * iter2 ) = actual_value2 - approximation2;
#ifdef DEBUG
					{
						boost::mutex::scoped_lock lock( text_output_mutex );
						boost::posix_time::ptime time_of_log ( boost::posix_time::second_clock::local_time() );
						math_model::debug::DebugLog ( 6 ) << time_of_log << " " << _segment_size << " " << count_segments << " " << counter_of_position_in_segment << " " << actual_value << " " << approximation << " " << actual_value - approximation << std::endl;
					}
#endif
					if ( _output_stream.is_open ( ) )
					{
						_output_stream << _segment_size << " " << count_segments << " " << counter_of_position_in_segment << " " << actual_value << " " << approximation << " " << actual_value - approximation << std::endl;
					}
					counter_of_position_in_segment ++;
				}
				count_segments ++;

#ifdef DEBUG
				{
					boost::mutex::scoped_lock lock( text_output_mutex );
					boost::posix_time::ptime time_of_log = boost::posix_time::second_clock::local_time();
					std::cerr << time_of_log << " Regression in segments finished with segment size " << _segment_size << std::endl;
				}
#endif

// calculation of fluctuation function
				double segment_variance = 0;
				double segment_variance2 = 0;
				double segment_covariance = 0;
				double segment_abs_covariance = 0;

				iter2 = new_segment2.begin();
				for ( auto iter = new_segment.begin() ; iter != end_segment_iter ; ++ iter , ++ iter2 )
				{
					double actual_value = ( * iter );
					double actual_value2 = ( * iter2 );

					segment_variance += ( actual_value * actual_value );
					segment_variance2 += ( actual_value2 * actual_value2 );
					segment_covariance += ( actual_value * actual_value2 );
					segment_abs_covariance += ( fabs ( actual_value ) * fabs ( actual_value2 ) );
				}

// normalization to get sampled variance or variance
				if ( _sample_variance == false )
				{
					segment_variance /= _segment_size;
					segment_variance2 /= _segment_size;
					segment_covariance /= _segment_size;
					segment_abs_covariance /= _segment_size;
				}
				else
				{
					segment_variance /= ( _segment_size - 1 );
					segment_variance2 /= ( _segment_size - 1 );
					segment_covariance /= ( _segment_size - 1 );
					segment_abs_covariance /= ( _segment_size - 1 );
				}

				for ( unsigned int i = 0 ; i < _q_th_order_fluctuation . size ( ) ; i ++ )
				{
					if ( _q_th_order_fluctuation [ i ] != 0.0 )
					{
						fluctuation_function [ i ] += pow ( segment_variance , _q_th_order_fluctuation [ i ] / 2.0 );
						fluctuation_function2 [ i ] += pow ( segment_variance2 , _q_th_order_fluctuation [ i ] / 2.0 );
						crosscorrelation_function [ i ] += pow ( fabs ( segment_covariance ) , _q_th_order_fluctuation [ i ] / 2.0 );
						abs_crosscorrelation_function [ i ] += pow ( fabs ( segment_abs_covariance ) , _q_th_order_fluctuation [ i ] / 2.0 );
					}
					else
					{
						fluctuation_function [ i ] += log ( segment_variance );
						fluctuation_function2 [ i ] += log ( segment_variance2 );
						crosscorrelation_function [ i ] += log ( fabs ( segment_covariance ) );
						abs_crosscorrelation_function [ i ] += log ( fabs ( segment_abs_covariance ) );
					}
				}
			}

// shifting iterator of start of new segment
			for ( unsigned int i = 0 ; ( ( i < moving_segment_shift ) && ( iter_data != end_iter ) ) ;  iter_data ++ , iter_data2 ++ , i ++ );
		}

// averaging data
		for ( unsigned int i = 0 ; i < _q_th_order_fluctuation . size ( ) ; i ++ )
		{
			if ( ( ( * _DCCA_result ) [ i ] [ 0 ] == nullptr ) || ( ( * _DCCA_result ) [ i ] [ 1 ] == nullptr ) || ( ( * _DCCA_result ) [ i ] [ 2 ] == nullptr ) || ( ( * _DCCA_result ) [ i ] [ 3 ] == nullptr ) )
			{
				std::cerr << "SegmentDCCAWorker::operator: Error: Nullptr encoutered" << std::endl;
				std::cerr << "Parameters are: Order: " << _order_of_method << " , Segment size: " << _segment_size << " , i: " << i << " , affiliated q: " << _q_th_order_fluctuation [ i ] << " , Pointers are: " << ( * _DCCA_result ) [ i ] [ 0 ] << " , " << ( * _DCCA_result ) [ i ] [ 1 ] << " , " << ( * _DCCA_result ) [ i ] [ 2 ] << " , " << ( * _DCCA_result ) [ i ] [ 3 ] << std::endl;
				exit ( EXIT_FAILURE );
			}

			if ( _q_th_order_fluctuation [ i ] != 0.0 )
			{
				double data1 = fluctuation_function [ i ] / count_segments;
				double data2 = fluctuation_function2 [ i ] / count_segments;
				double data3 = crosscorrelation_function [ i ] / count_segments;
				double data4 = abs_crosscorrelation_function [ i ] / count_segments;

				( * ( * _DCCA_result ) [ i ] [ 0 ] ) = pow ( data1 , 1.0 / _q_th_order_fluctuation [ i ] );
				( * ( * _DCCA_result ) [ i ] [ 1 ] ) = pow ( data2 , 1.0 / _q_th_order_fluctuation [ i ] );
				( * ( * _DCCA_result ) [ i ] [ 2 ] ) = pow ( fabs ( data3 ) , 1.0 / _q_th_order_fluctuation [ i ]  );
				( * ( * _DCCA_result ) [ i ] [ 3 ] ) = pow ( fabs ( data4 ) , 1.0 / _q_th_order_fluctuation [ i ]  );
			}
			else
			{
				double data1 = fluctuation_function [ i ] / ( 2 * count_segments );
				double data2 = fluctuation_function2 [ i ] / ( 2 * count_segments );
				double data3 = crosscorrelation_function [ i ] / ( 2 * count_segments );
				double data4 = abs_crosscorrelation_function [ i ] / ( 2 * count_segments );

				( * ( * _DCCA_result ) [ i ] [ 0 ] ) = exp ( data1 );
				( * ( * _DCCA_result ) [ i ] [ 1 ] ) = exp ( data2 );
				( * ( * _DCCA_result ) [ i ] [ 2 ] ) = exp ( data3 );
				( * ( * _DCCA_result ) [ i ] [ 3 ] ) = exp ( data4 );
			}
		}

#ifdef DEBUG
		{
			boost::mutex::scoped_lock lock( text_output_mutex );
			boost::posix_time::ptime time_of_log = boost::posix_time::second_clock::local_time();
			std::cerr << time_of_log << " Thread finished with segment size " << _segment_size << std::endl;
		}
#endif
	}
};

void IntegrationOfDataset ( std::list<double> & data )
{
// determination  of mean
	double complete_sum = 0;
	for ( std::list<double>::iterator iter = data.begin() ; iter != data.end() ; ++ iter )
	{
		double & value = ( * iter );
		complete_sum += value;
	}
	double mean = complete_sum / data.size();

// integration of dataset
	double sum = 0;
	for ( std::list<double>::iterator iter = data.begin() ; iter != data.end() ; ++ iter )
	{
		double & value = ( * iter );
		sum += ( value - mean );
		value = sum;
	}
}


void DetrendedFluctuationAnalysisLinearlyScaledDomain ( std::list<double> & data , std::vector<unsigned int> & segment_sizes , std::map<std::tuple<double, unsigned int, unsigned int> , double > & DFA , std::vector<double> & q_th_order_fluctuation , unsigned int start_index , unsigned int end_index , unsigned segment_step , const std::vector<long> & order_of_method , std::string output_analysis , bool moving_segment_adjustment , unsigned int moving_segment_shift , bool sample_variance , bool double_integration
#ifdef THREADING
, unsigned int number_threads
#ifdef MICROTHREADING
, unsigned int microthreads
#endif
#endif
, void ( * fluctuation_analyzer ) ( std::vector<double> &  , unsigned int , std::list<double> &  , double , math_functions::real_matrix * ) , double ( * polynomial ) ( double  , std::list<double> & , double  ) )
{
// loop over segment sizes that are going to be iterated
	unsigned int max_segment_size = end_index;
	unsigned int segment_size = start_index;

	while ( segment_size <= max_segment_size )
	{
		segment_sizes.push_back ( segment_size );
		segment_size += segment_step;
	}

	DetrendedFluctuationAnalysis ( data , segment_sizes , DFA , q_th_order_fluctuation , order_of_method , output_analysis , moving_segment_adjustment , moving_segment_shift , sample_variance , double_integration
#ifdef THREADING
, number_threads
#ifdef MICROTHREADING
, microthreads
#endif
#endif
, fluctuation_analyzer , polynomial );
}


void DetrendedFluctuationAnalysis ( std::list<double> & data , std::vector<unsigned int> & segment_sizes , std::map<std::tuple<double, unsigned int, unsigned int> , double > & DFA , std::vector<double> & q_th_order_fluctuation , const std::vector<long> & order_of_method , std::string output_analysis , bool moving_segment_adjustment , unsigned int moving_segment_shift , bool sample_variance , bool double_integration
#ifdef THREADING
, unsigned int number_threads
#ifdef MICROTHREADING
, unsigned int microthreads
#endif
#endif
, void ( * fluctuation_analyzer ) ( std::vector<double> &  , unsigned int , std::list<double> &  , double , math_functions::real_matrix * ) , double ( * polynomial ) ( double  , std::list<double> & , double  ) )
{
#ifdef DEBUG
	math_model::debug::OpenDebugStream ();
#endif

// open output file
	std::ofstream output_analysis_file;
	if ( ! output_analysis.empty() )
	{
		output_analysis_file.open ( output_analysis.c_str() );
	}

	std::list<double> integrated_dataset ( data );
	IntegrationOfDataset ( integrated_dataset );

	if ( double_integration == true )
	{
		IntegrationOfDataset ( integrated_dataset );
	}

	std::vector<unsigned int>::iterator actual_segment_iter = segment_sizes.begin ();
	std::vector<long>::const_iterator order_of_method_iter = order_of_method.cbegin ();

#ifdef THREADING
#ifdef QTHREADS
	std::list<SegmentDFAWorker *> DFA_threads;
#else
	boost::thread_group * DFA_threads = new boost::thread_group;
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
			std::cout << actual_time << " Start to submit job of the segment size " << ( * actual_segment_iter ) << std::endl;
		}

		const unsigned int number_of_q_elements = q_th_order_fluctuation.size();
		std::shared_ptr< std::vector< double * > > DFA_variance_results ( new std::vector< double * > ( number_of_q_elements ) );

		for ( unsigned int i = 0 ; i < number_of_q_elements ; i ++ )
		{
			DFA [ std::make_tuple ( q_th_order_fluctuation [ i ] , ( * order_of_method_iter ) , ( * actual_segment_iter ) ) ] = static_cast<std::map<std::tuple<double, unsigned int, unsigned int> , double >::mapped_type> ( 0.0 );

// allocation of data for results
			double & results = DFA [ std::make_tuple ( q_th_order_fluctuation [ i ] , ( * order_of_method_iter ) , ( * actual_segment_iter ) ) ];

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

		SegmentDFAWorker * thread_worker = new SegmentDFAWorker ( integrated_dataset , q_th_order_fluctuation , DFA_variance_results , output_analysis_file , ( * order_of_method_iter ) , ( * actual_segment_iter ) , moving_segment_adjustment , moving_segment_shift , sample_variance
#ifdef MICROTHREADING
 , microthreads
#endif
, fluctuation_analyzer , polynomial
		);

		{
			size_t number_of_threads;
			size_t number_finished_jobs;

			{
				boost::mutex::scoped_lock lock(add_jobs_mutex);
				boost::thread * new_thread = new boost::thread ( * thread_worker );

				DFA_threads -> add_thread ( new_thread );

				number_of_threads = DFA_threads->size ( );
				number_finished_jobs = finished_jobs.size ( );

				if ( ( number_of_threads - number_finished_jobs ) >= number_threads )
				{
					// wait for a thread to wake up
					add_jobs_CV.wait ( lock );
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
	finished_jobs.clear ();

#else

// adding list for results
	for ( unsigned int j = 0 ; j < q_th_order_fluctuation . size ( ) ; j ++ )
	{
		std::list<double> list_to_add;
		DFAvariance.push_back(list_to_add);
	}

// loop over segment size
	for ( int count_segments = 0 ; actual_segment_iter != segment_sizes.end () ; ++ actual_segment_iter , ++ count_segments)
	{
		std::list< std::vector<double> * > segments;
		std::vector<double> * new_segment = new std::vector<double> ( * actual_segment_iter );

		unsigned int counter = 0;
		unsigned int total_counter = 0;

		// preparation of segments starting from first data to last
		for ( std::list<double>::iterator iter_data = data.begin() ; iter_data != data.end() ; iter_data ++ )
		{
			double actual_data = ( * iter_data );
			total_counter++;

			if ( counter < * actual_segment_iter )
			{
				( * new_segment ) [ counter ] = actual_data;
				counter ++;
			}
			else
			{
				counter = 0;
				segments . push_back( new_segment );
				if ( data.size() - total_counter >= * actual_segment_iter )
				{
					new_segment = new std::vector<double> ( * actual_segment_iter );
					( * new_segment ) [ counter ] = actual_data;
					counter ++;
				}
				else
				{
					break;
				}
			}
		}

// preparation of segments starting from last data to first
		counter = 0;
		total_counter = 0;
		new_segment = new std::vector<double> ( * actual_segment_iter );
		for ( std::list<double>::reverse_iterator iter_data = data.rbegin() ; iter_data != data.rend() ; iter_data ++ )
		{
			double actual_data = ( * iter_data );
			total_counter++;

			if ( counter < * actual_segment_iter )
			{
				( * new_segment ) [ * actual_segment_iter - 1 - counter ] = actual_data;
				counter ++;
			}
			else
			{
				counter = 0;
				segments . push_back( new_segment );
				if ( data.size() - total_counter >= * actual_segment_iter )
				{
					new_segment = new std::vector <double> ( * actual_segment_iter );
					( * new_segment ) [ * actual_segment_iter - 1 - counter ] = actual_data;
					counter ++;
				}
				else
				{
					break;
				}
			}
		}

		total_counter = 0;
		for ( std::list< std::vector<double> * >::iterator segment_iter = segments.begin() ; segment_iter != segments.end() ; ++ segment_iter )
		{
			std::vector<double> * actual_segment = ( * segment_iter );
			std::list<double> coefficients;
			PolynomialRegression<std::vector<double>,std::vector<double>::iterator> ( * actual_segment, order_of_method , coefficients );

			unsigned int count_i = 0;
			for ( std::vector<double>::iterator iter = actual_segment ->begin() ; iter != actual_segment->end() ; ++ iter )
			{
				double actual_value = ( * iter );
				double approximation = Polynomial ( count_i , coefficients );
				( * iter ) = actual_value - approximation;
#ifdef DEBUG
				std::cout << * actual_segment_iter << " " << total_counter << " " << count_i << " " << actual_value << " " << approximation << " " << actual_value - approximation << std::endl;
#endif
				if ( output_analysis_file.is_open ( ) )
				{
					output_analysis_file << * actual_segment_iter << " " << total_counter << " " << count_i << " " << actual_value << " " << approximation << " " << actual_value - approximation << std::endl;
				}
				count_i ++;
				total_counter ++;
			}
		}

// allocation space for fluctuation fuction data and setup of initial values
		std::vector<double> fluctuation_function ( q_th_order_fluctuation.size() );
		for ( unsigned int j = 0 ; j < q_th_order_fluctuation . size ( ) ; j ++ )
		{
			fluctuation_function [ j ] = 0;
		}

		for ( std::list< std::vector<double> * >::iterator segment_iter = segments.begin() ; segment_iter != segments.end() ; ++ segment_iter )
		{
			std::vector<double> * actual_segment = ( * segment_iter );

			unsigned int count_i = 0;
			double segment_variance = 0;
			for ( std::vector<double>::iterator iter = actual_segment ->begin() ; iter != actual_segment->end() ; ++ iter )
			{
				double actual_value = ( * iter );

				segment_variance += ( actual_value * actual_value );
			}
			segment_variance /= * actual_segment_iter;

			for ( unsigned int i = 0 ; i < q_th_order_fluctuation . size ( ) ; i ++ )
			{
				if ( q_th_order_fluctuation [ i ] != 0.0 )
				{
					fluctuation_function [ i ] += pow ( segment_variance , q_th_order_fluctuation [ i ] / 2.0 );
				}
				else
				{
					fluctuation_function [ i ] += log ( segment_variance );
				}
			}
		}

		for ( unsigned int i = 0 ; i < q_th_order_fluctuation . size ( ) ; i ++ )
		{
			if ( q_th_order_fluctuation [ i ] != 0.0 )
			{
				fluctuation_function [ i ] /= segments.size();
				double data = pow ( fluctuation_function [ i ] , 1.0 / q_th_order_fluctuation [ i ] );

				DFAvariance [ count_segments ] .push_back ( data );
			}
			else
			{
				fluctuation_function [ i ] /= 2 * segments.size();
				double data = exp ( fluctuation_function [ i ] );

				DFAvariance [ count_segments ] .push_back ( data );
			}
		}

// clean up structures prepared
		for ( std::list< std::vector<double> * >::iterator segment_iter = segments.begin() ; segment_iter != segments.end() ; ++ segment_iter )
		{
			std::vector<double> * actual_list = ( * segment_iter );
			delete actual_list;
		}
		segments.clear();
	}
#endif
	if ( output_analysis_file.is_open() )
	{
		output_analysis_file.close ();
	}
}


void DetrendedCrossCorrelationAnalysis ( std::list<double> & data1 , std::list<double> & data2 , std::vector<unsigned int> & segment_sizes , std::map<std::tuple<double, unsigned int, unsigned int> , std::array<double,4> > & DCCA , std::vector<double> & q_th_order_fluctuation , std::vector<unsigned int> & order_of_method , std::string output_analysis , bool moving_segment_adjustment , unsigned int moving_segment_shift , bool sample_variance , bool double_integration
#ifdef THREADING
, unsigned int number_threads
#endif
, void ( * fluctuation_analyzer ) ( std::vector<double> &  , unsigned int , std::list<double> &  , double , math_functions::real_matrix * ) , double ( * polynomial ) ( double  , std::list<double> & , double  ) )
{
#ifdef DEBUG
	math_model::debug::OpenDebugStream ();
#endif

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
	std::list<double> integrated_dataset_1 ( data1 );
	IntegrationOfDataset ( integrated_dataset_1 );

	if ( double_integration == true )
	{
		IntegrationOfDataset ( integrated_dataset_1 );
	}

	// dataset 2
	std::list<double> integrated_dataset_2 ( data2 );
	IntegrationOfDataset ( integrated_dataset_2 );

	if ( double_integration == true )
	{
		IntegrationOfDataset ( integrated_dataset_2 );
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
		const unsigned int number_of_q_elements = q_th_order_fluctuation.size();
		std::shared_ptr<std::vector< std::array < double * , 4 > > > DCCA_variance_results ( new std::vector< std::array < double * , 4 > > ( number_of_q_elements ) );

		for ( unsigned int i = 0 ; i < number_of_q_elements ; i ++ )
		{
// loop over all fluctuations asked, F^{a}_q (s) , F^{b}_q (s) , XF^{a,b}_q (s) , CF^{a,b}_q (s)

			std::array < double , 4 > zero;
			zero.fill ( 0 );
			DCCA [ std::make_tuple ( q_th_order_fluctuation [ i ] , ( * order_of_method_iter ) , ( * actual_segment_iter ) ) ] = zero;

// allocation of data for results
			std::array < double , 4 > & results = DCCA [ std::make_tuple ( q_th_order_fluctuation [ i ] , ( * order_of_method_iter ) , ( * actual_segment_iter ) ) ];

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

		SegmentDCCAWorker * thread_worker = new SegmentDCCAWorker ( integrated_dataset_1 , integrated_dataset_2 , q_th_order_fluctuation , DCCA_variance_results , output_analysis_file , ( * order_of_method_iter ) , ( * actual_segment_iter ) , moving_segment_adjustment , moving_segment_shift , sample_variance
#ifdef MICROTHREADING
 , microthreads
#endif
, fluctuation_analyzer , polynomial
		);
#ifdef QTHREADS
		DCCA_threads . push_back ( thread_worker );
		thread_worker->start();
#else
		{
			boost::mutex::scoped_lock lock(add_jobs_mutex);
			boost::thread * new_thread = new boost::thread ( * thread_worker );
			DCCA_threads -> add_thread ( new_thread );
#endif
// increment job position
			++ actual_segment_iter;

			if ( ( actual_segment_iter == segment_sizes.end () ) && ( order_of_method_iter != order_of_method.end () ) )
			{
				std::cout << "We finished to enqueue jobs of the order " << ( * order_of_method_iter ) << std::endl;
				actual_segment_iter = segment_sizes.begin ();
				++ order_of_method_iter;
			}

#ifdef QTHREADS
			if ( DCCA_threads.size ( ) > number_threads )
			{
				SegmentDFAWorker * actual_thread = DCCA_threads.back();
				actual_thread->wait();
				DCCA_threads . pop_back ( );
				delete actual_thread;
			}
#else
			if ( DCCA_threads->size ( ) - finished_jobs.size ( ) > number_threads )
			{
// 			boost::posix_time::seconds workTime ( 10 );
// 			boost::this_thread::sleep(workTime);
				add_jobs_CV.wait ( lock );
// 			DCCA_threads -> join_all();
// 			delete DFA_threads;
// 			DFA_threads = new boost::thread_group;
			}
		}
	}
// 	std::cout << DCCA [ std::make_tuple ( 2. , 2 , 8 )] [ 0 ] << " " << DCCA [ std::make_tuple ( 2. , 2 , 8 )] [ 1 ] << " " << DCCA [ std::make_tuple ( 2. , 2 , 8 )] [ 2 ] << std::endl;
#endif

#ifdef QTHREADS
	while ( DCCA_threads.size ( ) )
	{
#else
	{
#endif

#ifdef QTHREADS
		SegmentDCCAWorker * actual_thread = DCCA_threads.back();
		actual_thread->wait();
		DCCA_threads . pop_back ( );
		delete actual_thread;
#else
// 		boost::posix_time::seconds workTime ( 10 );
// 		boost::this_thread::sleep(workTime);

		DCCA_threads -> join_all();
		delete DCCA_threads;
		finished_jobs.clear ();
#endif
	}
#else


// adding list for results
	for ( unsigned int j = 0 ; j < q_th_order_fluctuation . size ( ) ; j ++ )
	{
		std::list<double> list_to_add;
		DFAvariance.push_back(list_to_add);
	}

// loop over segment size
	for ( int count_segments = 0 ; actual_segment_iter != segment_sizes.end () ; ++ actual_segment_iter , ++ count_segments )
	{
		std::list< std::vector<double> * > segments;
		std::vector<double> * new_segment = new std::vector<double> ( * actual_segment_iter );

		unsigned int counter = 0;
		unsigned int total_counter = 0;

		// preparation of segments starting from first data to last
		for ( std::list<double>::iterator iter_data = data.begin() ; iter_data != data.end() ; iter_data ++ )
		{
			double actual_data = ( * iter_data );
			total_counter++;

			if ( counter < * actual_segment_iter )
			{
				( * new_segment ) [ counter ] = actual_data;
				counter ++;
			}
			else
			{
				counter = 0;
				segments . push_back( new_segment );
				if ( data.size() - total_counter >= * actual_segment_iter )
				{
					new_segment = new std::vector<double> ( * actual_segment_iter );
					( * new_segment ) [ counter ] = actual_data;
					counter ++;
				}
				else
				{
					break;
				}
			}
		}

// preparation of segments starting from last data to first
		counter = 0;
		total_counter = 0;
		new_segment = new std::vector<double> ( * actual_segment_iter );
		for ( std::list<double>::reverse_iterator iter_data = data.rbegin() ; iter_data != data.rend() ; iter_data ++ )
		{
			double actual_data = ( * iter_data );
			total_counter++;

			if ( counter < * actual_segment_iter )
			{
				( * new_segment ) [ * actual_segment_iter - 1 - counter ] = actual_data;
				counter ++;
			}
			else
			{
				counter = 0;
				segments . push_back( new_segment );
				if ( data.size() - total_counter >= * actual_segment_iter )
				{
					new_segment = new std::vector <double> ( * actual_segment_iter );
					( * new_segment ) [ * actual_segment_iter - 1 - counter ] = actual_data;
					counter ++;
				}
				else
				{
					break;
				}
			}
		}

		total_counter = 0;
		for ( std::list< std::vector<double> * >::iterator segment_iter = segments.begin() ; segment_iter != segments.end() ; ++ segment_iter )
		{
			std::vector<double> * actual_segment = ( * segment_iter );
			std::list<double> coefficients;
			PolynomialRegression<std::vector<double>,std::vector<double>::iterator> ( * actual_segment, order_of_method , coefficients );

			unsigned int count_i = 0;
			for ( std::vector<double>::iterator iter = actual_segment ->begin() ; iter != actual_segment->end() ; ++ iter )
			{
				double actual_value = ( * iter );
				double approximation = Polynomial ( count_i , coefficients );
				( * iter ) = actual_value - approximation;
#ifdef DEBUG
				std::cout << * actual_segment_iter << " " << total_counter << " " << count_i << " " << actual_value << " " << approximation << " " << actual_value - approximation << std::endl;
#endif
				if ( output_analysis_file.is_open ( ) )
				{
					output_analysis_file << * actual_segment_iter << " " << total_counter << " " << count_i << " " << actual_value << " " << approximation << " " << actual_value - approximation << std::endl;
				}
				count_i ++;
				total_counter ++;
			}
		}

// allocation space for fluctuation fuction data and setup of initial values
		std::vector<double> fluctuation_function ( q_th_order_fluctuation.size() );
		for ( unsigned int j = 0 ; j < q_th_order_fluctuation . size ( ) ; j ++ )
		{
			fluctuation_function [ j ] = 0;
		}

		for ( std::list< std::vector<double> * >::iterator segment_iter = segments.begin() ; segment_iter != segments.end() ; segment_++ iter )
		{
			std::vector<double> * actual_segment = ( * segment_iter );

			unsigned int count_i = 0;
			double segment_variance = 0;
			for ( std::vector<double>::iterator iter = actual_segment ->begin() ; iter != actual_segment->end() ; ++ iter )
			{
				double actual_value = ( * iter );

				segment_variance += ( actual_value * actual_value );
			}
			segment_variance /= * actual_segment_iter;

			for ( unsigned int i = 0 ; i < q_th_order_fluctuation . size ( ) ; i ++ )
			{
				if ( q_th_order_fluctuation [ i ] != 0.0 )
				{
					fluctuation_function [ i ] += pow ( segment_variance , q_th_order_fluctuation [ i ] / 2.0 );
				}
				else
				{
					fluctuation_function [ i ] += log ( segment_variance );
				}
			}
		}

		for ( unsigned int i = 0 ; i < q_th_order_fluctuation . size ( ) ; i ++ )
		{
			if ( q_th_order_fluctuation [ i ] != 0.0 )
			{
				fluctuation_function [ i ] /= segments.size();
				double data = pow ( fluctuation_function [ i ] , 1.0 / q_th_order_fluctuation [ i ] );

				DFAvariance [ count_segments ] .push_back ( data );
			}
			else
			{
				fluctuation_function [ i ] /= 2 * segments.size();
				double data = exp ( fluctuation_function [ i ] );

				DFAvariance [ count_segments ] .push_back ( data );
			}
		}

// clean up structures prepared
		for ( std::list< std::vector<double> * >::iterator segment_iter = segments.begin() ; segment_iter != segments.end() ; ++ segment_iter )
		{
			std::vector<double> * actual_list = ( * segment_iter );
			delete actual_list;
		}
		segments.clear();
	}
#endif
	if ( output_analysis_file.is_open() )
	{
		output_analysis_file.close ();
	}
}


void DetrendedFluctuationAnalysisExponentiallyScaledDomain ( std::list<double> & data , std::vector<unsigned int> & segment_sizes , std::map<std::tuple<double, unsigned int, unsigned int> , double > & DFA , std::vector<double> & q_th_order_fluctuation , unsigned int start_index , unsigned int end_index , unsigned segments , const std::vector<long> & order_of_method , std::string output_analysis , bool moving_segment_adjustment , unsigned int moving_segment_shift , bool sample_variance , bool double_integration
#ifdef THREADING
, unsigned int number_threads
#ifdef MICROTHREADING
, unsigned int microthreads
#endif
#endif
, void ( * fluctuation_analyzer ) ( std::vector<double> &  , unsigned int , std::list<double> &  , double , math_functions::real_matrix * ) , double ( * polynomial ) ( double  , std::list<double> & , double  ) )
{
// loop over segment sizes that are going to be iterated
	unsigned int max_segment_size = end_index;
	unsigned int segment_size = start_index;
	unsigned int counter = 0;

	double exponent = exp ( log ( double(end_index) / double(start_index) ) / double(segments) );

	while ( segment_size <= max_segment_size )
	{
		segment_sizes.push_back ( segment_size );

		counter ++;
		unsigned int new_segment_size = start_index * pow ( exponent , counter );
		if ( segment_size >= new_segment_size )
		{
			segment_size += 1;
		}
		else
		{
			segment_size = new_segment_size;
		}
	}

	DetrendedFluctuationAnalysis ( data , segment_sizes , DFA , q_th_order_fluctuation , order_of_method , output_analysis , moving_segment_adjustment , moving_segment_shift , sample_variance , double_integration
#ifdef THREADING
, number_threads
#ifdef MICROTHREADING
, microthreads
#endif
#endif
, fluctuation_analyzer , polynomial );

}


// void DetrendedFluctuationAnalysisExponentiallyScaledDomainSampleVarianceTest ( std::list<double> & data , std::list<unsigned int> & segment_sizes , std::vector<std::list<double> > & DFAvariance , std::vector<double> & q_th_order_fluctuation , unsigned int start_index , unsigned int end_index , unsigned segments , unsigned int order_of_method , std::string output_analysis , bool moving_segment_adjustment , unsigned int moving_segment_shift , bool sample_variance , bool double_integration
// #ifdef THREADING
// , unsigned int number_threads
// #ifdef MICROTHREADING
// , unsigned int microthreads
// #endif
// #endif
// , void ( * fluctuation_analyzer ) ( std::vector<double> &  , unsigned int , std::list<double> &  , double , math_functions::real_matrix * ) , double ( * polynomial ) ( double  , std::list<double> & , double  ) )
// {
// #ifdef DEBUG
// 	math_model::debug::OpenDebugStream ();
// #endif
// 
// 	if ( data.size() <= end_index )
// 	{
// 		std::cerr << "DetrendedFluctuationAnalysisExponentiallyScaledDomain: End index of DFA to large. Higher values are irrlelevant due to bad statistics" << std::endl;
// 		return;
// 	}
// 
// 	if ( order_of_method + 1 > start_index )
// 	{
// 		std::cerr << "DetrendedFluctuationAnalysisExponentiallyScaledDomain: Possibly wrong argument. Initial index should be at least order of the method + 1." << std::endl;
// 	}
// 
// // open output file
// 	std::ofstream output_analysis_file;
// 	if ( ! output_analysis.empty() )
// 	{
// 		output_analysis_file.open ( output_analysis.c_str() );
// 	}
// 
// 	std::list<double> integrated_dataset ( data );
// 	IntegrationOfDataset ( integrated_dataset );
// 
// 	if ( double_integration == true )
// 	{
// 		IntegrationOfDataset ( integrated_dataset );
// 	}
// 
// 	double exponent = exp ( log ( double(end_index) / double(start_index) ) / double(segments) );
// 
// #ifdef THREADING
// 	unsigned int max_segment_size = end_index;
// 	unsigned int segment_size = start_index;
// 	unsigned int counter = 0;
// #ifdef QTHREADS
// 	std::list<SegmentDFAWorker *> DFA_threads;
// #else
// 	boost::thread_group * DFA_threads = new boost::thread_group;
// #endif
// 
// // preparing space for output data with different q
// 	for ( unsigned int i = 0 ; i < q_th_order_fluctuation.size() ; i ++ )
// 	{
// 		std::list<double> * new_list = new std::list<double>;
// 		DFAvariance.push_back ( * new_list );
// 	}
// 
// // loop over segment sizes that are going to be iterated
// 	while ( segment_size <= max_segment_size )
// 	{
// // praparation of storage for the data and preparation of structure that will used for access
// 		std::vector<double *> DFAvariance_results;
// 		for ( unsigned int i = 0 ; i < q_th_order_fluctuation.size() ; i ++ )
// 		{
// 			DFAvariance [ i ] . push_back ( 0.0 );
// 			double * DFAvariance_result = & ( DFAvariance [ i ] . back ( ) );
// 			DFAvariance_results.push_back ( DFAvariance_result );
// 		}
// 
// 		segment_sizes . push_back ( segment_size );
// 		SegmentDFAWorker * new_thread = new SegmentDFAWorker ( integrated_dataset , q_th_order_fluctuation , DFAvariance_results , output_analysis_file , order_of_method , segment_size , moving_segment_adjustment , moving_segment_shift , sample_variance
// #ifdef MICROTHREADING
//  , microthreads
// #endif
// , fluctuation_analyzer , polynomial
// 		);
// #ifdef QTHREADS
// 		DFA_threads .push_back ( new_thread );
// 		new_thread->start();
// #else
// 		DFA_threads -> add_thread ( new boost::thread ( * new_thread ) );
// #endif
// 
// 		counter ++;
// 		unsigned int new_segment_size = start_index * pow ( exponent , counter );
// 		if ( segment_size >= new_segment_size )
// 		{
// 			segment_size += 1;
// 		}
// 		else
// 		{
// 			segment_size = new_segment_size;
// 		}
// 
// #ifdef QTHREADS
// 		if ( DFA_threads.size ( ) > number_threads )
// 		{
// #else
// 		if ( DFA_threads->size ( ) > number_threads )
// 		{
// #endif
// #ifdef QTHREADS
// 			SegmentDFAWorker * actual_thread = DFA_threads.back();
// 			actual_thread->wait();
// 			DFA_threads . pop_back ( );
// 			delete actual_thread;
// #else
// // 			boost::posix_time::seconds workTime ( 10 );
// // 			boost::this_thread::sleep(workTime);
// 
// 			DFA_threads -> join_all();
// 			delete DFA_threads;
// 			DFA_threads = new boost::thread_group;
// #endif
// 		}
// 	}
// 
// #ifdef QTHREADS
// 	while ( DFA_threads.size ( ) )
// 	{
// #else
// 	{
// #endif
// 
// #ifdef QTHREADS
// 		SegmentDFAWorker * actual_thread = DFA_threads.back();
// 		actual_thread->wait();
// 		DFA_threads . pop_back ( );
// 		delete actual_thread;
// #else
// // 		boost::posix_time::seconds workTime ( 10 );
// // 		boost::this_thread::sleep(workTime);
// 
// 		DFA_threads -> join_all();
// 		delete DFA_threads;
// #endif
// 	}
// #else
// 
// // adding list for results
// 	for ( unsigned int j = 0 ; j < q_th_order_fluctuation . size ( ) ; j ++ )
// 	{
// 		std::list<double> list_to_add;
// 		DFAvariance.push_back(list_to_add);
// 	}
// 
// // loop over segment size
// 	unsigned int max_segment_size = end_index;
// 	unsigned int segment_size = start_index - 1;
// 	for ( unsigned int i = 0 ; i <= segments ; i ++ )
// 	{
// 		unsigned int new_segment_size = start_index * pow ( exponent , i );
// 		if ( segment_size >= new_segment_size )
// 		{
// 			segment_size += 1;
// 		}
// 		else
// 		{
// 			segment_size = new_segment_size;
// 		}
// 		segment_sizes.push_back( segment_size );
// 
// 		std::list< std::vector<double> * > segments;
// 		std::vector<double> * new_segment = new std::vector<double> ( segment_size );
// 
// 		unsigned int counter = 0;
// 		unsigned int total_counter = 0;
// 
// // preparation of cached matrices for regression procedure
// 		math_functions::real_matrix transformation_y_data_cached;
// 		PrepareMatricesForRegression ( segment_size  , order_of_method , transformation_y_data_cached );
// 
// 		// preparation of segments starting from first data to last
// 		for ( std::list<double>::iterator iter_data = data.begin() ; iter_data != data.end() ; iter_data ++ )
// 		{
// 			double actual_data = ( * iter_data );
// 			total_counter++;
// 
// 			if ( counter < segment_size )
// 			{
// 				( * new_segment ) [ counter ] = actual_data;
// 				counter ++;
// 			}
// 			else
// 			{
// 				counter = 0;
// 				segments . push_back( new_segment );
// 				if ( data.size() - total_counter >= segment_size )
// 				{
// 					new_segment = new std::vector<double> ( segment_size );
// 					( * new_segment ) [ counter ] = actual_data;
// 					counter ++;
// 				}
// 				else
// 				{
// 					break;
// 				}
// 			}
// 		}
// 
// 		// preparation of segments starting from last data to first
// 		counter = 0;
// 		total_counter = 0;
// 		new_segment = new std::vector<double> ( segment_size );
// 		for ( std::list<double>::reverse_iterator iter_data = data.rbegin() ; iter_data != data.rend() ; iter_data ++ )
// 		{
// 			double actual_data = ( * iter_data );
// 			total_counter++;
// 
// 			if ( counter < segment_size )
// 			{
// 				( * new_segment ) [ segment_size - 1 - counter ] = actual_data;
// 				counter ++;
// 			}
// 			else
// 			{
// 				counter = 0;
// 				segments . push_back( new_segment );
// 				if ( data.size() - total_counter >= segment_size )
// 				{
// 					new_segment = new std::vector <double> ( segment_size );
// 					( * new_segment ) [ segment_size - 1 - counter ] = actual_data;
// 					counter ++;
// 				}
// 				else
// 				{
// 					break;
// 				}
// 			}
// 		}
// 
// 		total_counter = 0;
// 		for ( std::list< std::vector<double> * >::iterator segment_iter = segments.begin() ; segment_iter != segments.end() ; segment_++ iter )
// 		{
// 			std::vector<double> * actual_segment = ( * segment_iter );
// 			std::list<double> coefficients;
// 			fluctuation_analyzer (  * actual_segment, order_of_method , coefficients , 1.0 , & transformation_y_data_cached );
// //			PolynomialRegression<std::vector<double>,std::vector<double>::iterator> ( * actual_segment, order_of_method , coefficients , & transformation_y_data_cached );
// 
// 			unsigned int count_i = 0;
// 			for ( std::vector<double>::iterator iter = actual_segment ->begin() ; iter != actual_segment->end() ; ++ iter )
// 			{
// 				double actual_value = ( * iter );
// // 				double approximation = Polynomial ( count_i , coefficients );
// 				double approximation = polynomial ( count_i , coefficients , 1.0 );
// 				( * iter ) = actual_value - approximation;
// #ifdef DEBUG
// 				std::cout << segment_size << " " << total_counter << " " << count_i << " " << actual_value << " " << approximation << " " << actual_value - approximation << std::endl;
// #endif
// 				if ( output_analysis_file.is_open ( ) )
// 				{
// 					output_analysis_file << segment_size << " " << total_counter << " " << count_i << " " << actual_value << " " << approximation << " " << actual_value - approximation << std::endl;
// 				}
// 				count_i ++;
// 				total_counter ++;
// 			}
// 		}
// 
// // allocation space for fluctuation fuction data and setup of initial values
// 		std::vector<double> fluctuation_function ( q_th_order_fluctuation.size() );
// 		for ( unsigned int j = 0 ; j < q_th_order_fluctuation . size ( ) ; j ++ )
// 		{
// 			fluctuation_function [ j ] = 0;
// 		}
// 
// 		for ( std::list< std::vector<double> * >::iterator segment_iter = segments.begin() ; segment_iter != segments.end() ; segment_++ iter )
// 		{
// 			std::vector<double> * actual_segment = ( * segment_iter );
// 
// 			unsigned int count_i = 0;
// 			double segment_variance = 0;
// 			for ( std::vector<double>::iterator iter = actual_segment ->begin() ; iter != actual_segment->end() ; ++ iter )
// 			{
// 				double actual_value = ( * iter );
// 
// 				segment_variance += ( actual_value * actual_value );
// 			}
// 			segment_variance /= ( segment_size - 1 );
// 
// 			for ( unsigned int j = 0 ; j < q_th_order_fluctuation . size ( ) ; j ++ )
// 			{
// 				if ( q_th_order_fluctuation [ j ] != 0.0 )
// 				{
// 					fluctuation_function [ j ] += pow ( sqrt ( segment_variance ) , q_th_order_fluctuation [ j ] );
// 				}
// 				else
// 				{
// 					fluctuation_function [ j ] = log ( segment_variance );
// 				}
// 			}
// 		}
// 
// 		for ( unsigned int j = 0 ; j < q_th_order_fluctuation . size ( ) ; j ++ )
// 		{
// 			if ( q_th_order_fluctuation [ j ] != 0.0 )
// 			{
// 				fluctuation_function [ j ] /= segments.size();
// 				double data = pow ( fluctuation_function [ j ] , 1.0 / q_th_order_fluctuation [ j ] );
// 
// 				DFAvariance [ j ] . push_back ( data );
// 			}
// 			else
// 			{
// 				fluctuation_function [ j ] /= 2 * segments.size();
// 				double data = exp ( fluctuation_function [ j ] );
// 
// 				DFAvariance [ j ] . push_back ( data );
// 			}
// 		}
// 
// // clean up structures prepared
// 		for ( std::list< std::vector<double> * >::iterator segment_iter = segments.begin() ; segment_iter != segments.end() ; segment_++ iter )
// 		{
// 			std::vector<double> * actual_list = ( * segment_iter );
// 			delete actual_list;
// 		}
// 		segments.clear();
// 	}
// #endif
// 	if ( output_analysis_file.is_open() )
// 	{
// 		output_analysis_file.close ();
// 	}
// }


void DetrendedHarmonicFluctuationAnalysis ( std::list<double> & data , std::list<double> & DFAvariance , double q_th_order_fluctuation , unsigned int order_of_method , unsigned segment_step , unsigned int min_number_of_segments , std::string output_analysis )
{
	if ( data.size() <= min_number_of_segments * ( order_of_method + 1 ) )
	{
		std::cerr << "DetrendedFluctuationAnalysis: Wrong number of arguments" << std::endl;
		return;
	}

	std::ofstream output_analysis_file;
	if ( ! output_analysis.empty() )
	{
		output_analysis_file.open ( output_analysis.c_str() );
	}

	unsigned int max_segment_size = data.size() / min_number_of_segments;
	for ( unsigned int segment_size = order_of_method + 2 ; segment_size <= max_segment_size ; segment_size += segment_step )
	{
		std::list< std::list<double> * > segments;
		segments . push_back( new std::list<double> );
		unsigned int counter = 0;
		unsigned int total_counter = 0;

		// preparation of segments starting from first data to last
		for ( std::list<double>::iterator iter_data = data.begin() ; iter_data != data.end() ; iter_data ++ )
		{
			double actual_data = ( * iter_data );
			total_counter++;

			if ( counter < segment_size )
			{
				segments.back()->push_back( actual_data );
				counter ++;
			}
			else
			{
				counter = 0;
				if ( data.size() - total_counter >= segment_size )
				{
					segments . push_back( new std::list<double> );
					segments.back()->push_back( actual_data );
					counter ++;
				}
				else
				{
					break;
				}
			}
		}

		// preparation of segments starting from last data to first
		counter = 0;
		total_counter = 0;
		segments . push_back( new std::list<double> );
		for ( std::list<double>::reverse_iterator iter_data = data.rbegin() ; iter_data != data.rend() ; iter_data ++ )
		{
			double actual_data = ( * iter_data );
			total_counter++;

			if ( counter < segment_size )
			{
				segments.back()->push_front( actual_data );
				counter ++;
			}
			else
			{
				counter = 0;
				if ( data.size() - total_counter >= segment_size )
				{
					segments . push_back( new std::list<double> );
					segments.back()->push_front( actual_data );
					counter ++;
				}
				else
				{
					break;
				}
			}
		}

		total_counter = 0;
		for ( std::list< std::list<double> * >::iterator segment_iter = segments.begin() ; segment_iter != segments.end() ; ++ segment_iter )
		{
			std::list<double> * actual_segment = ( * segment_iter );
			std::list<double> coefficients;

			HarmonicRegression<std::list<double>,std::list<double>::iterator> ( * actual_segment, order_of_method , coefficients , segment_size );

			unsigned int count_i = 0;
			for ( std::list<double>::iterator iter = actual_segment ->begin() ; iter != actual_segment->end() ; ++ iter )
			{
				double actual_value = ( * iter );

				double approximation = HarmonicSerie ( count_i , coefficients , segment_size );

				( * iter ) = actual_value - approximation;
#ifdef DEBUG
				std::cout << segment_size << " " << total_counter << " " << count_i << " " << actual_value << " " << approximation << " " << actual_value - approximation << std::endl;
#endif
				if ( output_analysis_file.is_open ( ) )
				{
					output_analysis_file << segment_size << " " << total_counter << " " << count_i << " " << actual_value << " " << approximation << " " << actual_value - approximation << std::endl;
				}

				count_i ++;
				total_counter ++;
			}
		}

		double fluctuation_function = 0;
		for ( std::list< std::list<double> * >::iterator segment_iter = segments.begin() ; segment_iter != segments.end() ; ++ segment_iter )
		{
			std::list<double> * actual_segment = ( * segment_iter );

			double segment_variance = 0;
			for ( std::list<double>::iterator iter = actual_segment ->begin() ; iter != actual_segment->end() ; ++ iter )
			{
				double actual_value = ( * iter );

				segment_variance += ( actual_value * actual_value );
			}
			segment_variance /= segment_size;

			if ( q_th_order_fluctuation != 0.0 )
			{
				fluctuation_function += pow ( segment_variance , q_th_order_fluctuation / 2.0 );
			}
			else
			{
				fluctuation_function += log ( segment_variance );
			}
		}
		if ( q_th_order_fluctuation != 0.0 )
		{
			fluctuation_function /= segments.size();
			double data = pow ( fluctuation_function , 1.0 / q_th_order_fluctuation );

			DFAvariance.push_back( data );
		}
		else
		{
			fluctuation_function /= 2 * segments.size();
			double data = exp ( fluctuation_function );

			DFAvariance.push_back( data );
		}

// clean up structures prepared
		for ( std::list< std::list<double> * >::iterator segment_iter = segments.begin() ; segment_iter != segments.end() ; ++ segment_iter )
		{
			std::list<double> * actual_list = ( * segment_iter );
			delete actual_list;
		}
		segments.clear();
	}
	if ( output_analysis_file.is_open() )
	{
		output_analysis_file.close ();
	}
}

void Autocorrelation ( std::vector<double> data , std::list<double> & autocorrelation , unsigned segment_step , unsigned int  min_number_of_segments)
{
	if ( data.size() == 0 )
	{
		std::cerr << "Autocorrelation: no data to calculate" << std::endl;
	}

	double mean = 0;
	double sqr_x = 0;
	for ( unsigned int i = 0 ; i < data.size() ; i ++ )
	{
		mean += data [ i ];
		sqr_x += data [ i ] * data [ i ];
	}
	mean /= data.size();
	double variance = sqr_x / ( data.size() ) - mean * mean * (  data.size() ) / (  data.size() );

	unsigned int limit_step_size = data.size() / min_number_of_segments;
	for ( unsigned int step_size = 0 ; step_size <= limit_step_size ; step_size ++ )
	{
		double actual_autocorrelation = 0;
		for ( unsigned int final_index = step_size ; final_index < data.size() ; final_index ++ )
		{
			unsigned int initial_index = final_index - step_size;

			actual_autocorrelation += ( ( data [ final_index ] - mean ) * ( data [ initial_index ] - mean ) );
		}
		actual_autocorrelation /= variance * ( data.size() - step_size );

		autocorrelation.push_back(actual_autocorrelation);
	}
}

void AutocorrelationCyclicData ( std::vector<double> data , std::list<double> & autocorrelation , unsigned segment_step , unsigned int  min_number_of_segments)
{
	if ( data.size() == 0 )
	{
		std::cerr << "Autocorrelation: no data to calculate" << std::endl;
	}

	double mean = 0;
	double sqr_x = 0;
	for ( unsigned int i = 0 ; i < data.size() ; i ++ )
	{
		mean += data [ i ];
		sqr_x += data [ i ] * data [ i ];
	}
	mean /= data.size();
	double variance = sqr_x / ( data.size() ) - mean * mean * (  data.size() ) / (  data.size() );

	unsigned int limit_step_size = data.size() / min_number_of_segments;
	for ( unsigned int step_size = 0 ; step_size <= limit_step_size ; step_size ++ )
	{
		double actual_autocorrelation = 0;
		for ( unsigned int final_index = 0 ; final_index < data.size() ; final_index ++ )
		{
			int initial_index = final_index - step_size;
			if ( initial_index < 0 )
			{
				initial_index += data.size();
			}

			actual_autocorrelation += ( ( data [ final_index ] - mean ) * ( data [ initial_index ] - mean ) );
		}
		actual_autocorrelation /= variance * ( data.size() );

		autocorrelation.push_back(actual_autocorrelation);
	}
}

void StructureFunction ( std::vector<double> data , std::list<double> & timelagfunction , double q_th_order , unsigned segment_step , unsigned int min_number_of_segments )
{
	if ( data.size() == 0 )
	{
		std::cerr << "Autocorrelation: no data to calculate" << std::endl;
	}

	unsigned int limit_step_size = data.size() / min_number_of_segments;
	for ( unsigned int step_size = 0 ; step_size <= limit_step_size ; step_size ++ )
	{
		double actual_autocorrelation = 0;
		for ( unsigned int final_index = step_size ; final_index < data.size() ; final_index ++ )
		{
			unsigned int initial_index = final_index - step_size;

			double data_1 = data [ initial_index ];
			double data_2 = data [ final_index ];
			actual_autocorrelation += pow ( fabs ( data_2 - data_1 ) , q_th_order );
		}
		actual_autocorrelation /= ( data.size() - step_size );

		timelagfunction.push_back( pow ( actual_autocorrelation , 1.0 / q_th_order ) );
	}
}

void StructureFunctionCyclicData ( std::vector<double> data , std::list<double> & timelagfunction , double q_th_order , unsigned segment_step , unsigned int min_number_of_segments )
{
	if ( data.size() == 0 )
	{
		std::cerr << "Autocorrelation: no data to calculate" << std::endl;
	}

	unsigned int limit_step_size = data.size() / min_number_of_segments;
	for ( unsigned int step_size = 0 ; step_size <= limit_step_size ; step_size ++ )
	{
		double actual_autocorrelation = 0;
		for ( unsigned int final_index = step_size ; final_index < data.size() ; final_index ++ )
		{
			int initial_index = final_index - step_size;
			if ( initial_index < 0 )
			{
				initial_index += data.size();
			}

			double data_1 = data [ initial_index ];
			double data_2 = data [ final_index ];
			actual_autocorrelation += pow ( fabs ( data_2 - data_1 ) , q_th_order );
		}
		actual_autocorrelation /= ( data.size() - step_size );

		timelagfunction.push_back( pow ( actual_autocorrelation , 1.0 / q_th_order ) );
	}
}

/*
struct extremum
{
	int _position;
	double _value;

	extremum ( int position , double value )
	: _position( position ) , _value ( value )
	{
	}
};
*/

bool IsLocalMaximum ( double y1 , double y2 , double y3 )
{
	if ( ( y2 > y1 ) && ( y2 > y3 ) )
	{
		return true;
	}
	else
	{
		return false;
	}
}


bool IsLocalMinimum ( double y1 , double y2 , double y3 )
{
	if ( ( y2 < y1 ) && ( y2 < y3 ) )
	{
		return true;
	}
	else
	{
		return false;
	}
}

/*
bool SiftingOfDataset ( std::vector<double> & data , std::vector<double>::iterator data_begin , std::vector<double>::iterator data_end , std::vector<double> & sifted_dataset , std::vector<double> & fluctuation_part_of_dataset , bool extended_algorithm , double level_of_difference )
{
// searching for maxima and minima of the dataset
	std::list<extremum *> maxima;
	std::list<extremum *> minima;

// possible enhancement regarding maximum and minimum from previous segment
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
	}
	else
	{
		std::vector<double>::iterator before_iter = data_begin;
		before_iter --;
		std::vector<double>::iterator iter = data_begin;
		std::vector<double>::iterator after_iter = data_begin;
		++ after_iter;
		int datapoint = 0;

		for (  ; before_iter != data.begin() ; )
		{
			if ( IsLocalMinimum ( ( * before_iter ) , ( * iter ) , ( * after_iter ) ) )
			{
				minima . push_back ( new extremum ( datapoint , ( * iter ) ) );
#ifdef DEBUG
				std::cout << "Minimum: " << datapoint  << " " << ( * iter ) << std::endl;
#endif
				break;
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

		before_iter = data_begin;
		-- before_iter;
		iter = data_begin;
		after_iter = data_begin;
		++ after_iter;
		datapoint = 0;

		for ( ; before_iter != data.begin() ; )
		{
			if ( IsLocalMaximum ( ( * before_iter ) , ( * iter ) , ( * after_iter ) ) )
			{
				maxima . push_back ( new extremum ( datapoint , ( * iter ) ) );
#ifdef DEBUG
				std::cout << "Maximum: " << datapoint  << " " << ( * iter ) << std::endl;
#endif
				break;
			}

			after_iter = iter;
			iter = before_iter;
			before_iter --;
			datapoint --;
		}
		if ( before_iter == data.begin() )
		{
			maxima . push_back ( new extremum ( 0 , ( * data_begin ) ) );
#ifdef DEBUG
			std::cout << "Maximum: " << 0 << " " << ( * data_begin ) << std::endl;
#endif
		}
	}

	unsigned int datapoint = 1;
	std::vector<double>::iterator end_iter = data_end;
	std::vector<double>::iterator iter = data_begin;
	++ iter;
	std::vector<double>::iterator before_iter = data_begin;
	std::vector<double>::iterator after_iter = iter;
	++ after_iter;

// searching for local extrema in the interval
	for ( ; after_iter != end_iter ; )
	{
		double actual_datapoint = ( * iter );

		if ( IsLocalMinimum ( ( * before_iter ) , ( * iter ) , ( * after_iter ) ) )
		{
			minima . push_back ( new extremum ( datapoint , ( * iter ) ) );
#ifdef DEBUG
			std::cout << "Minimum: " << datapoint  << " " << ( * iter ) << std::endl;
#endif
		}
		else
		{
			if ( IsLocalMaximum ( ( * before_iter ) , ( * iter ) , ( * after_iter ) ) )
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

	if ( ( extended_algorithm == false ) || ( data_end == data . end ( ) ) )
	{
		std::vector<double>::iterator temp_iter = data_end;
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

		for ( ; after_iter != data.end() ; )
		{
			if ( IsLocalMinimum ( ( * before_iter ) , ( * iter ) , ( * after_iter ) ) )
			{
				minima . push_back ( new extremum ( temp_datapoint , ( * iter ) ) );
#ifdef DEBUG
				std::cout << "Additional minimum: " << temp_datapoint  << " " << ( * iter ) << std::endl;
#endif
				break;
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

		temp_datapoint = datapoint;
		iter = data_end;
		iter --;
		before_iter = data_end;
		before_iter --;
		before_iter --;
		after_iter = data_end;

		for ( ; after_iter != data.end() ; )
		{
			if ( IsLocalMaximum ( ( * before_iter ) , ( * iter ) , ( * after_iter ) ) )
			{
				maxima . push_back ( new extremum ( temp_datapoint , ( * iter ) ) );
#ifdef DEBUG
				std::cout << "Additional maximum: " << temp_datapoint  << " " << ( * iter ) << std::endl;
#endif
				break;
			}
			before_iter = iter;
			iter = after_iter;
			++ after_iter;
			++ temp_datapoint;
		}
		if ( after_iter == data.end() )
		{
			maxima . push_back ( new extremum ( temp_datapoint , ( * iter ) ) );
#ifdef DEBUG
			std::cout << "Additional maximum: " << temp_datapoint  << " " << ( * iter ) << std::endl;
#endif
		}
	}

	int number_of_extrema = 0;

// interpolation of maxima and minima
// initialization of the library
	gsl_interp_accel * acc_maxima = gsl_interp_accel_alloc ();
	gsl_interp_accel * acc_minima = gsl_interp_accel_alloc ();
	const gsl_interp_type * type_of_interpolation = gsl_interp_cspline; 
	gsl_spline * spline_maxima = gsl_spline_alloc ( type_of_interpolation , maxima.size() );
	gsl_spline * spline_minima = gsl_spline_alloc ( type_of_interpolation , minima.size() );

// preparation of data
	double * positions_maxima = ( double * ) malloc ( maxima.size() * sizeof ( double ) );
	double * dataset_maxima = ( double * ) malloc ( maxima.size() * sizeof ( double ) );
	std::list<extremum *>::iterator iter2 = maxima.begin();
	for ( unsigned int i = 0 ; i < maxima.size() ; i ++ )
	{
		positions_maxima [ i ] = ( * iter2 ) ->_position;
		dataset_maxima [ i ] = ( * iter2 ) ->_value;
		++ iter2;

		if ( ( positions_maxima [ i ] > 0 ) && ( positions_maxima [ i ] < datapoint ) )
		{
			number_of_extrema ++;
		}
	}

	double * positions_minima = ( double * ) malloc ( minima.size() * sizeof ( double ) );
	double * dataset_minima = ( double * ) malloc ( minima.size() * sizeof ( double ) );

	std::list<extremum *>::iterator iter3 = minima . begin ( );
	for ( unsigned int i = 0 ; i < minima . size ( ) ; i ++ )
	{
		positions_minima [ i ] = ( * iter3 ) ->_position;
		dataset_minima [ i ] = ( * iter3 ) ->_value;
		++ iter3;

		if ( ( positions_minima [ i ] > 0 ) && ( positions_minima [ i ] < datapoint ) )
		{
			number_of_extrema ++;
		}
	}

// interpolation
	gsl_spline_init ( spline_maxima , positions_maxima , dataset_maxima , maxima . size ( ) );
	gsl_spline_init ( spline_minima , positions_minima , dataset_minima , minima . size ( ) );

// filtration
#ifdef DEBUG
	for ( unsigned int xi = 0 ; xi <= datapoint ; xi += 1 )
	{
		std::cout << xi << " max: " << gsl_spline_eval (spline_maxima, xi, acc_maxima) << " min: " << gsl_spline_eval (spline_minima, xi, acc_minima) << " mean: " << ( gsl_spline_eval (spline_maxima, xi, acc_maxima) + gsl_spline_eval (spline_minima, xi, acc_minima) ) / 2.0  << std::endl;
	}
#endif

// test
	int crossings = 0;
	double maxima_mean = 0;
	double minima_mean = 0;
	std::vector<double>::iterator iter4 = data_begin;
	iter4 ++;
	bool above = ( ( gsl_spline_eval (spline_maxima, 1, acc_maxima ) + gsl_spline_eval (spline_minima, 1, acc_minima ) ) / 2.0 > ( * iter4 ) ) ? false : true;
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
	if ( ( number_of_extrema + 1 >= crossings ) && ( crossings >= number_of_extrema - 1 ) && ( fabs ( maxima_mean + minima_mean ) <= level_of_difference ) )
	{
#ifdef DEBUG
		std::cout << "Sifting accepted" << std::endl;
#endif
		std::vector<double>::iterator dataset_iter = data_begin;
		std::vector<double>::iterator sift_end = sifted_dataset.end();
		std::vector<double>::iterator sift_iter = sifted_dataset.begin();
		std::vector<double>::iterator fluctuation_iter = fluctuation_part_of_dataset.begin();
		std::vector<double>::iterator fluctuation_end = fluctuation_part_of_dataset.end();
		int position = 0;
		for ( ; dataset_iter != data_end ; )
		{
			double approximation = ( gsl_spline_eval (spline_maxima, position, acc_maxima) + gsl_spline_eval (spline_minima, position, acc_minima) ) / 2.0;
			double & sift_data = ( * sift_iter );
			double & fluc_data = ( * fluctuation_iter );
			double & original_data = ( * dataset_iter );
			if ( sift_iter != sift_end)
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
		return_value = true;
	}
	else
	{
#ifdef DEBUG
		std::cout << "Sifting is not accepted" << std::endl;
#endif
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
*/

void EmpiricalModeDecomposionDetrendedFluctuationAnalysisExponentiallyScaledDomain ( std::list<double> & data , std::list<unsigned int> & segment_sizes , std::vector<std::list<double> > & DFAvariance , std::vector<double> & q_th_order_fluctuation , unsigned int start_index , unsigned int end_index , unsigned segments , std::string output_analysis , bool moving_segment_adjustment , unsigned int moving_segment_shift , bool sample_variance , bool double_integration , int minimal_sifting
#ifdef THREADING
, unsigned int number_threads
#ifdef MICROTHREADING
, unsigned int microthreads
#endif
#endif
, bool ( * sifting_method )  ( std::list<double> & , std::list<double>::iterator , std::list<double>::iterator , std::list<double> & , std::list<double> & , std::list<double>::iterator & , std::list<double>::iterator & , bool , bool , bool , double ) )
{
#ifdef DEBUG
	math_model::debug::OpenDebugStream ();
#endif

	if ( data.size() <= end_index )
	{
		std::cerr << "EmpiricalModeDecomposionDetrendedFluctuationAnalysisExponentiallyScaledDomain: End index of DFA to large. Higher values are irrlelevant due to bad statistics" << std::endl;
		return;
	}

// open output file
	std::ofstream output_analysis_file;
	if ( ! output_analysis.empty() )
	{
		output_analysis_file.open ( output_analysis.c_str() );
	}

	std::list<double> integrated_dataset ( data );
	IntegrationOfDataset ( integrated_dataset );

	if ( double_integration == true )
	{
		IntegrationOfDataset ( integrated_dataset );
	}

	double exponent = exp ( log ( double(end_index) / double(start_index) ) / double(segments) );

#ifdef THREADING
	unsigned int max_segment_size = end_index;
	unsigned int segment_size = start_index;
	unsigned int counter = 0;
#ifdef QTHREADS
	std::list<SegmentEMD_DFAWorker *> DFA_threads;
#else
	boost::thread_group * DFA_threads = new boost::thread_group;
#endif

// preparing space for output data with different q
	for ( unsigned int i = 0 ; i < q_th_order_fluctuation.size() ; i ++ )
	{
		std::list<double> * new_list = new std::list<double>;
		DFAvariance.push_back ( * new_list );
	}

// loop over segment sizes that are going to be iterated
	while ( segment_size <= max_segment_size )
	{
// praparation of storage for the data and preparation of structure that will used for access
		std::vector<double *> DFAvariance_results;
		for ( unsigned int i = 0 ; i < q_th_order_fluctuation.size() ; i ++ )
		{
			DFAvariance [ i ] . push_back ( 0.0 );
			double * DFAvariance_result = & ( DFAvariance [ i ] . back ( ) );
			DFAvariance_results.push_back ( DFAvariance_result );
		}

		segment_sizes . push_back ( segment_size );
		SegmentEMD_DFAWorker * new_thread = new SegmentEMD_DFAWorker ( integrated_dataset , q_th_order_fluctuation , DFAvariance_results , output_analysis_file ,  segment_size , moving_segment_adjustment , moving_segment_shift , sample_variance , minimal_sifting , sifting_method 
#ifdef MICROTHREADING
 , microthreads
#endif
	, 10
//  , floor( log2 ( segment_size ) ) + 1
		);
#ifdef QTHREADS
		DFA_threads .push_back ( new_thread );
		new_thread->start();
#else
		DFA_threads -> add_thread ( new boost::thread ( * new_thread ) );
#endif

		counter ++;
		unsigned int new_segment_size = start_index * pow ( exponent , counter );
		if ( segment_size >= new_segment_size )
		{
			segment_size += 1;
		}
		else
		{
			segment_size = new_segment_size;
		}

#ifdef QTHREADS
		if ( DFA_threads.size ( ) > number_threads )
		{
#else
		if ( DFA_threads->size ( ) > number_threads )
		{
#endif
#ifdef QTHREADS
			SegmentDFAWorker * actual_thread = DFA_threads.back();
			actual_thread->wait();
			DFA_threads . pop_back ( );
			delete actual_thread;
#else
// 			boost::posix_time::seconds workTime ( 10 );
// 			boost::this_thread::sleep(workTime);

			DFA_threads -> join_all();
			delete DFA_threads;
			DFA_threads = new boost::thread_group;
#endif
		}
	}

#ifdef QTHREADS
	while ( DFA_threads.size ( ) )
	{
#else
	{
#endif

#ifdef QTHREADS
		SegmentDFAWorker * actual_thread = DFA_threads.back();
		actual_thread->wait();
		DFA_threads . pop_back ( );
		delete actual_thread;
#else
// 		boost::posix_time::seconds workTime ( 10 );
// 		boost::this_thread::sleep(workTime);

		DFA_threads -> join_all();
		delete DFA_threads;
#endif
	}
#endif
}

void EmpiricalModeDecomposionDetrendedFluctuationAnalysis ( std::list<double> & data , std::list<unsigned int> & segment_sizes , std::vector<std::list<double> > & DFAvariance , std::vector<double> & q_th_order_fluctuation , unsigned int start_index , unsigned int end_index , unsigned segment_step , std::string output_analysis , bool moving_segment_adjustment , unsigned int moving_segment_shift , bool sample_variance , bool double_integration , int minimal_sifting
#ifdef THREADING
, unsigned int number_threads
#ifdef MICROTHREADING
, unsigned int microthreads
#endif
#endif
, bool ( * sifting_method )  ( std::list<double> & , std::list<double>::iterator , std::list<double>::iterator , std::list<double> & , std::list<double> & , std::list<double>::iterator & , std::list<double>::iterator & , bool , bool , bool , double ) )
{
#ifdef DEBUG
	math_model::debug::OpenDebugStream ();
#endif

	if ( data.size() <= end_index * ( 4 ) )
	{
		std::cerr << "EmpiricalModeDecomposionDetrendedFluctuationAnalysis: End index of DFA to large. Higher values are irrlelevant due to bad statistics" << std::endl;
		return;
	}

// open output file
	std::ofstream output_analysis_file;
	if ( ! output_analysis.empty() )
	{
		output_analysis_file.open ( output_analysis.c_str() );
	}

	std::list<double> integrated_dataset ( data );
	IntegrationOfDataset ( integrated_dataset );

	if ( double_integration == true )
	{
		IntegrationOfDataset ( integrated_dataset );
	}

#ifdef THREADING
	unsigned int max_segment_size = end_index;
	unsigned int segment_size = start_index;
#ifdef QTHREADS
	std::list<SegmentEMD_DFAWorker *> DFA_threads;
#else
	boost::thread_group * DFA_threads = new boost::thread_group;
#endif

// preparing space for output data with different q
	for ( unsigned int i = 0 ; i < q_th_order_fluctuation.size() ; i ++ )
	{
		std::list<double> new_list;
		DFAvariance.push_back ( new_list );
	}

// loop over segment sizes that are going to be iterated
	while ( segment_size <= max_segment_size )
	{
		std::vector<double *> DFAvariance_results;
		for ( unsigned int i = 0 ; i < q_th_order_fluctuation.size() ; i ++ )
		{
			DFAvariance [ i ] . push_back ( 0.0 );
			double * DFAvariance_result = & ( DFAvariance [ i ] . back ( ) );
			DFAvariance_results.push_back ( DFAvariance_result );
		}

		segment_sizes.push_back ( segment_size );
		SegmentEMD_DFAWorker * new_thread = new SegmentEMD_DFAWorker ( data , q_th_order_fluctuation , DFAvariance_results , output_analysis_file , segment_size , moving_segment_adjustment , moving_segment_shift , sample_variance , minimal_sifting , sifting_method
#ifdef MICROTHREADING
 , microthreads
#endif
		);
#ifdef QTHREADS
		DFA_threads .push_back ( new_thread );
		new_thread->start();
#else
		DFA_threads -> add_thread ( new boost::thread ( * new_thread ) );
#endif
		segment_size += segment_step;

#ifdef QTHREADS
		if ( DFA_threads.size ( ) > number_threads )
		{
#else
		if ( DFA_threads->size ( ) > number_threads )
		{
#endif
#ifdef QTHREADS
			SegmentDFAWorker * actual_thread = DFA_threads.back();
			actual_thread->wait();
			DFA_threads . pop_back ( );
			delete actual_thread;
#else
// 			boost::posix_time::seconds workTime ( 10 );
// 			boost::this_thread::sleep(workTime);

			DFA_threads -> join_all();
			delete DFA_threads;
			DFA_threads = new boost::thread_group;
#endif
		}
	}

#ifdef QTHREADS
	while ( DFA_threads.size ( ) )
	{
#else
	{
#endif

#ifdef QTHREADS
		SegmentDFAWorker * actual_thread = DFA_threads.back();
		actual_thread->wait();
		DFA_threads . pop_back ( );
		delete actual_thread;
#else
// 		boost::posix_time::seconds workTime ( 10 );
// 		boost::this_thread::sleep(workTime);

		DFA_threads -> join_all();
		delete DFA_threads;
#endif
	}
#endif
}

}
