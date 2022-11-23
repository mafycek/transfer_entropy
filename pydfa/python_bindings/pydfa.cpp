/***************************************************************************
 *   Copyright (C) 2014 by Hynek Lavička   *
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
#include <memory>
#include <vector>
#include <string>

#include <boost/python.hpp>
#include <boost/python/raw_function.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/concept_check.hpp>

#ifdef THREADING
#include <boost/thread.hpp>
#include <boost/date_time.hpp>
#include <boost/thread/shared_mutex.hpp>
#include <boost/thread/condition_variable.hpp>
#include <boost/thread/detail/thread.hpp>
#endif

#include "fluctuation_analyzer.h"
#include "detrended_fluctuation_analysis.cpp"


namespace py = boost::python;
using namespace boost::python;

template<int ...> struct seq{};
template<int N, int ...S> struct gens : gens<N-1, N-1, S...>{};
template<int ...S> struct gens<0, S...> {typedef seq<S...> type;};

template <typename ...Args>
struct cpptuple2pytuple_wrapper
{
	std::tuple<Args...> params;
	cpptuple2pytuple_wrapper(const std::tuple<Args...>& _params):params(_params){}
 
	py::tuple delayed_dispatch() 
	{
		return callFunc(typename gens<sizeof...(Args)>::type());
	}
 
	template<int ...S>
	py::tuple callFunc(seq<S...>)
	{
		return py::make_tuple(std::get<S>(params) ...);
	}
};

template <typename ...Args>
struct pytuple2cpptuple_wrapper
{
	py::tuple params;
	pytuple2cpptuple_wrapper(const py::tuple& _params):params(_params){}
 
	std::tuple<Args...> delayed_dispatch()
	{
		return callFunc(typename gens<sizeof...(Args)>::type());
	}
 
	template<int ...S>
	std::tuple<Args...> callFunc(seq<S...>)
	{
		return std::make_tuple((static_cast<Args>(py::extract<Args>(params[S])))...);
	}
};

// Convert (C++) tuple to (Python) tuple as PyObject*.
template<typename ... Args> 
PyObject* cpptuple2pytuple(const std::tuple<Args...>& t)
{
	cpptuple2pytuple_wrapper<Args...> wrapper(t);
	py::tuple bpt = wrapper.delayed_dispatch();
	return py::incref(py::object(bpt).ptr());
}

// Convert (Python) tuple to (C++) tuple.
template<typename ... Args> std::tuple<Args...> pytuple2cpptuple(PyObject* obj)
{
	py::tuple tup(py::borrowed(obj));
	pytuple2cpptuple_wrapper<Args...> wrapper(tup);
	std::tuple<Args...> bpt = wrapper.delayed_dispatch();
	return bpt;
}

template<typename ... Args>
struct cpptuple_to_python_tuple
{
	static PyObject* convert(const std::tuple<Args...>& t)
	{
		return cpptuple2pytuple<Args...>(t);
	}
};

template<typename ... Args>
struct cpptuple_from_python_tuple
{
	cpptuple_from_python_tuple()
	{
		py::converter::registry::push_back(&convertible, &construct, py::type_id<std::tuple<Args...> >());
	}
 
	static void * convertible ( PyObject * obj_ptr )
	{
		if (!PyTuple_CheckExact(obj_ptr))
		{
			return 0;
		}
		return obj_ptr;
	}

	static void construct ( PyObject * obj_ptr , py::converter::rvalue_from_python_stage1_data * data )
	{
		void * storage = ( ( py::converter::rvalue_from_python_storage<std::tuple<Args...> > * ) data ) -> storage.bytes;
		new (storage) std::tuple<Args...>(pytuple2cpptuple<Args...>(obj_ptr));
		data->convertible = storage;
	}
};

template<typename ...Args>
void create_tuple_converter()
{
	py::to_python_converter<std::tuple<Args...>, cpptuple_to_python_tuple<Args...> > ( );
	cpptuple_from_python_tuple<Args...> ( );
}


template<typename TYPE_INCOMMING , typename TYPE_OUTGOING , typename CARRIED_TYPE>
int convertor_to_CXX ( TYPE_INCOMMING income , TYPE_OUTGOING & outcome)
{
	int err_count = 0;
	std::size_t size_of_incomming_object = py::len (income);

	for ( std::size_t count = 0 ; count < size_of_incomming_object ; count ++ )
	{
		// conversion of each element to type CARRIED_TYPE
		py::extract<CARRIED_TYPE> element ( income[count] );
		if ( ! element.check() )
		{
			std::cout << "Dropped element at " << count << std::endl;
			err_count ++;
			continue;
		}

		outcome.push_back ( element() );
	}

	return err_count;
}

template<typename TYPE_INCOMMING>
int convertor_to_python ( TYPE_INCOMMING & income , py::list & outcome)
{
	typename TYPE_INCOMMING::iterator iter_end = income.end();
	for ( typename TYPE_INCOMMING::iterator iter = income.begin() ; iter != iter_end ; ++ iter )
	{
		outcome.append ( ( * iter ) );
	}

	return 0;
}

template<typename TYPE_INCOMMING>
std :: ostream & operator<< ( std :: ostream & aStream , std::vector<TYPE_INCOMMING> & storage )
{
	for ( typename std::vector<TYPE_INCOMMING>::iterator iter = storage.begin() ; iter != storage.end() ; ++ iter )
	{
		std::cout << ( * iter ) << " ";
	}

	return aStream;
}

void myprint ( double x )
{
	std::cout << x << " ";
}


template <class K, class V>
py::dict toPythonDict1 ( std::map<K, V> & map )
{
	py::dict dictionary;
	for (auto iter = map.begin(); iter != map.end(); ++iter)
	{
		dictionary [ py::make_tuple ( std::get<0> ( iter -> first ) , std::get<1> ( iter -> first ) , std::get<2> ( iter -> first ) ) ] = ( iter -> second );
	}
	return dictionary;
}


template <class K, class V>
py::dict toPythonDict4 ( std::map<K, V> & map )
{
	py::dict dictionary;
	for (auto iter = map.begin(); iter != map.end(); ++iter)
	{
		dictionary [ py::make_tuple ( std::get<0> ( iter -> first ) , std::get<1> ( iter -> first ) , std::get<2> ( iter -> first ) ) ] = py::make_tuple ( ( iter -> second ) [ 0 ] , ( iter -> second ) [ 1 ] , ( iter -> second ) [ 2 ] , ( iter -> second ) [ 3 ] );
	}
	return dictionary;
}

template<class TYPE1, class TYPE2, class TYPE3, class TYPE4>
tuple DFA ( tuple args , dict kargs )
{
	std::locale::global(std::locale(""));
	std::cout.imbue(std::locale());

	py::extract<py::list> data_list ( args [0] );
	py::extract<py::tuple> segment_size_tuple ( args [1] );
	py::extract<py::list> q_values_list ( args [2] );
	py::extract<long> order_number ( args [3] );
	py::extract<py::list> order_list ( args [3] );
	py::extract<bool> moving_window ( args [4] );
	py::extract<int> moving_window_size ( args [5] );
	py::extract<bool> double_integration ( args [6] );
	py::extract<bool> sample_variance ( args [7] );
	py::extract<unsigned int> number_of_threads ( args [8] );

	std::vector<double> q_values;
	std::vector<long> orders_of_method;

	if ( ( !data_list.check() ) && ( !segment_size_tuple.check() ) && ( !q_values_list.check() ) && ( ! ( order_number.check() || order_list.check() ) ) && ( !moving_window.check() ) && ( !moving_window_size.check() ) && ( !double_integration.check() ) && ( !sample_variance.check() ) )
	{
		std::cout << "Bad argument! Unable to convert incomming variables." << std::endl;
		PyErr_BadArgument();
	}
	else
	{
		if ( order_number.check() )
		{
			orders_of_method.push_back( order_number );
		}
		else
		{
			convertor_to_CXX<py::list , std::vector<long> , long> ( order_list () , orders_of_method );
		}
	}

	std::list<double> dataset;
	convertor_to_CXX<py::list , std::list<double> , double> ( data_list () , dataset );

	std::vector<unsigned int> segment_sizes;
	convertor_to_CXX<py::tuple , std::vector<unsigned int> , unsigned int> ( segment_size_tuple () , segment_sizes );

	std::map<std::tuple<double, unsigned int, unsigned int> , double > DFA_results;

	convertor_to_CXX<py::list , std::vector<double> , double> ( q_values_list () , q_values );
	std::cout << "Method start." << std::endl;
// 	for_each (dataset.begin(), dataset.end(), myprint );
	std::string output_filename_debug;

#ifdef THREADING
	std::cout << "Calculationg orders: " << orders_of_method << std::endl;

// 	statistical_functions::DetrendedFluctuationAnalysis ( dataset , segment_sizes , DFA_results , q_values , orders_of_method , output_filename_debug.c_str() , moving_window () , moving_window_size () , sample_variance () , double_integration () , number_of_threads () );

	statistical_functions::detrended_fluctuation_analysis::DetrendedFluctuationAnalysisTemplate <TYPE1, TYPE2 , TYPE3 , TYPE4> ( dataset , segment_sizes , DFA_results , q_values , orders_of_method , output_filename_debug.c_str() , moving_window () , moving_window_size () , sample_variance () , double_integration () , number_of_threads ()  );

	std::cout << "Method finished." << std::endl;
#endif

	py::dict results_dictionary;
	if ( DFA_results.empty() == true )
	{
		std::cerr << "The dictionaty is empty." << std::endl;

		return make_tuple ( results_dictionary );
	}
	else
	{
		// conversion of C++ dictionary to python dictionary - change neccessary in case of update of output 
		results_dictionary = toPythonDict1 ( DFA_results );

		DFA_results.clear ();

		return make_tuple ( results_dictionary );
	}
}

template<class TYPE1, class TYPE2, class TYPE3, class TYPE4>
tuple DCCA ( tuple args , dict kargs )
{
	std::locale::global(std::locale(""));
	std::cout.imbue(std::locale());

	py::extract<py::list> data_list1 ( args [0] );
	py::extract<py::list> data_list2 ( args [1] );
	py::extract<py::tuple> segment_size_tuple ( args [2] );
	py::extract<py::list> q_values_list ( args [3] );
	py::extract<long> order_number ( args [4] );
	py::extract<py::list> order_list ( args [4] );
	py::extract<bool> moving_window ( args [5] );
	py::extract<int> moving_window_size ( args [6] );
	py::extract<bool> double_integration ( args [7] );
	py::extract<bool> sample_variance ( args [8] );
	py::extract<unsigned int> number_of_threads ( args [9] );

	std::vector<double> q_values;
	std::vector<unsigned int> orders_of_method;

	if ( ( !data_list1.check() ) && ( !data_list2.check() ) && ( !segment_size_tuple.check() ) && ( !q_values_list.check() ) && ( ! ( order_number.check() || order_list.check() ) ) && ( !moving_window.check() ) && ( !moving_window_size.check() ) && ( !double_integration.check() ) && ( !sample_variance.check() ) )
	{
		std::cout << "Bad argument! Unable to convert incomming variables." << std::endl;
		PyErr_BadArgument();
	}
	else
	{
		if ( order_number.check() )
		{
			orders_of_method.push_back( order_number );
		}
		else
		{
			convertor_to_CXX<py::list , std::vector<unsigned int> , unsigned int> ( order_list () , orders_of_method );
		}
	}

	std::list<double> dataset1;
	convertor_to_CXX<py::list , std::list<double> , double> ( data_list1 () , dataset1 );

	std::list<double> dataset2;
	convertor_to_CXX<py::list , std::list<double> , double> ( data_list2 () , dataset2 );

	std::vector<unsigned int> segment_sizes;
	convertor_to_CXX<py::tuple , std::vector<unsigned int> , unsigned int> ( segment_size_tuple () , segment_sizes );

	std::list<std::shared_ptr< std::vector<std::list<double> > > > DFA_complete_results;

	convertor_to_CXX<py::list , std::vector<double> , double> ( q_values_list () , q_values );
	std::cout << "Method start." << std::endl;
// 	for_each (dataset.begin(), dataset.end(), myprint );
	std::string output_filename_debug;

#ifdef THREADING
	std::cout << "Calculationg orders: " << orders_of_method << std::endl;
	std::map<std::tuple<double, unsigned int, unsigned int> , std::array<double,4> > DCCA_results;

// 	statistical_functions::DetrendedCrossCorrelationAnalysis ( dataset1 , dataset2 , segment_sizes , DCCA_results , q_values , orders_of_method , output_filename_debug.c_str() , moving_window () , moving_window_size () , sample_variance () , double_integration () , number_of_threads () );
	statistical_functions::detrended_fluctuation_analysis::DetrendedCrossCorrelationAnalysisTemplate <TYPE1, TYPE2, TYPE3 , TYPE4> ( dataset1 , dataset2 , segment_sizes , DCCA_results , q_values , orders_of_method ,  output_filename_debug.c_str() , moving_window () , moving_window_size () , sample_variance () , double_integration () , number_of_threads () );

	std::cout << "Method finished." << std::endl;

	py::dict results_dictionary;
	if ( DCCA_results.empty() == true )
	{
		std::cerr << "The dictionaty is empty." << std::endl;

		return make_tuple ( results_dictionary );
	}
	else
	{
		// conversion of C++ dictionary to python dictionary - change neccessary in case of update of output 
		results_dictionary = toPythonDict4 ( DCCA_results );

		DCCA_results.clear ();

		return make_tuple ( results_dictionary );
	}
#endif
}

BOOST_PYTHON_MODULE(pydfa)
{
// exposed functions
	py::def("DFA", DFA<double,double,double,double> );
	py::def("DFA_100_30", py::raw_function ( DFA<double,boost::multiprecision::number<boost::multiprecision::cpp_bin_float<100> >,boost::multiprecision::number<boost::multiprecision::cpp_bin_float<30> >,double> ) );
	py::def("DFA_150_30", py::raw_function ( DFA<double,boost::multiprecision::number<boost::multiprecision::cpp_bin_float<150> >,boost::multiprecision::number<boost::multiprecision::cpp_bin_float<30> >,double> ) );
	py::def("DFA_200_30", py::raw_function ( DFA<double,boost::multiprecision::number<boost::multiprecision::cpp_bin_float<200> >,boost::multiprecision::number<boost::multiprecision::cpp_bin_float<30> >,double> ) );
	py::def("DFA_250_30", py::raw_function ( DFA<double,boost::multiprecision::number<boost::multiprecision::cpp_bin_float<250> >,boost::multiprecision::number<boost::multiprecision::cpp_bin_float<30> >,double> ) );
	py::def("DFA_300_30", py::raw_function ( DFA<double,boost::multiprecision::number<boost::multiprecision::cpp_bin_float<300> >,boost::multiprecision::number<boost::multiprecision::cpp_bin_float<30> >,double> ) );

	py::def("DCCA", py::raw_function ( DCCA<double,double,double,double> ) );
	py::def("DCCA_100_30", py::raw_function ( DCCA<double,boost::multiprecision::number<boost::multiprecision::cpp_bin_float<100> >,boost::multiprecision::number<boost::multiprecision::cpp_bin_float<30> >,double> ) );
	py::def("DCCA_150_30", py::raw_function ( DCCA<double,boost::multiprecision::number<boost::multiprecision::cpp_bin_float<150> >,boost::multiprecision::number<boost::multiprecision::cpp_bin_float<30> >,double> ) );
	py::def("DCCA_200_30", py::raw_function ( DCCA<double,boost::multiprecision::number<boost::multiprecision::cpp_bin_float<200> >,boost::multiprecision::number<boost::multiprecision::cpp_bin_float<30> >,double> ) );
	py::def("DCCA_250_30", py::raw_function ( DCCA<double,boost::multiprecision::number<boost::multiprecision::cpp_bin_float<250> >,boost::multiprecision::number<boost::multiprecision::cpp_bin_float<30> >,double> ) );
	py::def("DCCA_300_30", py::raw_function ( DCCA<double,boost::multiprecision::number<boost::multiprecision::cpp_bin_float<300> >,boost::multiprecision::number<boost::multiprecision::cpp_bin_float<30> >,double> ) );
}
