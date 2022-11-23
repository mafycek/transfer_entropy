#include <list>
#include <tuple>

#include <exception>
#include <stdexcept>

#include <algorithm>
#include <chrono>
#include <iostream>
#include <functional>

#include "math_functions.h"
#include "detrended_fluctuation_analysis.h"

namespace statistical_functions
{

namespace detrended_fluctuation_analysis
{

void enqueue_thread_for_cleanup_template ( )
{
	boost::mutex::scoped_lock lock(DFASynchronization::GetDFASynchronization().GetAddJobsMutex () );

	DFASynchronization::GetDFASynchronization().GetFinishedJobs().push_back ( boost::this_thread::get_id() );

	DFASynchronization::GetDFASynchronization().GetAddJobsCV().notify_one ();
}

}

}
