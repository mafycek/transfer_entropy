
#include <iostream>

#include <msgpack.hpp>
#include <matplotlibcpp.h>

#include <boost/stacktrace.hpp>
#include <boost/exception/all.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include <cpptrace/from_current.hpp>
#include <cpptrace/formatting.hpp>
#include <cpptrace/cpptrace.hpp>

#include "cpptrace_helper.h"
#include "renyi_entropy.h"

namespace po = boost::program_options;
namespace plt = matplotlibcpp;

int main ( int argc, char *argv[] )
{
	typedef double calculation_type;
	cpptrace::absorb_trace_exceptions(false);
	cpptrace::use_default_stderr_logger();
	cpptrace::register_terminate_handler();
	warmup_cpptrace();
	segfault_handler_cpptrace ();

	po::options_description desc ( "Allowed options" );
	desc.add_options() ( "help, h", "produce help message" )
	("directory, d", po::value<std::string>()->default_value("."), "Folder to export results" )
	("file,f", po::value<std::string>()->default_value("CRE.bin"), "Output file" );

	po::variables_map vm;
	po::store ( po::parse_command_line ( argc, argv, desc ), vm );
	po::notify ( vm );

	CPPTRACE_TRY
	{
		std::string output, directory;


		if ( vm.count ( "help" ) )
		{
				std::cout << desc << std::endl;
				cpptrace::generate_trace().print();
				return 1;
		}
		if ( vm.count ( "directory" ) )
		{
				directory = vm["directory"].as<std::string>();
		}
		output = vm["file"].as<std::string>();

		std::string str;
		boost::filesystem::path output_file =
		boost::filesystem::path(directory) / boost::filesystem::path(output);
		boost::filesystem::ifstream input_file_handler ( output_file );

		if (input_file_handler.is_open())
		{
			std::string buffer(std::istreambuf_iterator<char>{input_file_handler}, {});
			input_file_handler >> buffer;
			msgpack::object_handle object_handle_buffer = msgpack::unpack(buffer.data(), buffer.size());

			// deserialized object is valid during the msgpack::object_handle instance is alive.
			msgpack::object deserialized = object_handle_buffer.get();

			std::tuple<renyi_entropy::renyi_entropy<calculation_type>::collection_result_conditional_information_transfer_type, renyi_entropy::renyi_entropy<calculation_type>::type_average_result_conditional_information_transfer_type> result_conditional_information_transfer;

			// deserialized.convert(result_conditional_information_transfer);
			// auto processed_results = std::get<1>(result_conditional_information_transfer);

			plt::InitializeInterpreter();
			plt::useBackend("TkCairo");
			std::cout << "Backend: " << plt::getUsedBackend() << std::endl;
			plt::rcParams(std::string("text.usetex"), true);
			plt::xlabel("$Time\\ \\alpha \\beta$");
			plt::ylabel("Student confusion");
			plt::title("Standard usage"); // set a title
			plt::legend();
			plt::show();
			//plt::savefig("minimal.pdf");
		}
		else
		{
			std::cerr << "Cannot open input file" << std::endl;
		}
	}
	CPPTRACE_CATCH (std::exception& exc)
	{
		std::cerr << "Exception: "<< exc.what() << std::endl;
		cpptrace::from_current_exception().print();
	}

}
