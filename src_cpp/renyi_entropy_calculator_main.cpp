
#include <iostream>
#include <ranges>
#include <unordered_map>
#include <stacktrace>
#include <algorithm>
#include <functional>
#include <print>
#include <type_traits>

#include <cpptrace/from_current.hpp>
#include <cpptrace/formatting.hpp>
#include <cpptrace/cpptrace.hpp>

#include <boost/stacktrace.hpp>
#include <boost/exception/all.hpp>
#include <boost/program_options.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/p_square_quantile.hpp>
#include <boost/accumulators/statistics/tail_quantile.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/tail_quantile.hpp>
#include <boost/accumulators/framework/accumulator_set.hpp>
#include <boost/accumulators/statistics/tail.hpp>
#include <boost/accumulators/statistics_fwd.hpp>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/zstd.hpp>
#include <boost/iostreams/filter/lzma.hpp>
#include <boost/iostreams/filter/bzip2.hpp>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <Python.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "msgpack.hpp"

#include "serialize_tuple.h"
#include "random_samples.h"
#include "renyi_entropy.h"
#include "utils.h"
#include "cpptrace_helper.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
#include "pandas_wrapper.h"
#pragma GCC diagnostic pop

// required to stacktrace
typedef boost::error_info<struct tag_stacktrace, boost::stacktrace::stacktrace> traced;

template <class E>
void throw_with_trace ( const E& e )
{
    throw boost::enable_error_info ( e )
            << traced ( boost::stacktrace::stacktrace() );
}

namespace po = boost::program_options;
namespace bio = boost::iostreams;
namespace py = pybind11;

// compression filter
struct zstd_ostream : boost::iostreams::filtering_ostream
{
    zstd_ostream ( std::ostream& os )
    {
        bio::zstd_params zstd_params ( 13 );
        bio::filtering_ostream::push ( bio::zstd_compressor{zstd_params} );
        bio::filtering_ostream::push ( os );
    }
};

const auto microseconds_in_second =
    static_cast<double> ( std::chrono::duration_cast<std::chrono::microseconds> ( 1s ).count() );

int main ( int argc, char *argv[] )
{
    cpptrace::absorb_trace_exceptions ( false );
    cpptrace::use_default_stderr_logger();
    cpptrace::register_terminate_handler();
    warmup_cpptrace();
    segfault_handler_cpptrace ();

    std::string directory, output, model_name;
    unsigned int maximal_neighborhood, number_of_datapoints;
    unsigned int start_neighbor = 5; // initial neighbour for averaging over neighbours
    unsigned int runs, runs_surrogate, runs_shuffle;
    unsigned int number_of_threads;
    double mean1, mean2, mean3, std1, std2, std3, alpha_interaction, beta_interaction, eta_interaction, theta_interaction, gamma_interaction;
    std::vector<unsigned int> selection_x, selection_y, selection_z;

    typedef double calculation_type;
    std::vector<std::vector<unsigned int>> history_firsts, future_firsts,
        history_seconds, history_thirds;
    std::vector<calculation_type> alpha;
    std::string alphas_string, history_first_string, future_first_string,
        history_second_string, history_third_string;
    po::options_description desc ( "Allowed options" );
    desc.add_options() ( "help,h", "produce help message" )
    ( "directory,d", po::value<std::string>()->default_value ( "." ), "Folder to export results" )
    ( "file,f", po::value<std::string>()->default_value ( "CRE.bin.zstd" ), "Output file" )
    ( "history_first", po::value<std::string>()->composing(), "History of the first timeries X" )
    ( "future_first", po::value<std::string>()->composing(), "Future of the first timeries X" )
    ( "history_second", po::value<std::string>()->composing(), "History of the second timeseries Y" )
    ( "history_third", po::value<std::string>()->composing(), "History of the third timeseries Z" )
    ( "runs,r", po::value<unsigned int>()->default_value ( 1U ), "Number of runs" )
    ( "runs_shuffle", po::value<unsigned int>(), "Number shuffle runs" )
    ( "runs_surrogate", po::value<unsigned int>(), "Number surrogate runs" )
    ( "maximal_neighborhood", po::value<unsigned int>(), "Maximal neighborhood" )
    ( "alpha", po::value<std::string>()->composing(), "Renyi entropy parameter" )
    ( "multithreading", "Multithreading" )
    ( "number_of_threads", po::value<unsigned int>()->default_value ( std::thread::hardware_concurrency() ), "Number of thread to user" )
    ( "model", po::value<std::string>()->default_value ( "gaussian" ), "Model to investigate [gaussian, AB, CB, ADB, ACB]" )
    ( "number_of_datapoints", po::value<unsigned int>()->default_value ( 1000U ), "Size of timeseries to analyze" )

    ( "mean1", po::value<double>()->default_value ( 0. ), "Mean 1 of random noise" )
    ( "mean2", po::value<double>()->default_value ( 0. ), "Mean 2 of random noise" )
    ( "mean3", po::value<double>()->default_value ( 0. ), "Mean 3 of random noise" )
    ( "std1", po::value<double>()->default_value ( 1. ), "Standard deviation 1 of random noise" )
    ( "std2", po::value<double>()->default_value ( 1. ), "Standard deviation 2 of random noise" )
    ( "std3", po::value<double>()->default_value ( 0.5 ), "Standard deviation 3 of random noise" )
    ( "alpha_interaction", po::value<double>()->default_value ( 0.5 ), "Alpha interaction" )
    ( "beta_interaction", po::value<double>()->default_value ( 0.2 ), "Beta interaction" )
    ( "gamma_interaction", po::value<double>()->default_value ( 0.5 ), "Gamma interaction" )
    ( "eta_interaction", po::value<double>()->default_value ( 0.5 ), "Eta interaction" )
    ( "theta_interaction", po::value<double>()->default_value ( 0.2 ), "Theta interaction" )

    ( "selection_x", po::value<std::vector<unsigned int>> ( &selection_x )->multitoken(), "Selection of X" )
    ( "selection_y", po::value<std::vector<unsigned int>> ( &selection_y )->multitoken(), "Selection of Y" )
    ( "selection_z", po::value<std::vector<unsigned int>> ( &selection_z )->multitoken(), "Selection of Z" )
    ;

    po::variables_map vm;
    po::store ( po::parse_command_line ( argc, argv, desc ), vm );
    po::notify ( vm );

    CPPTRACE_TRY
    {
        if ( vm.count ( "help" ) )
        {
            std::cout << desc << std::endl;
            cpptrace::generate_trace().print();
            return 1;
        }

        bool multithreading = false;
        if ( vm.count ( "multithreading" ) )
        {
            multithreading = true;
        }
        if ( vm.count ( "directory" ) )
        {
            directory = vm["directory"].as<std::string>();
        }
        output = vm["file"].as<std::string>();
        if ( vm.count ( "history_first" ) )
        {
            history_first_string = vm["history_first"].as<std::string>();
            chop_string_arrays<std::string, unsigned int> (
                history_first_string, ',', history_firsts,
                [] ( const std::string str, std::size_t *pos ) -> unsigned int
            {
                return std::stoul ( str, pos, 10 );
            } );
        }
        else
        {
            std::cerr << "history_first option must be used" << std::endl;
        }
        if ( vm.count ( "future_first" ) )
        {
            future_first_string = vm["future_first"].as<std::string>();
            chop_string_arrays<std::string, unsigned int> (
                future_first_string, ',', future_firsts,
                [] ( const std::string str, std::size_t *pos ) -> unsigned int
            {
                return std::stoul ( str, pos, 10 );
            } );
        }
        else
        {
            std::cerr << "future_first option must be used" << std::endl;
        }
        if ( vm.count ( "history_second" ) )
        {
            history_second_string = vm["history_second"].as<std::string>();
            chop_string_arrays<std::string, unsigned int> (
                history_second_string, ',', history_seconds,
                [] ( const std::string str, std::size_t *pos ) -> unsigned int
            {
                return std::stoul ( str, pos, 10 );
            } );
        }
        else
        {
            std::cerr << "history_second option must be used" << std::endl;
        }
        if ( vm.count ( "history_third" ) )
        {
            history_third_string = vm["history_third"].as<std::string>();
            chop_string_arrays<std::string, unsigned int> (
                history_third_string, ',', history_thirds,
                [] ( const std::string str, std::size_t *pos ) -> unsigned int
            {
                return std::stoul ( str, pos, 10 );
            } );
        }
        else
        {
            std::cerr << "history_third option must be used" << std::endl;
            history_thirds.push_back({{}});
        }
        if ( vm.count ( "maximal_neighborhood" ) )
        {
            maximal_neighborhood = vm["maximal_neighborhood"].as<unsigned int>();
        }
        if ( vm.count ( "alpha" ) )
        {
            alphas_string = vm["alpha"].as<std::string>();
            convert_array<std::string, calculation_type> ( alphas_string, alpha, std::stod );
        }
        runs = vm["runs"].as<unsigned int>();
        number_of_datapoints = vm["number_of_datapoints"].as<unsigned int>();
        mean1 = vm["mean1"].as<double>();
        mean2 = vm["mean2"].as<double>();
        mean3 = vm["mean3"].as<double>();
        std1 = vm["std1"].as<double>();
        std2 = vm["std2"].as<double>();
        std3 = vm["std3"].as<double>();
        alpha_interaction = vm["alpha_interaction"].as<double>();
        beta_interaction = vm["beta_interaction"].as<double>();
        eta_interaction = vm["eta_interaction"].as<double>();
        theta_interaction = vm["theta_interaction"].as<double>();
        gamma_interaction = vm["gamma_interaction"].as<double>();

        if ( vm.count ( "runs_shuffle" ) )
        {
            runs_shuffle = vm["runs_shuffle"].as<unsigned int>();
        }
        else
        {
            runs_shuffle = runs;
        }

        if ( vm.count ( "runs_surrogate" ) )
        {
            runs_surrogate = vm["runs_surrogate"].as<unsigned int>();
        }
        else
        {
            runs_surrogate = runs;
        }
        number_of_threads = vm["number_of_threads"].as<unsigned int>();
        model_name = vm["model"].as<std::string>();

        if ( vm.count ( "selection_x" ) )
        {
            if ( selection_x.size() == 0 )
            {
                std::cerr << "selection_x: nonzero values needed\n";
                return 1;
            }
            std::cout << "selection x: ";
            for(const auto & item: selection_x)
            {
                std::cout << item << " ";
            }
            std::cout << std::endl;
        }
        else
        {
            std::cerr << "selection_x: nonzero values needed\n";
            return 1;
        }

        if ( vm.count ( "selection_y" ) )
        {
            if ( selection_y.size() == 0 )
            {
                std::cerr << "selection_y: nonzero values needed\n";
                return 1;
            }
            std::cout << "selection y: ";
            for(const auto & item: selection_y)
            {
                std::cout << item << " ";
            }
            std::cout << std::endl;
        }
        else
        {
            std::cerr << "selection_y: nonzero values needed\n";
            return 1;
        }

        {
            std::cout << "selection z: ";
            for(const auto & item: selection_z)
            {
                std::cout << item << " ";
            }
            std::cout << std::endl;
        }

        unsigned int postselection_X_future;
        unsigned int postselection_X_history;
        unsigned int postselection_Y_history;
        auto indices_to_use_range = std::ranges::iota_view{1U, maximal_neighborhood};

        renyi_entropy::renyi_entropy<calculation_type>::collection_result_conditional_information_transfer_type collection_result_RTE;
        renyi_entropy::renyi_entropy<calculation_type>::result_RTE_t collection_result_RTEs;
        renyi_entropy::renyi_entropy<calculation_type>::type_average_result_conditional_information_transfer_type processing_RTE;

        std::cout << "model: ";
        for ( auto &item : alpha )
        {
            std::cout << item << " ";
        }
        std::cout << std::endl;

#ifdef PYTHON_EXPORT
        python_wrappers::pandas_wrapper pandas;
        std::list<py::tuple> columns_of_table;
        const auto column_multiindex = pandas.pandas_multiindex_from_tuples ( "tuples"_a=columns_of_table, "names"_a=std::list<py::str>{py::str ( "swap_datasets" ), py::str ( "shuffle_indicator" ), py::str ( "surrogate_indicator" ), py::str ( "sample" ), py::str ( "future_first" ), py::str ( "histories_first" ), py::str ( "histories_second" ), py::str ( "renyi_type" ), py::str ( "neighbor" ) } );

        auto result_dataframe = pandas.pandas_dataframe ( "columns"_a=column_multiindex );
#endif
        Eigen::MatrixXd dataset;
        auto start_timeseries_construction = std::chrono::high_resolution_clock::now();
        {
            std::vector<std::pair<std::string, double>> output{{"mean1", mean1}, {"mean2", mean2}, {"mean3", mean3}, {"std1", std1}, {"std2", std2}, {"std3", std3}, {"alpha_interaction", alpha_interaction}, {"beta_interaction", beta_interaction}, {"gamma_interaction", gamma_interaction}, {"eta_interaction", eta_interaction}, {"theta_interaction", theta_interaction}};
            std::cout << "----------- " << std::endl;
            std::cout << "model: " << model_name << std::endl;
            for(auto &item: output)
            {
                std::cout << item.first << " " << item.second << std::endl;
            }
             std::cout << "----------- " << std::endl;

            std::random_device random_device{};
            std::mt19937 generator{random_device() };
            std::normal_distribution distribution1{mean1, std1};
            std::normal_distribution distribution2{mean2, std2};
            std::normal_distribution distribution3{mean3, std3};
            auto distribution_generator1 = [&]()
            {
                return distribution1 ( generator );
            };
            auto distribution_generator2 = [&]()
            {
                return distribution2 ( generator );
            };
            auto distribution_generator3 = [&]()
            {
                return distribution3 ( generator );
            };

            std::vector<std::tuple<std::string, std::function<void () >>> random_model{
                {
                    "gaussian", [&]()
                    {
                        random_samples::samples_normal_distribution_uncorrelated ( dataset, mean1, std1, number_of_datapoints, 2 );
                    }
                },
                {
                    "AB", [&]()
                    {
                        random_samples::AB_process ( dataset, alpha_interaction, beta_interaction, eta_interaction, number_of_datapoints, std::make_tuple ( distribution_generator1, distribution_generator2 ) );
                    }
                },
                {
                    "CB", [&]()
                    {
                        random_samples::CB_process ( dataset, alpha_interaction, beta_interaction, theta_interaction, number_of_datapoints, std::make_tuple ( distribution_generator1, distribution_generator2 ) );
                    }
                },
                {
                    "ADB", [&]()
                    {
                        random_samples::ADB_process ( dataset, alpha_interaction, beta_interaction, gamma_interaction, eta_interaction, theta_interaction, number_of_datapoints, std::make_tuple ( distribution_generator1, distribution_generator2, distribution_generator3 ) );
                    }
                },
                {
                    "ACB", [&]()
                    {
                        random_samples::ACB_process ( dataset, alpha_interaction, beta_interaction, gamma_interaction, eta_interaction, theta_interaction, number_of_datapoints, std::make_tuple ( distribution_generator1, distribution_generator2, distribution_generator3 ) );
                    }
                }
            };

            auto result = std::find_if (
                              random_model.begin(),
                              random_model.end(),
                              [&] ( const auto &item )
            {
                if ( std::get<0> ( item ) == model_name )
                {
                    std::get<1> ( item ) ();
                    return true;
                }
                else return false;
            } );
            if ( result == random_model.end() )
            {
                std::cout << "Unable to find model" << std::endl;
                exit ( 1 );
            }
        }
        auto end_timeseries_construction =
            std::chrono::high_resolution_clock::now();
        auto elapsed_timeseries_construction =
            end_timeseries_construction - start_timeseries_construction;
        auto microseconds_timeseries_construction =
            std::chrono::duration_cast<std::chrono::microseconds> ( elapsed_timeseries_construction ).count();
        BOOST_LOG_TRIVIAL ( trace ) << "Timeseries construction " << microseconds_timeseries_construction / microseconds_in_second << " seconds";

        auto start_RTE_calculation = std::chrono::high_resolution_clock::now();
        // calculation of raw values
        for ( auto swap_datasets :
    {
        false, true
    } )
        {
            for ( auto & type_of_average:
                {
                    std::tuple<bool,bool> ( {false, false} ), std::tuple<bool,bool> ( {true, false} ), std::tuple<bool,bool> ( {false, true} )
                } )
            {
                auto [shuffle_indicator, surrogate_indicator] = type_of_average;

                unsigned int maximal_samples = ( ( shuffle_indicator == false ) && ( surrogate_indicator == false ) ) ? 1U : ( ( shuffle_indicator == true ) && ( surrogate_indicator == false ) ? runs_shuffle : runs_surrogate );
                // calculation of shuffled datasets
                for ( unsigned int sample = 0; sample < maximal_samples; ++ sample )
                {
                    std::cout << std::format ( "Swapped: {} Shuffled: {} Surrogate: {} Run: {}", swap_datasets, shuffle_indicator, surrogate_indicator, sample ) << std::endl;
#ifndef NDEBUG
                    std::cout << "Original dataset: " << std::endl << dataset << std::endl;
#endif
                    auto start_dataset_preparation = std::chrono::high_resolution_clock::now();
                    auto [dataset1, dataset2, dataset3] =
                        renyi_entropy::renyi_entropy<calculation_type>::prepare_dataset (
                            dataset, swap_datasets, shuffle_indicator, surrogate_indicator, selection_x, selection_y, selection_z );
#ifndef NDEBUG
                    std::cout << "Original dataset: " << std::endl << dataset << std::endl << "dataset1: " << std::endl << dataset1 << std::endl << "dataset2: " << std::endl << dataset2 << std::endl << "dataset3: " << dataset3 << std::endl;
#endif
                    auto end_dataset_preparation = std::chrono::high_resolution_clock::now();
                    auto elapsed_dataset_preparation = end_dataset_preparation - start_dataset_preparation;
                    auto microseconds_dataset_preparation = std::chrono::duration_cast<std::chrono::microseconds> ( elapsed_dataset_preparation ).count();
                    BOOST_LOG_TRIVIAL ( trace ) << "Dataset preparation " << microseconds_dataset_preparation / microseconds_in_second << " seconds";

                    for ( auto &future_first : future_firsts )
                    {
                        for ( auto &histories_first : history_firsts )
                        {
                            for ( auto &histories_second : history_seconds )
                            {
                                for ( auto &histories_third : history_thirds )
                                {
                                    renyi_entropy::renyi_entropy<calculation_type>::collection_conditional_information_transfer_key_type collection_key = std::make_tuple ( future_first, histories_first, histories_second );
                                    std::map<std::string, std::any> configuration_prepare_dataset =
                                    {
                                        std::pair<std::string, std::any> ( "transpose", true ),
                                        std::pair<std::string, std::any> ( "history_index_x", histories_first ),
                                        std::pair<std::string, std::any> ( "history_index_y", histories_second ),
                                        std::pair<std::string, std::any> ( "history_index_z", histories_third ),
                                        std::pair<std::string, std::any> ( "future_index_x", future_first ),
                                        std::pair<std::string, std::any> ( "postselection_y_fut", postselection_X_future ),
                                        std::pair<std::string, std::any> ( "postselection_z_hist", postselection_Y_history ),
                                        std::pair<std::string, std::any> ( "postselection_y_hist", postselection_X_history )
                                    };
                                    auto [y_future, y_history, x_history, z_history] = renyi_entropy::
                                                                            renyi_entropy<calculation_type>::PrepareDatasetForTransferEntropy (
                                                                                    dataset1, dataset2, dataset3, configuration_prepare_dataset );
#ifndef NDEBUG
                                    std::cout << "y_fut: " << y_future << std::endl << "y_hist: " << y_history << std::endl << "x_hist: " << x_history << std::endl << "z_hist: " << z_history << std::endl;
#endif
                                    std::vector<unsigned int> indices_to_use ( indices_to_use_range.begin(), indices_to_use_range.end() );
                                    std::map<std::string, std::any> configuration_renyi_entropy =
                                    {
                                        std::pair<std::string, std::any> ( "transpose", true ),
                                        std::pair<std::string, std::any> ( "axis_to_join", 0 ),
                                        std::pair<std::string, std::any> ( "method", std::string("LeonenkoProzanto") ),
                                        std::pair<std::string, std::any> ( "alphas", alpha ),
                                        std::pair<std::string, std::any> ( "enhanced_calculation", true ),
                                        std::pair<std::string, std::any> ( "indices", indices_to_use ),
                                        std::pair<std::string, std::any> ( "multithreading", multithreading ) // true
                                    };

                                    renyi_entropy::renyi_entropy<calculation_type>::renyi_conditional_information_transfer (
                                        collection_result_RTEs,
                                        y_future,
                                        y_history,
                                        x_history,
                                        z_history,
                                        configuration_renyi_entropy,
                                        [&] ( std::string type_of_calculation, unsigned int neighbor )->renyi_entropy::renyi_entropy<calculation_type>::result_RTE_key
                                    {
                                        return std::make_tuple ( swap_datasets, shuffle_indicator, surrogate_indicator, sample, future_first, histories_first, histories_second, type_of_calculation, neighbor );
                                    }
                                    );

#ifdef PYTHON_EXPORT
                                    for ( const auto & RTE_type : renyi_entropy::RTE_types )
                                    {
                                        auto new_column_specifier = py::make_tuple ( py::bool_ ( swap_datasets ), py::bool_ ( shuffle_indicator ), py::bool_ ( surrogate_indicator ), py::int_ ( sample ), py::tuple ( py::cast ( future_first ) ), py::tuple ( py::cast ( histories_first ) ), py::tuple ( py::cast ( histories_second ) ), py::str ( *const_cast<std::string *> ( RTE_type ) ) );
                                        columns_of_table.push_back ( new_column_specifier );

                                        auto dictionary = py::dict();
                                        for ( auto &[key, value]: ( * ( transfer_entropy_results.begin() ) ).second )
                                        {
                                            auto [neighbor, alpha] = key;

                                            dictionary[py::make_tuple ( py::int_ ( neighbor ), py::float_ ( alpha ) )] = py::float_ ( value );
                                        }
                                        auto empty_series = pandas.pandas_series ( dictionary );
                                        result_dataframe[new_column_specifier] = empty_series;
                                        py::print ( result_dataframe );
                                        //auto empty_dataframe = pandas.pandas_dataframe();
                                        //empty_dataframe[new_column_specifier] = ;
                                    }
#endif
                                }
                            }
                        }
                    }
                }
            }
        }
        auto end_RTE_calculation =
            std::chrono::high_resolution_clock::now();
        auto elapsed_RTE_calculation =
            end_RTE_calculation - start_RTE_calculation;
        auto microseconds_RTE_calculation =
            std::chrono::duration_cast<std::chrono::microseconds> ( elapsed_RTE_calculation ).count();
        BOOST_LOG_TRIVIAL ( trace ) << "Result calculation " << microseconds_RTE_calculation / microseconds_in_second << " seconds";

        auto start_result_processing = std::chrono::high_resolution_clock::now();

#ifdef PYTHON_EXPORT
        std::list<py::tuple> indices_of_table;
        const auto &single_result = ( * ( * ( * collection_result_RTE.begin() ).second.begin() ).second.begin() ).second;
        for ( const auto &[key, value]: single_result )
        {
            const auto &[nn_index, alpha] = key;
            indices_of_table.push_back ( py::make_tuple ( py::int_ ( nn_index ), py::float_ ( alpha ) ) );
        }
        const auto index_multiindex = pandas.pandas_multiindex_from_tuples ( "tuples"_a=indices_of_table, "names"_a=std::list<py::str>{py::str ( "index_of_neighbor" ), py::str ( "alpha" ) } );
        auto empty_dataframe = pandas.pandas_dataframe ( "index"_a=index_multiindex, "columns"_a=column_multiindex );
        py::print ( empty_dataframe );
#endif
        auto end_result_processing =
            std::chrono::high_resolution_clock::now();
        auto elapsed_result_processing =
            end_result_processing - start_result_processing;
        auto microseconds_result_processing =
            std::chrono::duration_cast<std::chrono::microseconds> ( elapsed_result_processing ).count();
        BOOST_LOG_TRIVIAL ( trace ) << "Result processing " << microseconds_result_processing / microseconds_in_second << " seconds";

        boost::filesystem::path output_file =
            boost::filesystem::path ( directory ) / boost::filesystem::path ( output );
        boost::filesystem::ofstream output_file_handler ( output_file );
        zstd_ostream zstd_compression_stream{output_file_handler};

        bool msg_pack = false;
        // save results
        if ( msg_pack == true )
        {
            std::stringstream ss;
            msgpack::pack ( ss, collection_result_RTEs ); //renyi_entropy::renyi_entropy<calculation_type>::storage_RTE(collection_result_RTE, processing_RTE)
            std::cout << ss.str().size() << std::endl;
            zstd_compression_stream << ss.str();
        }
        else
        {
            boost::archive::binary_oarchive output_archive(zstd_compression_stream);
            output_archive << collection_result_RTEs;
        }

        auto end_result_storing =
            std::chrono::high_resolution_clock::now();
        auto elapsed_result_storing =
            end_result_storing - end_result_processing;
        auto microseconds_result_storing =
            std::chrono::duration_cast<std::chrono::microseconds> ( elapsed_result_storing ).count();
        BOOST_LOG_TRIVIAL ( trace ) << "Result storing " << microseconds_result_storing / microseconds_in_second << " seconds";
    }
    CPPTRACE_CATCH ( std::exception& exc )
    {
        std::cerr << "Exception: "<< exc.what() << std::endl;
        cpptrace::from_current_exception().print();
    }
}
// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
