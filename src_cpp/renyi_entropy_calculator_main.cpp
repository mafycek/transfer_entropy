
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

#include <Python.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "msgpack.hpp"

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
void throw_with_trace(const E& e) {
    throw boost::enable_error_info(e)
        << traced(boost::stacktrace::stacktrace());
}

namespace po = boost::program_options;
namespace bio = boost::iostreams;
namespace py = pybind11;

struct zstd_ostream : boost::iostreams::filtering_ostream
{
    zstd_ostream(std::ostream& os)
    {
        bio::zstd_params zstd_params ( 13 );
        bio::filtering_ostream::push( bio::zstd_compressor{zstd_params} );
        bio::filtering_ostream::push( os );
    }
};

const auto microseconds_in_second =
        static_cast<double> (std::chrono::duration_cast<std::chrono::microseconds> (1s).count());

int main ( int argc, char *argv[] )
{
    cpptrace::absorb_trace_exceptions(false);
    cpptrace::use_default_stderr_logger();
    cpptrace::register_terminate_handler();
    warmup_cpptrace();
    segfault_handler_cpptrace ();

    std::string directory, output, model_name;
    unsigned int maximal_neighborhood, number_of_datapoints;
    unsigned int start_neighbor = 5; // initial neighbour for averaging over neighbours
    unsigned int runs, runs_surrogate, runs_shuffle;
    unsigned int number_of_threads;
    double mean1, mean2, mean3, std1, std2, std3;
    typedef double calculation_type;
    std::vector<std::vector<unsigned int>> history_firsts, future_firsts,
        history_seconds;
    std::vector<calculation_type> alpha;
    std::string alphas_string, history_first_string, future_first_string,
        history_second_string;
    po::options_description desc ( "Allowed options" );
    desc.add_options() ( "help, h", "produce help message" ) 
    ("directory,d", po::value<std::string>()->default_value("."), "Folder to export results" )
    ("file,f", po::value<std::string>()->default_value("CRE.bin.zstd"), "Output file" )
    ("history_first", po::value<std::string>()->composing(), "History of the first timeries" ) 
    ("future_first", po::value<std::string>()->composing(), "Future of the first timeries" ) 
    ("history_second", po::value<std::string>()->composing(), "History of the second timeseries" )
    ("runs,r", po::value<unsigned int>()->default_value(1U), "Number of runs" )
    ("runs_shuffle", po::value<unsigned int>(), "Number shuffle runs" )
    ("runs_surrogate", po::value<unsigned int>(), "Number surrogate runs" )
    ("maximal_neighborhood", po::value<unsigned int>(), "Maximal neighborhood" )
    ("alpha", po::value<std::string>()->composing(), "Renyi entropy parameter" )
    ("multithreading", "Multithreading" )
    ("number_of_threads", po::value<unsigned int>()->default_value(std::thread::hardware_concurrency()), "Number of thread to user" )
    ("model", po::value<std::string>()->default_value("gaussian"), "Model to investigate [gaussian, AB, CB, ADB, ACB]" )
    ("number_of_datapoints", po::value<unsigned int>()->default_value(1000U), "Size of timeseries to analyze" )
    ("mean1", po::value<double>()->default_value(0.), "Mean 1 of random noise" )
    ("mean2", po::value<double>()->default_value(0.), "Mean 2 of random noise" )
    ("mean3", po::value<double>()->default_value(0.), "Mean 3 of random noise" )
    ("std1", po::value<double>()->default_value(1.), "Standard deviation 1 of random noise" )
    ("std2", po::value<double>()->default_value(1.), "Standard deviation 2 of random noise" )
    ("std3", po::value<double>()->default_value(0.25), "Standard deviation 3 of random noise" )
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

        unsigned int postselection_X_future;
        unsigned int postselection_X_history;
        unsigned int postselection_Y_history;
        auto indices_to_use_range = std::ranges::iota_view{1U, maximal_neighborhood};

        renyi_entropy::renyi_entropy<calculation_type>::collection_result_conditional_information_transfer_type collection_result_RTE;
        renyi_entropy::renyi_entropy<calculation_type>::result_RTE_t collection_result_RTEs;
        renyi_entropy::renyi_entropy<calculation_type>::type_average_result_conditional_information_transfer_type processing_RTE;

        for ( auto &item : alpha )
        {
            std::cout << item << " ";
        }
        std::cout << std::endl;

#ifdef PYTHON_EXPORT
        python_wrappers::pandas_wrapper pandas;
        std::list<py::tuple> columns_of_table;
        const auto column_multiindex = pandas.pandas_multiindex_from_tuples("tuples"_a=columns_of_table, "names"_a=std::list<py::str>{py::str("swap_datasets"), py::str("shuffle_indicator"), py::str("surrogate_indicator"), py::str("sample"), py::str("future_first"), py::str("histories_first"), py::str("histories_second"), py::str("renyi_type"), py::str("neighbor")});

        auto result_dataframe = pandas.pandas_dataframe("columns"_a=column_multiindex);
#endif
        Eigen::MatrixXd dataset;
        auto start_timeseries_construction = std::chrono::high_resolution_clock::now();
        {
            std::random_device random_device{};
            std::mt19937 generator{random_device()};
            std::normal_distribution distribution1{mean1, std1};
            std::normal_distribution distribution2{mean2, std2};
            std::normal_distribution distribution3{mean3, std3};
            auto distrbution_generator1 = [&]()
            {
                return distribution1(generator);
            };
            auto distrbution_generator2 = [&]()
            {
                return distribution2(generator);
            };
            auto distrbution_generator3 = [&]()
            {
                return distribution3(generator);
            };

            std::vector<std::tuple<std::string, std::function<void ()>>> random_model{
                {"gaussian", [&](){ random_samples::samples_normal_distribution_uncorrelated ( dataset, mean1, std1, number_of_datapoints, 2 ); } },
                {"AB", [&](){ random_samples::AB_process(dataset, 0.5, 0.2, 1, number_of_datapoints, std::make_tuple(distrbution_generator1, distrbution_generator2));}},
                {"CB", [&](){ random_samples::CB_process(dataset, 0.5, 0.2, 1, number_of_datapoints, std::make_tuple(distrbution_generator1, distrbution_generator2));}},
                {"ADB", [&](){ random_samples::ADB_process(dataset, 0.5, 0.5, 0.2, 0.2, 0.5, number_of_datapoints, std::make_tuple(distrbution_generator1, distrbution_generator2, distrbution_generator3));}},
                {"ACB", [&](){ random_samples::ACB_process(dataset, 0.5, 0.5, 0.2, 0.2, 0.5, number_of_datapoints, std::make_tuple(distrbution_generator1, distrbution_generator2, distrbution_generator3));}}
            };

            auto result = std::find_if(
                random_model.begin(),
                random_model.end(),
                [&](const auto &item)
                {
                    if (std::get<0>(item) == model_name )
                    {
                        std::get<1>(item)();
                        return true;
                    }
                    else return false;
                });
            if (result == random_model.end())
            {
                std::cout << "Unable to find model" << std::endl;
                exit(1);
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
                true, false
            } )
        {
            for (auto & type_of_average: {std::tuple<bool,bool>({false, false}), std::tuple<bool,bool>({true, false}), std::tuple<bool,bool>({false, true}) })
            {
                auto [shuffle_indicator, surrogate_indicator] = type_of_average;

                unsigned int maximal_samples = (( shuffle_indicator == false) && (surrogate_indicator == false )) ? 1U : ( ( shuffle_indicator == true) && (surrogate_indicator == false ) ? runs_shuffle : runs_surrogate );
                // calculation of shuffled datasets
                for (unsigned int sample = 0; sample < maximal_samples; ++ sample)
                {
                    auto start_dataset_preparation = std::chrono::high_resolution_clock::now();
                    auto [dataset1, dataset2] =
                        renyi_entropy::renyi_entropy<calculation_type>::prepare_dataset (
                            dataset, swap_datasets, shuffle_indicator, surrogate_indicator, 1, 1 );
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
                                renyi_entropy::renyi_entropy<calculation_type>::collection_conditional_information_transfer_key_type collection_key = std::make_tuple(future_first, histories_first, histories_second);
                                std::map<std::string, std::any> configuration_prepare_dataset =
                                {
                                    std::pair<std::string, std::any> ( "transpose", true ),
                                    std::pair<std::string, std::any> ( "history_index_x", histories_first ),
                                    std::pair<std::string, std::any> ( "history_index_y", histories_second ),
                                    std::pair<std::string, std::any> ( "future_index_x", future_first ),
                                    std::pair<std::string, std::any> ( "postselection_y_fut", postselection_X_future ),
                                    std::pair<std::string, std::any> ( "postselection_z_hist", postselection_Y_history ),
                                    std::pair<std::string, std::any> ( "postselection_y_hist", postselection_X_history )
                                };
                                auto [y_future, y_history, z_history] = renyi_entropy::
                                                                        renyi_entropy<calculation_type>::PrepareDatasetForTransferEntropy (
                                                                                dataset1, dataset2, configuration_prepare_dataset );

                                std::vector<unsigned int> indices_to_use (indices_to_use_range.begin(), indices_to_use_range.end());
                                std::map<std::string, std::any> configuration_renyi_entropy =
                                {
                                    std::pair<std::string, std::any> ( "transpose", true ),
                                    std::pair<std::string, std::any> ( "axis_to_join", 0 ),
                                    std::pair<std::string, std::any> ( "method", "LeonenkoProzanto" ),
                                    std::pair<std::string, std::any> ( "alphas", alpha ),
                                    std::pair<std::string, std::any> ( "enhanced_calculation", true ),
                                    std::pair<std::string, std::any> ( "indices", indices_to_use ),
                                    std::pair<std::string, std::any> ( "multithreading", true )
                                };

                                renyi_entropy::renyi_entropy<calculation_type>::renyi_conditional_information_transfer (
                                    collection_result_RTEs,
                                    y_future,
                                    y_history,
                                    z_history,
                                    configuration_renyi_entropy,
                                    [&](std::string type_of_calculation, unsigned int neighbor)->renyi_entropy::renyi_entropy<calculation_type>::result_RTE_key
                                    {
                                        return std::make_tuple(swap_datasets, shuffle_indicator, surrogate_indicator, sample, future_first, histories_first, histories_second, type_of_calculation, neighbor);
                                    }
                                );

#ifdef PYTHON_EXPORT
                                for (const auto & RTE_type : renyi_entropy::RTE_types)
                                {
                                    auto new_column_specifier = py::make_tuple(py::bool_(swap_datasets), py::bool_(shuffle_indicator), py::bool_(surrogate_indicator), py::int_(sample), py::tuple(py::cast(future_first)), py::tuple(py::cast(histories_first)), py::tuple(py::cast(histories_second)), py::str(*const_cast<std::string *>(RTE_type)) );
                                    columns_of_table.push_back(new_column_specifier);

                                    auto dictionary = py::dict();
                                    for (auto &[key, value]: (*(transfer_entropy_results.begin())).second)
                                    {
                                        auto [neighbor, alpha] = key;

                                        dictionary[py::make_tuple(py::int_(neighbor), py::float_(alpha))] = py::float_(value);
                                    }
                                    auto empty_series = pandas.pandas_series(dictionary);
                                    result_dataframe[new_column_specifier] = empty_series;
                                    py::print(result_dataframe);
                                    //auto empty_dataframe = pandas.pandas_dataframe();
                                    //empty_dataframe[new_column_specifier] = ;
                                }
#endif
                            }
                        }
                    }
                    std::cout << std::format("Sw: {} Sh:{} Sur:{} R:{}", swap_datasets, shuffle_indicator, surrogate_indicator, sample) << std::endl;
                }
            }
        }
        auto end_RTE_calculation =
            std::chrono::high_resolution_clock::now();
        auto elapsed_RTE_calculation =
            end_RTE_calculation - start_RTE_calculation;
        auto microseconds_RTE_calculation =
            std::chrono::duration_cast<std::chrono::microseconds> ( elapsed_RTE_calculation ).count();
        BOOST_LOG_TRIVIAL ( trace ) << "Result colculation " << microseconds_RTE_calculation / microseconds_in_second << " seconds";

#ifdef PYTHON_EXPORT
        std::list<py::tuple> indices_of_table;
        const auto &single_result = (*(*(* collection_result_RTE.begin()).second.begin()).second.begin()).second;
        for(const auto &[key, value]: single_result)
        {
            const auto &[nn_index, alpha] = key;
            indices_of_table.push_back(py::make_tuple(py::int_(nn_index), py::float_(alpha)));
        }
        const auto index_multiindex = pandas.pandas_multiindex_from_tuples("tuples"_a=indices_of_table, "names"_a=std::list<py::str>{py::str("index_of_neighbor"), py::str("alpha")} );
        auto empty_dataframe = pandas.pandas_dataframe("index"_a=index_multiindex, "columns"_a=column_multiindex);
        py::print(empty_dataframe);
#endif

        auto start_result_processing = std::chrono::high_resolution_clock::now();
        std::set<std::tuple<bool, bool, bool, std::vector<unsigned int>, std::vector<unsigned int>, std::vector<unsigned int>>> main_RTE_titles;
        std::set<std::tuple<bool, bool, bool, std::vector<unsigned int>, std::vector<unsigned int>, std::vector<unsigned int>>> ballance_RTE_titles;
        std::set<std::tuple<bool, bool, bool, std::vector<unsigned int>, std::vector<unsigned int>, std::vector<unsigned int>>> ballance_effective_RTE_titles;

        {
            auto collection_key_views = std::views::keys(collection_result_RTEs);
            std::vector<decltype(collection_result_RTEs)::key_type> collection_keys{ collection_key_views.begin(), collection_key_views.end() };
            for(auto & collection_key: collection_keys)
            {
                const auto &[swap_datasets, shuffle_indicator, surrogate_indicator, sample, future_first, histories_first, histories_second, type_of_calculation, neighbor] = collection_key;
                if (type_of_calculation == renyi_entropy::conditional_renyi_entropy_label)
                {
                    if (((shuffle_indicator == true ) && (surrogate_indicator == false)) || ((shuffle_indicator == false ) && (surrogate_indicator == true)))
                    {
                        main_RTE_titles.insert(std::make_tuple(swap_datasets, shuffle_indicator, surrogate_indicator, future_first, histories_first, histories_second));

                    }
                    if (((shuffle_indicator == false ) && (surrogate_indicator == false)) || ((shuffle_indicator == false ) && (surrogate_indicator == false)))
                    {
                        if (swap_datasets == false)
                        {
                            ballance_RTE_titles.insert(std::make_tuple(swap_datasets, shuffle_indicator, surrogate_indicator, future_first, histories_first, histories_second));
                        }
                    }
                    if (((shuffle_indicator == true ) && (surrogate_indicator == false)) || ((shuffle_indicator == false ) && (surrogate_indicator == true)))
                    {
                        if (swap_datasets == false)
                        {
                            ballance_effective_RTE_titles.insert(std::make_tuple(swap_datasets, shuffle_indicator, surrogate_indicator, future_first, histories_first, histories_second));
                        }
                    }
                }
            }
        }

        // effective RTE calculation, run statistics
        for(auto & main_RTE_title: main_RTE_titles)
        {
            const auto [swap_datasets, shuffle_indicator, surrogate_indicator, future_first, histories_first, histories_second] = main_RTE_title;
            const auto type_of_calculation = renyi_entropy::conditional_renyi_entropy_label;
            for (const auto neighbor: indices_to_use_range)
            {
                const auto collection_key = std::make_tuple(swap_datasets, false, false, 0, future_first, histories_first, histories_second, type_of_calculation, neighbor);
                auto & main_conditional_renyi_entropy_results = collection_result_RTEs[collection_key];

                const auto CRE_key_views = std::views::keys(main_conditional_renyi_entropy_results);
                std::vector<renyi_entropy::renyi_entropy<double>::renyi_key_type> CRE_keys{ CRE_key_views.begin(), CRE_key_views.end() };
                for (auto & CRE_key: CRE_keys)
                {
                    using stats_accumulators = boost::accumulators::stats<boost::accumulators::tag::mean, boost::accumulators::tag::variance, boost::accumulators::tag::median, boost::accumulators::tag::tail_quantile<boost::accumulators::left>>;
                    boost::accumulators::accumulator_set< double, stats_accumulators > accumulator_statistics(boost::accumulators::tag::tail<boost::accumulators::left>::cache_size = 100);

                    for ( unsigned int run = 0; run < runs; ++ run )
                    {
                        const auto collection_key = std::make_tuple(swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, type_of_calculation, neighbor);
                        auto & counterpart_conditional_renyi_entropy_results = collection_result_RTEs[collection_key];

                        auto entropy_difference = main_conditional_renyi_entropy_results[CRE_key] - counterpart_conditional_renyi_entropy_results[CRE_key];
                        accumulator_statistics ( entropy_difference );
                    }
                    auto mean = boost::accumulators::mean(accumulator_statistics);
                    auto variance = boost::accumulators::variance(accumulator_statistics);
                    auto median = boost::accumulators::median(accumulator_statistics);
                    auto quantile_01 = boost::accumulators::quantile(accumulator_statistics, boost::accumulators::quantile_probability = 0.1);
                    auto quantile_02 = boost::accumulators::quantile(accumulator_statistics, boost::accumulators::quantile_probability = 0.2);
                    auto quantile_03 = boost::accumulators::quantile(accumulator_statistics, boost::accumulators::quantile_probability = 0.3);
                    auto quantile_04 = boost::accumulators::quantile(accumulator_statistics, boost::accumulators::quantile_probability = 0.4);
                    auto quantile_06 = boost::accumulators::quantile(accumulator_statistics, boost::accumulators::quantile_probability = 0.6);
                    auto quantile_07 = boost::accumulators::quantile(accumulator_statistics, boost::accumulators::quantile_probability = 0.7);
                    auto quantile_08 = boost::accumulators::quantile(accumulator_statistics, boost::accumulators::quantile_probability = 0.8);
                    auto quantile_09 = boost::accumulators::quantile(accumulator_statistics, boost::accumulators::quantile_probability = 0.9);
                    const auto final_CRE_key_mean = std::make_tuple(false, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "mean", renyi_entropy::average_runs);
                    const auto final_CRE_key_standard_deviation = std::make_tuple(false, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "standard_deviation", renyi_entropy::average_runs);
                    const auto final_CRE_key_median = std::make_tuple(false, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "median", renyi_entropy::average_runs);
                    const auto final_CRE_key_quantile_01 = std::make_tuple(false, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_01", renyi_entropy::average_runs);
                    const auto final_CRE_key_quantile_02 = std::make_tuple(false, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_02", renyi_entropy::average_runs);
                    const auto final_CRE_key_quantile_03 = std::make_tuple(false, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_03", renyi_entropy::average_runs);
                    const auto final_CRE_key_quantile_04 = std::make_tuple(false, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_04", renyi_entropy::average_runs);
                    const auto final_CRE_key_quantile_06 = std::make_tuple(false, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_06", renyi_entropy::average_runs);
                    const auto final_CRE_key_quantile_07 = std::make_tuple(false, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_07", renyi_entropy::average_runs);
                    const auto final_CRE_key_quantile_08 = std::make_tuple(false, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_08", renyi_entropy::average_runs);
                    const auto final_CRE_key_quantile_09 = std::make_tuple(false, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_09", renyi_entropy::average_runs);

                    processing_RTE[final_CRE_key_mean][CRE_key] = mean;
                    processing_RTE[final_CRE_key_standard_deviation][CRE_key] = variance;
                    processing_RTE[final_CRE_key_median][CRE_key] = median;
                    processing_RTE[final_CRE_key_quantile_01][CRE_key] = quantile_01;
                    processing_RTE[final_CRE_key_quantile_02][CRE_key] = quantile_02;
                    processing_RTE[final_CRE_key_quantile_03][CRE_key] = quantile_03;
                    processing_RTE[final_CRE_key_quantile_04][CRE_key] = quantile_04;
                    processing_RTE[final_CRE_key_quantile_06][CRE_key] = quantile_06;
                    processing_RTE[final_CRE_key_quantile_07][CRE_key] = quantile_07;
                    processing_RTE[final_CRE_key_quantile_08][CRE_key] = quantile_08;
                    processing_RTE[final_CRE_key_quantile_09][CRE_key] = quantile_09;
                }
            }
        }

        // ballance calculation, run statistics
        for(auto & ballance_RTE_title: ballance_RTE_titles)
        {
            const auto [swap_datasets, shuffle_indicator, surrogate_indicator, future_first, histories_first, histories_second] = ballance_RTE_title;
            const auto type_of_calculation = renyi_entropy::conditional_renyi_entropy_label;

            for (const auto neighbor: indices_to_use_range)
            {
                const auto collection_key = std::make_tuple(swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, type_of_calculation, neighbor);
                auto & main_conditional_renyi_entropy_results = collection_result_RTEs[collection_key];

                const auto CRE_key_views = std::views::keys(main_conditional_renyi_entropy_results);
                std::vector<renyi_entropy::renyi_entropy<double>::renyi_key_type> CRE_keys{ CRE_key_views.begin(), CRE_key_views.end() };
                for (auto & CRE_key: CRE_keys)
                {
                    using stats_accumulators = boost::accumulators::stats<boost::accumulators::tag::mean, boost::accumulators::tag::variance, boost::accumulators::tag::median, boost::accumulators::tag::tail_quantile<boost::accumulators::left>>;
                    boost::accumulators::accumulator_set< double, stats_accumulators > accumulator_statistics(boost::accumulators::tag::tail<boost::accumulators::left>::cache_size = 100);

                    for ( unsigned int run = 0; run < runs; ++ run )
                    {
                        const auto collection_key = std::make_tuple(swap_datasets, shuffle_indicator, surrogate_indicator, run, future_first, histories_first, histories_second, type_of_calculation, neighbor);
                        auto & main_conditional_renyi_entropy_results = collection_result_RTEs[collection_key];

                        const auto collection_reverse_direction_key = std::make_tuple(!swap_datasets, shuffle_indicator, surrogate_indicator, run, future_first, histories_first, histories_second, type_of_calculation, neighbor);
                        auto & main_reverse_conditional_renyi_entropy_results = collection_result_RTEs[collection_reverse_direction_key];

                        auto entropy_difference = main_conditional_renyi_entropy_results[CRE_key] - main_reverse_conditional_renyi_entropy_results[CRE_key];
                        accumulator_statistics ( entropy_difference );
                    }

                    auto mean = boost::accumulators::mean(accumulator_statistics);
                    auto variance = boost::accumulators::variance(accumulator_statistics);
                    auto median = boost::accumulators::median(accumulator_statistics);
                    auto quantile_01 = boost::accumulators::quantile(accumulator_statistics, boost::accumulators::quantile_probability = 0.1);
                    auto quantile_02 = boost::accumulators::quantile(accumulator_statistics, boost::accumulators::quantile_probability = 0.2);
                    auto quantile_03 = boost::accumulators::quantile(accumulator_statistics, boost::accumulators::quantile_probability = 0.3);
                    auto quantile_04 = boost::accumulators::quantile(accumulator_statistics, boost::accumulators::quantile_probability = 0.4);
                    auto quantile_06 = boost::accumulators::quantile(accumulator_statistics, boost::accumulators::quantile_probability = 0.6);
                    auto quantile_07 = boost::accumulators::quantile(accumulator_statistics, boost::accumulators::quantile_probability = 0.7);
                    auto quantile_08 = boost::accumulators::quantile(accumulator_statistics, boost::accumulators::quantile_probability = 0.8);
                    auto quantile_09 = boost::accumulators::quantile(accumulator_statistics, boost::accumulators::quantile_probability = 0.9);
                    const auto final_CRE_key_mean = std::make_tuple(true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "mean", renyi_entropy::average_runs);
                    const auto final_CRE_key_standard_deviation = std::make_tuple(true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "standard_deviation", renyi_entropy::average_runs);
                    const auto final_CRE_key_median = std::make_tuple(true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "median", renyi_entropy::average_runs);
                    const auto final_CRE_key_quantile_01 = std::make_tuple(true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_01", renyi_entropy::average_runs);
                    const auto final_CRE_key_quantile_02 = std::make_tuple(true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_02", renyi_entropy::average_runs);
                    const auto final_CRE_key_quantile_03 = std::make_tuple(true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_03", renyi_entropy::average_runs);
                    const auto final_CRE_key_quantile_04 = std::make_tuple(true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_04", renyi_entropy::average_runs);
                    const auto final_CRE_key_quantile_06 = std::make_tuple(true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_06", renyi_entropy::average_runs);
                    const auto final_CRE_key_quantile_07 = std::make_tuple(true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_07", renyi_entropy::average_runs);
                    const auto final_CRE_key_quantile_08 = std::make_tuple(true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_08", renyi_entropy::average_runs);
                    const auto final_CRE_key_quantile_09 = std::make_tuple(true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_09", renyi_entropy::average_runs);

                    processing_RTE[final_CRE_key_mean][CRE_key] = mean;
                    processing_RTE[final_CRE_key_standard_deviation][CRE_key] = variance;
                    processing_RTE[final_CRE_key_median][CRE_key] = median;
                    processing_RTE[final_CRE_key_quantile_01][CRE_key] = quantile_01;
                    processing_RTE[final_CRE_key_quantile_02][CRE_key] = quantile_02;
                    processing_RTE[final_CRE_key_quantile_03][CRE_key] = quantile_03;
                    processing_RTE[final_CRE_key_quantile_04][CRE_key] = quantile_04;
                    processing_RTE[final_CRE_key_quantile_06][CRE_key] = quantile_06;
                    processing_RTE[final_CRE_key_quantile_07][CRE_key] = quantile_07;
                    processing_RTE[final_CRE_key_quantile_08][CRE_key] = quantile_08;
                    processing_RTE[final_CRE_key_quantile_09][CRE_key] = quantile_09;
                }
            }
        }

        // ballance effective calculation, run statistics
        for(auto & ballance_effective_RTE_title: ballance_effective_RTE_titles)
        {
            const auto [swap_datasets, shuffle_indicator, surrogate_indicator, future_first, histories_first, histories_second] = ballance_effective_RTE_title;
            const auto type_of_calculation = renyi_entropy::conditional_renyi_entropy_label;

            for (const auto neighbor: indices_to_use_range)
            {
                const auto collection_key = std::make_tuple(swap_datasets, false, false, 0, future_first, histories_first, histories_second, type_of_calculation, neighbor);
                auto & main_conditional_renyi_entropy_results = collection_result_RTEs[collection_key];
                const auto swap_collection_key = std::make_tuple(!swap_datasets, false, false, 0, future_first, histories_first, histories_second, type_of_calculation, neighbor);
                auto & main_reverse_conditional_renyi_entropy_results = collection_result_RTEs[swap_collection_key];

                const auto CRE_key_views = std::views::keys(main_conditional_renyi_entropy_results);
                std::vector<renyi_entropy::renyi_entropy<double>::renyi_key_type> CRE_keys{ CRE_key_views.begin(), CRE_key_views.end() };
                for (auto & CRE_key: CRE_keys)
                {
                    using stats_accumulators = boost::accumulators::stats<boost::accumulators::tag::mean, boost::accumulators::tag::variance, boost::accumulators::tag::median, boost::accumulators::tag::tail_quantile<boost::accumulators::left>>;
                    boost::accumulators::accumulator_set< double, stats_accumulators > accumulator_statistics(boost::accumulators::tag::tail<boost::accumulators::left>::cache_size = 100);

                    for ( unsigned int run = 0; run < runs; ++ run )
                    {
                        const auto collection_key = std::make_tuple(swap_datasets, shuffle_indicator, surrogate_indicator, run, future_first, histories_first, histories_second, type_of_calculation, neighbor);
                        auto & counterpart_conditional_renyi_entropy_results = collection_result_RTEs[collection_key];

                        const auto collection_reverse_direction_key = std::make_tuple(!swap_datasets, shuffle_indicator, surrogate_indicator, run, future_first, histories_first, histories_second, type_of_calculation, neighbor);
                        auto & counterpart_reverse_conditional_renyi_entropy_results = collection_result_RTEs[collection_reverse_direction_key];

                        auto entropy_difference = main_conditional_renyi_entropy_results[CRE_key] - counterpart_conditional_renyi_entropy_results[CRE_key] - main_reverse_conditional_renyi_entropy_results[CRE_key] + counterpart_reverse_conditional_renyi_entropy_results[CRE_key];
                        accumulator_statistics ( entropy_difference );
                    }

                    auto mean = boost::accumulators::mean(accumulator_statistics);
                    auto variance = boost::accumulators::variance(accumulator_statistics);
                    auto median = boost::accumulators::median(accumulator_statistics);
                    auto quantile_01 = boost::accumulators::quantile(accumulator_statistics, boost::accumulators::quantile_probability = 0.1);
                    auto quantile_02 = boost::accumulators::quantile(accumulator_statistics, boost::accumulators::quantile_probability = 0.2);
                    auto quantile_03 = boost::accumulators::quantile(accumulator_statistics, boost::accumulators::quantile_probability = 0.3);
                    auto quantile_04 = boost::accumulators::quantile(accumulator_statistics, boost::accumulators::quantile_probability = 0.4);
                    auto quantile_06 = boost::accumulators::quantile(accumulator_statistics, boost::accumulators::quantile_probability = 0.6);
                    auto quantile_07 = boost::accumulators::quantile(accumulator_statistics, boost::accumulators::quantile_probability = 0.7);
                    auto quantile_08 = boost::accumulators::quantile(accumulator_statistics, boost::accumulators::quantile_probability = 0.8);
                    auto quantile_09 = boost::accumulators::quantile(accumulator_statistics, boost::accumulators::quantile_probability = 0.9);
                    const auto final_CRE_key_mean = std::make_tuple(true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "mean", renyi_entropy::average_runs);
                    const auto final_CRE_key_standard_deviation = std::make_tuple(true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "standard_deviation", renyi_entropy::average_runs);
                    const auto final_CRE_key_median = std::make_tuple(true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "median", renyi_entropy::average_runs);
                    const auto final_CRE_key_quantile_01 = std::make_tuple(true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_01", renyi_entropy::average_runs);
                    const auto final_CRE_key_quantile_02 = std::make_tuple(true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_02", renyi_entropy::average_runs);
                    const auto final_CRE_key_quantile_03 = std::make_tuple(true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_03", renyi_entropy::average_runs);
                    const auto final_CRE_key_quantile_04 = std::make_tuple(true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_04", renyi_entropy::average_runs);
                    const auto final_CRE_key_quantile_06 = std::make_tuple(true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_06", renyi_entropy::average_runs);
                    const auto final_CRE_key_quantile_07 = std::make_tuple(true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_07", renyi_entropy::average_runs);
                    const auto final_CRE_key_quantile_08 = std::make_tuple(true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_08", renyi_entropy::average_runs);
                    const auto final_CRE_key_quantile_09 = std::make_tuple(true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_09", renyi_entropy::average_runs);

                    processing_RTE[final_CRE_key_mean][CRE_key] = mean;
                    processing_RTE[final_CRE_key_standard_deviation][CRE_key] = variance;
                    processing_RTE[final_CRE_key_median][CRE_key] = median;
                    processing_RTE[final_CRE_key_quantile_01][CRE_key] = quantile_01;
                    processing_RTE[final_CRE_key_quantile_02][CRE_key] = quantile_02;
                    processing_RTE[final_CRE_key_quantile_03][CRE_key] = quantile_03;
                    processing_RTE[final_CRE_key_quantile_04][CRE_key] = quantile_04;
                    processing_RTE[final_CRE_key_quantile_06][CRE_key] = quantile_06;
                    processing_RTE[final_CRE_key_quantile_07][CRE_key] = quantile_07;
                    processing_RTE[final_CRE_key_quantile_08][CRE_key] = quantile_08;
                    processing_RTE[final_CRE_key_quantile_09][CRE_key] = quantile_09;
                }
            }
        }

        // effective RTE calculation, run and neighbor statistics
        for(auto & main_RTE_title: main_RTE_titles)
        {
            const auto [swap_datasets, shuffle_indicator, surrogate_indicator, future_first, histories_first, histories_second] = main_RTE_title;
            const auto type_of_calculation = renyi_entropy::conditional_renyi_entropy_label;

            const auto collection_key = std::make_tuple(swap_datasets, false, false, 0, future_first, histories_first, histories_second, type_of_calculation, indices_to_use_range[0]);
            auto & main_conditional_renyi_entropy_results = collection_result_RTEs[collection_key];
            const auto CRE_key_views = std::views::keys(main_conditional_renyi_entropy_results);
            std::vector<renyi_entropy::renyi_entropy<double>::renyi_key_type> CRE_keys{ CRE_key_views.begin(), CRE_key_views.end() };
            for (auto & CRE_key: CRE_keys)
            {
                using stats_accumulators = boost::accumulators::stats<boost::accumulators::tag::mean, boost::accumulators::tag::variance, boost::accumulators::tag::median, boost::accumulators::tag::tail_quantile<boost::accumulators::left>>;
                boost::accumulators::accumulator_set< double, stats_accumulators >
                //cache needs to contain complete dataset = #runs * #neighbors
                accumulator_statistics(boost::accumulators::tag::tail<boost::accumulators::left>::cache_size = 1090);

                int i = 0;
                for (const auto neighbor: indices_to_use_range)
                {
                    const auto collection_key = std::make_tuple(swap_datasets, false, false, 0, future_first, histories_first, histories_second, type_of_calculation, neighbor);
                    auto & main_conditional_renyi_entropy_results = collection_result_RTEs[collection_key];

                    for ( unsigned int run = 0; run < runs; ++ run )
                    {
                        const auto collection_key = std::make_tuple(swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, type_of_calculation, neighbor);
                        auto & counterpart_conditional_renyi_entropy_results = collection_result_RTEs[collection_key];

                        auto entropy_difference = main_conditional_renyi_entropy_results[CRE_key] - counterpart_conditional_renyi_entropy_results[CRE_key];
                        accumulator_statistics ( entropy_difference );
                        ++ i;
                    }
                }
                auto mean = boost::accumulators::mean(accumulator_statistics);
                auto variance = boost::accumulators::variance(accumulator_statistics);
                auto median = boost::accumulators::median(accumulator_statistics);
                auto quantile_01 = boost::accumulators::quantile(accumulator_statistics, boost::accumulators::quantile_probability = 0.1);
                auto quantile_02 = boost::accumulators::quantile(accumulator_statistics, boost::accumulators::quantile_probability = 0.2);
                auto quantile_03 = boost::accumulators::quantile(accumulator_statistics, boost::accumulators::quantile_probability = 0.3);
                auto quantile_04 = boost::accumulators::quantile(accumulator_statistics, boost::accumulators::quantile_probability = 0.4);
                auto quantile_06 = boost::accumulators::quantile(accumulator_statistics, boost::accumulators::quantile_probability = 0.6);
                auto quantile_07 = boost::accumulators::quantile(accumulator_statistics, boost::accumulators::quantile_probability = 0.7);
                auto quantile_08 = boost::accumulators::quantile(accumulator_statistics, boost::accumulators::quantile_probability = 0.8);
                auto quantile_09 = boost::accumulators::quantile(accumulator_statistics, boost::accumulators::quantile_probability = 0.9);
                const auto final_CRE_key_mean = std::make_tuple(false, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "mean", renyi_entropy::average_runs_neighbors);
                const auto final_CRE_key_standard_deviation = std::make_tuple(false, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "standard_deviation", renyi_entropy::average_runs_neighbors);
                const auto final_CRE_key_median = std::make_tuple(false, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "median", renyi_entropy::average_runs_neighbors);
                const auto final_CRE_key_quantile_01 = std::make_tuple(false, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "quantile_01", renyi_entropy::average_runs_neighbors);
                const auto final_CRE_key_quantile_02 = std::make_tuple(false, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "quantile_02", renyi_entropy::average_runs_neighbors);
                const auto final_CRE_key_quantile_03 = std::make_tuple(false, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "quantile_03", renyi_entropy::average_runs_neighbors);
                const auto final_CRE_key_quantile_04 = std::make_tuple(false, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "quantile_04", renyi_entropy::average_runs_neighbors);
                const auto final_CRE_key_quantile_06 = std::make_tuple(false, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "quantile_06", renyi_entropy::average_runs_neighbors);
                const auto final_CRE_key_quantile_07 = std::make_tuple(false, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "quantile_07", renyi_entropy::average_runs_neighbors);
                const auto final_CRE_key_quantile_08 = std::make_tuple(false, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "quantile_08", renyi_entropy::average_runs_neighbors);
                const auto final_CRE_key_quantile_09 = std::make_tuple(false, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "quantile_09", renyi_entropy::average_runs_neighbors);

                processing_RTE[final_CRE_key_mean][CRE_key] = mean;
                processing_RTE[final_CRE_key_standard_deviation][CRE_key] = variance;
                processing_RTE[final_CRE_key_median][CRE_key] = median;
                processing_RTE[final_CRE_key_quantile_01][CRE_key] = quantile_01;
                processing_RTE[final_CRE_key_quantile_02][CRE_key] = quantile_02;
                processing_RTE[final_CRE_key_quantile_03][CRE_key] = quantile_03;
                processing_RTE[final_CRE_key_quantile_04][CRE_key] = quantile_04;
                processing_RTE[final_CRE_key_quantile_06][CRE_key] = quantile_06;
                processing_RTE[final_CRE_key_quantile_07][CRE_key] = quantile_07;
                processing_RTE[final_CRE_key_quantile_08][CRE_key] = quantile_08;
                processing_RTE[final_CRE_key_quantile_09][CRE_key] = quantile_09;
            }
        }
        auto end_result_processing =
            std::chrono::high_resolution_clock::now();
        auto elapsed_result_processing =
            end_result_processing - start_result_processing;
        auto microseconds_result_processing =
            std::chrono::duration_cast<std::chrono::microseconds> ( elapsed_result_processing ).count();
        BOOST_LOG_TRIVIAL ( trace ) << "Result processing " << microseconds_result_processing / microseconds_in_second << " seconds";

        // save results
        std::stringstream ss;
        msgpack::pack(ss, std::make_tuple(collection_result_RTE, processing_RTE)); //renyi_entropy::renyi_entropy<calculation_type>::storage_RTE(collection_result_RTE, processing_RTE)
        std::cout << ss.str().size() << std::endl;
        boost::filesystem::path output_file =
        boost::filesystem::path(directory) / boost::filesystem::path(output);
        boost::filesystem::ofstream output_file_handler ( output_file );
        zstd_ostream zstd_compression_stream{output_file_handler};

        zstd_compression_stream << ss.str();

        auto end_result_storing =
            std::chrono::high_resolution_clock::now();
        auto elapsed_result_storing =
            end_result_storing - end_result_processing;
        auto microseconds_result_storing =
            std::chrono::duration_cast<std::chrono::microseconds> ( elapsed_result_storing ).count();
        BOOST_LOG_TRIVIAL ( trace ) << "Result storing " << microseconds_result_storing / microseconds_in_second << " seconds";
    }
    CPPTRACE_CATCH (std::exception& exc)
    {
        std::cerr << "Exception: "<< exc.what() << std::endl;
        cpptrace::from_current_exception().print();
    }
}
