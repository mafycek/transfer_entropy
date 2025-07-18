
#include <iostream>
#include <ranges>
#include <unordered_map>
#include <stacktrace>
#include <print>

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

#include "msgpack.hpp"

//#include "ql/quantlib.hpp"
//#include "ql/math/statistics/histogram.cpp"

#include "random_samples.h"
#include "renyi_entropy.h"
#include "utils.h"
#include "cpptrace_helper.h"

// required to stacktrace
typedef boost::error_info<struct tag_stacktrace, boost::stacktrace::stacktrace> traced;

template <class E>
void throw_with_trace(const E& e) {
    throw boost::enable_error_info(e)
        << traced(boost::stacktrace::stacktrace());
}

namespace po = boost::program_options;
namespace bio = boost::iostreams;

struct zstd_ostream : boost::iostreams::filtering_ostream
{
    zstd_ostream(std::ostream& os)
    {
        bio::zstd_params zstd_params ( 13 );
        bio::filtering_ostream::push( bio::zstd_compressor{zstd_params} );
        bio::filtering_ostream::push( os );
    }
};

int main ( int argc, char *argv[] )
{
    cpptrace::absorb_trace_exceptions(false);
    cpptrace::use_default_stderr_logger();
    cpptrace::register_terminate_handler();
    warmup_cpptrace();
    segfault_handler_cpptrace ();

    std::string directory, output;
    unsigned int maximal_neighborhood;
    unsigned int start_neighbor = 5; // initial neighbour for averaging over neighbours
    unsigned int runs, runs_surrogate, runs_shuffle;
    unsigned int number_of_threads;
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

        unsigned int postselection_X_future;
        unsigned int postselection_X_history;
        unsigned int postselection_Y_history;

        renyi_entropy::renyi_entropy<calculation_type>::collection_result_conditional_information_transfer_type collection_result_RTE;
        renyi_entropy::renyi_entropy<calculation_type>::type_average_result_conditional_information_transfer_type processing_RTE;

        for ( auto &item : alpha )
        {
            std::cout << item << " ";
        }
        std::cout << std::endl;
        Eigen::MatrixXd dataset;
        random_samples::samples_normal_distribution_uncorrelated ( dataset, 0, 1, 1000, 2 );

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
                    auto [dataset1, dataset2] =
                        renyi_entropy::renyi_entropy<calculation_type>::prepare_dataset (
                            dataset, swap_datasets, shuffle_indicator, surrogate_indicator, 1, 1 );

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

                                auto indices_to_use_range =
                                    std::ranges::iota_view{1U, maximal_neighborhood};
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

                                auto transfer_entropy_results = renyi_entropy::renyi_entropy<calculation_type>::
                                                        renyi_conditional_information_transfer (
                                                            y_future, y_history, z_history,
                                                            configuration_renyi_entropy );
                                collection_result_RTE[collection_key][std::make_tuple(swap_datasets, shuffle_indicator, surrogate_indicator, sample)] = transfer_entropy_results;
                            }
                        }
                    }
                    std::cout << std::format("Sw: {} Sh:{} Sur:{} R:{}", swap_datasets, shuffle_indicator, surrogate_indicator, sample) << std::endl;
                }
            }
        }

        auto collection_key_views = std::views::keys(collection_result_RTE);
        std::vector<renyi_entropy::renyi_entropy<calculation_type>::collection_conditional_information_transfer_key_type> collection_keys{ collection_key_views.begin(), collection_key_views.end() };
        for(auto & collection_key: collection_keys)
        {
            auto & result_conditional_information_transfer = collection_result_RTE[collection_key];

            // calculate directional RTE
            for (  auto swap_datasets :
                {
                    true, false
                })
            {
                auto &main_conditional_renyi_entropy_results = result_conditional_information_transfer[std::make_tuple(swap_datasets, false, false, 0)][renyi_entropy::conditional_renyi_entropy_label];
                std::set<unsigned int> neighbours;
                std::set<double> alphas;

                for (auto & type_of_collection: {std::tuple<bool,bool>({true, false}), std::tuple<bool,bool>({false, true}) })
                {
                    auto [shuffle_indicator, surrogate_indicator] = type_of_collection;

                    auto CRE_key_views = std::views::keys(main_conditional_renyi_entropy_results);
                    std::vector<renyi_entropy::renyi_entropy<calculation_type>::renyi_key_type> CRE_keys{ CRE_key_views.begin(), CRE_key_views.end() };
                    for (auto & CRE_key: CRE_keys)
                    {
                        auto [neighbor, alpha] = CRE_key;
                        neighbours.insert(neighbor);
                        alphas.insert(alpha);
                        // effective RTE statistical
                        //# boost::accumulators::stats<boost::accumulators::tag::count, boost::accumulators::tag::sum, boost::accumulators::tag::immediate_mean, boost::accumulators::tag::variance, boost::accumulators::tag::median, boost::accumulators::tag::p_square_quantile >
                        using stats_accumulators = boost::accumulators::stats<boost::accumulators::tag::mean, boost::accumulators::tag::variance, boost::accumulators::tag::median, boost::accumulators::tag::tail_quantile<boost::accumulators::left>>;
                        boost::accumulators::accumulator_set< double, stats_accumulators > accumulator_statistics(boost::accumulators::tag::tail<boost::accumulators::left>::cache_size = 100);

                        for ( unsigned int run = 0; run < runs; ++ run )
                        {
                            auto &counterpart_conditional_renyi_entropy_results = result_conditional_information_transfer[std::tuple(swap_datasets, shuffle_indicator, surrogate_indicator, run)][renyi_entropy::conditional_renyi_entropy_label];
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
                        auto final_CRE_key = std::make_tuple(false, swap_datasets, shuffle_indicator, surrogate_indicator, 0);
                        processing_RTE[renyi_entropy::average_runs]["mean"][collection_key][final_CRE_key][CRE_key] = mean;
                        processing_RTE[renyi_entropy::average_runs]["variance"][collection_key][final_CRE_key][CRE_key] = variance;
                        processing_RTE[renyi_entropy::average_runs]["median"][collection_key][final_CRE_key][CRE_key] = median;
                        processing_RTE[renyi_entropy::average_runs]["quantile_01"][collection_key][final_CRE_key][CRE_key] = quantile_01;
                        processing_RTE[renyi_entropy::average_runs]["quantile_02"][collection_key][final_CRE_key][CRE_key] = quantile_02;
                        processing_RTE[renyi_entropy::average_runs]["quantile_03"][collection_key][final_CRE_key][CRE_key] = quantile_03;
                        processing_RTE[renyi_entropy::average_runs]["quantile_04"][collection_key][final_CRE_key][CRE_key] = quantile_04;
                        processing_RTE[renyi_entropy::average_runs]["quantile_06"][collection_key][final_CRE_key][CRE_key] = quantile_06;
                        processing_RTE[renyi_entropy::average_runs]["quantile_07"][collection_key][final_CRE_key][CRE_key] = quantile_07;
                        processing_RTE[renyi_entropy::average_runs]["quantile_08"][collection_key][final_CRE_key][CRE_key] = quantile_08;
                        processing_RTE[renyi_entropy::average_runs]["quantile_09"][collection_key][final_CRE_key][CRE_key] = quantile_09;
                    }
                }

                for (auto & type_of_collection: {std::tuple<bool,bool>({true, false}), std::tuple<bool,bool>({false, true}) })
                {
                    auto [shuffle_indicator, surrogate_indicator] = type_of_collection;

                    for (auto alpha: alphas)
                    {
                        renyi_entropy::renyi_entropy<calculation_type>::renyi_key_type CRE_key(0, alpha);

                        // effective RTE
                        using stats_accumulators = boost::accumulators::stats<boost::accumulators::tag::mean, boost::accumulators::tag::variance, boost::accumulators::tag::median, boost::accumulators::tag::tail_quantile<boost::accumulators::left>>;
                        boost::accumulators::accumulator_set< double, stats_accumulators > accumulator_statistics(boost::accumulators::tag::tail<boost::accumulators::left>::cache_size = runs*(maximal_neighborhood- start_neighbor));

                        for ( unsigned int neighbour = start_neighbor; neighbour < maximal_neighborhood; ++neighbour)
                        {
                            for ( unsigned int run = 0; run < runs; ++ run )
                            {
                                renyi_entropy::renyi_entropy<calculation_type>::renyi_key_type data_CRE_key(neighbour, alpha);
                                auto &counterpart_conditional_renyi_entropy_results = result_conditional_information_transfer[std::tuple(swap_datasets, shuffle_indicator, surrogate_indicator, run)][renyi_entropy::conditional_renyi_entropy_label];
                                auto entropy_difference = main_conditional_renyi_entropy_results[data_CRE_key] - counterpart_conditional_renyi_entropy_results[data_CRE_key];
                                accumulator_statistics ( entropy_difference );
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
                        auto final_CRE_key = std::make_tuple(false, swap_datasets, shuffle_indicator, surrogate_indicator, 0);
                        processing_RTE[renyi_entropy::average_runs_neighbors]["mean"][collection_key][final_CRE_key][CRE_key] = mean;
                        processing_RTE[renyi_entropy::average_runs_neighbors]["variance"][collection_key][final_CRE_key][CRE_key] = variance;
                        processing_RTE[renyi_entropy::average_runs_neighbors]["median"][collection_key][final_CRE_key][CRE_key] = median;
                        processing_RTE[renyi_entropy::average_runs_neighbors]["quantile_01"][collection_key][final_CRE_key][CRE_key] = quantile_01;
                        processing_RTE[renyi_entropy::average_runs_neighbors]["quantile_02"][collection_key][final_CRE_key][CRE_key] = quantile_02;
                        processing_RTE[renyi_entropy::average_runs_neighbors]["quantile_03"][collection_key][final_CRE_key][CRE_key] = quantile_03;
                        processing_RTE[renyi_entropy::average_runs_neighbors]["quantile_04"][collection_key][final_CRE_key][CRE_key] = quantile_04;
                        processing_RTE[renyi_entropy::average_runs_neighbors]["quantile_06"][collection_key][final_CRE_key][CRE_key] = quantile_06;
                        processing_RTE[renyi_entropy::average_runs_neighbors]["quantile_07"][collection_key][final_CRE_key][CRE_key] = quantile_07;
                        processing_RTE[renyi_entropy::average_runs_neighbors]["quantile_08"][collection_key][final_CRE_key][CRE_key] = quantile_08;
                        processing_RTE[renyi_entropy::average_runs_neighbors]["quantile_09"][collection_key][final_CRE_key][CRE_key] = quantile_09;
                    }
                }
            }

            // ballance calculation
            auto swap_datasets = false;
            auto &main_conditional_renyi_entropy_results = result_conditional_information_transfer[std::make_tuple(swap_datasets, false, false, 0)][renyi_entropy::conditional_renyi_entropy_label];
            auto &main_reverse_conditional_renyi_entropy_results = result_conditional_information_transfer[std::make_tuple(!swap_datasets, false, false, 0)][renyi_entropy::conditional_renyi_entropy_label];
            std::set<unsigned int> neighbours;
            std::set<double> alphas;

            for (auto & type_of_collection: {std::tuple<bool,bool>({true, false}), std::tuple<bool,bool>({false, true}) })
            {
                auto [shuffle_indicator, surrogate_indicator] = type_of_collection;

                auto CRE_key_views = std::views::keys(main_conditional_renyi_entropy_results);
                std::vector<renyi_entropy::renyi_entropy<calculation_type>::renyi_key_type> CRE_keys{ CRE_key_views.begin(), CRE_key_views.end() };
                for (auto & CRE_key: CRE_keys)
                {
                    auto [neighbor, alpha] = CRE_key;
                    neighbours.insert(neighbor);
                    alphas.insert(alpha);
                    // effective RTE statistical
                    //# boost::accumulators::stats<boost::accumulators::tag::count, boost::accumulators::tag::sum, boost::accumulators::tag::immediate_mean, boost::accumulators::tag::variance, boost::accumulators::tag::median, boost::accumulators::tag::p_square_quantile >
                    using stats_accumulators = boost::accumulators::stats<boost::accumulators::tag::mean, boost::accumulators::tag::variance, boost::accumulators::tag::median, boost::accumulators::tag::tail_quantile<boost::accumulators::left>>;
                    boost::accumulators::accumulator_set< double, stats_accumulators > accumulator_statistics(boost::accumulators::tag::tail<boost::accumulators::left>::cache_size = 100);

                    for ( unsigned int run = 0; run < runs; ++ run )
                    {
                        auto &counterpart_conditional_renyi_entropy_results = result_conditional_information_transfer[std::tuple(swap_datasets, shuffle_indicator, surrogate_indicator, run)][renyi_entropy::conditional_renyi_entropy_label];
                        auto &reverse_counterpart_conditional_renyi_entropy_results = result_conditional_information_transfer[std::tuple(!swap_datasets, shuffle_indicator, surrogate_indicator, run)][renyi_entropy::conditional_renyi_entropy_label];

                        auto entropy_difference = main_conditional_renyi_entropy_results[CRE_key] - counterpart_conditional_renyi_entropy_results[CRE_key] - main_reverse_conditional_renyi_entropy_results[CRE_key] + reverse_counterpart_conditional_renyi_entropy_results[CRE_key];
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
                    auto final_CRE_key = std::make_tuple(true, swap_datasets, shuffle_indicator, surrogate_indicator, 0);
                    processing_RTE[renyi_entropy::average_runs]["mean"][collection_key][final_CRE_key][CRE_key] = mean;
                    processing_RTE[renyi_entropy::average_runs]["variance"][collection_key][final_CRE_key][CRE_key] = variance;
                    processing_RTE[renyi_entropy::average_runs]["median"][collection_key][final_CRE_key][CRE_key] = median;
                    processing_RTE[renyi_entropy::average_runs]["quantile_01"][collection_key][final_CRE_key][CRE_key] = quantile_01;
                    processing_RTE[renyi_entropy::average_runs]["quantile_02"][collection_key][final_CRE_key][CRE_key] = quantile_02;
                    processing_RTE[renyi_entropy::average_runs]["quantile_03"][collection_key][final_CRE_key][CRE_key] = quantile_03;
                    processing_RTE[renyi_entropy::average_runs]["quantile_04"][collection_key][final_CRE_key][CRE_key] = quantile_04;
                    processing_RTE[renyi_entropy::average_runs]["quantile_06"][collection_key][final_CRE_key][CRE_key] = quantile_06;
                    processing_RTE[renyi_entropy::average_runs]["quantile_07"][collection_key][final_CRE_key][CRE_key] = quantile_07;
                    processing_RTE[renyi_entropy::average_runs]["quantile_08"][collection_key][final_CRE_key][CRE_key] = quantile_08;
                    processing_RTE[renyi_entropy::average_runs]["quantile_09"][collection_key][final_CRE_key][CRE_key] = quantile_09;
                }
            }

            for (auto & type_of_collection: {std::tuple<bool,bool>({true, false}), std::tuple<bool,bool>({false, true}) })
            {
                auto [shuffle_indicator, surrogate_indicator] = type_of_collection;

                for (auto alpha: alphas)
                {
                    renyi_entropy::renyi_entropy<calculation_type>::renyi_key_type CRE_key(0, alpha);

                    // effective RTE
                    using stats_accumulators = boost::accumulators::stats<boost::accumulators::tag::mean, boost::accumulators::tag::variance, boost::accumulators::tag::median, boost::accumulators::tag::tail_quantile<boost::accumulators::left>>;
                    boost::accumulators::accumulator_set< double, stats_accumulators > accumulator_statistics(boost::accumulators::tag::tail<boost::accumulators::left>::cache_size = runs*(maximal_neighborhood- start_neighbor));

                    for ( unsigned int neighbour = start_neighbor; neighbour < maximal_neighborhood; ++neighbour)
                    {
                        for ( unsigned int run = 0; run < runs; ++ run )
                        {
                            renyi_entropy::renyi_entropy<calculation_type>::renyi_key_type data_CRE_key(neighbour, alpha);
                            auto &counterpart_conditional_renyi_entropy_results = result_conditional_information_transfer[std::tuple(swap_datasets, shuffle_indicator, surrogate_indicator, run)][renyi_entropy::conditional_renyi_entropy_label];
                            auto &reverse_counterpart_conditional_renyi_entropy_results = result_conditional_information_transfer[std::tuple(!swap_datasets, shuffle_indicator, surrogate_indicator, run)][renyi_entropy::conditional_renyi_entropy_label];

                            auto entropy_difference = main_conditional_renyi_entropy_results[data_CRE_key] - counterpart_conditional_renyi_entropy_results[data_CRE_key] - main_reverse_conditional_renyi_entropy_results[CRE_key] + reverse_counterpart_conditional_renyi_entropy_results[CRE_key];
                            accumulator_statistics ( entropy_difference );
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
                    auto final_CRE_key = std::make_tuple(true, swap_datasets, shuffle_indicator, surrogate_indicator, 0);
                    processing_RTE[renyi_entropy::average_runs_neighbors]["mean"][collection_key][final_CRE_key][CRE_key] = mean;
                    processing_RTE[renyi_entropy::average_runs_neighbors]["variance"][collection_key][final_CRE_key][CRE_key] = variance;
                    processing_RTE[renyi_entropy::average_runs_neighbors]["median"][collection_key][final_CRE_key][CRE_key] = median;
                    processing_RTE[renyi_entropy::average_runs_neighbors]["quantile_01"][collection_key][final_CRE_key][CRE_key] = quantile_01;
                    processing_RTE[renyi_entropy::average_runs_neighbors]["quantile_02"][collection_key][final_CRE_key][CRE_key] = quantile_02;
                    processing_RTE[renyi_entropy::average_runs_neighbors]["quantile_03"][collection_key][final_CRE_key][CRE_key] = quantile_03;
                    processing_RTE[renyi_entropy::average_runs_neighbors]["quantile_04"][collection_key][final_CRE_key][CRE_key] = quantile_04;
                    processing_RTE[renyi_entropy::average_runs_neighbors]["quantile_06"][collection_key][final_CRE_key][CRE_key] = quantile_06;
                    processing_RTE[renyi_entropy::average_runs_neighbors]["quantile_07"][collection_key][final_CRE_key][CRE_key] = quantile_07;
                    processing_RTE[renyi_entropy::average_runs_neighbors]["quantile_08"][collection_key][final_CRE_key][CRE_key] = quantile_08;
                    processing_RTE[renyi_entropy::average_runs_neighbors]["quantile_09"][collection_key][final_CRE_key][CRE_key] = quantile_09;
                }
            }
        }

        // save results
        std::stringstream ss;
        msgpack::pack(ss, processing_RTE); //renyi_entropy::renyi_entropy<calculation_type>::storage_RTE(collection_result_RTE, processing_RTE)
        std::cout << ss.str().size() << std::endl;
        boost::filesystem::path output_file =
        boost::filesystem::path(directory) / boost::filesystem::path(output);
        boost::filesystem::ofstream output_file_handler ( output_file );
        zstd_ostream zstd_compression_stream{output_file_handler};

        output_file_handler << ss.str();
    }
    CPPTRACE_CATCH (std::exception& exc)
    {
        std::cerr << "Exception: "<< exc.what() << std::endl;
        cpptrace::from_current_exception().print();
    }
}
