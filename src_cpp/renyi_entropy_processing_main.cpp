
#include <iostream>
#include <map>
#include <array>
#include <variant>
#include <ranges>
#include <sstream>
#include <vector>

#include <cpptrace/from_current.hpp>
#include <cpptrace/formatting.hpp>
#include <cpptrace/cpptrace.hpp>

#include <boost/program_options.hpp>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/zstd.hpp>
#include <boost/iostreams/filter/lzma.hpp>
#include <boost/iostreams/filter/bzip2.hpp>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

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
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/kurtosis.hpp>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <msgpack.hpp>

#include "cpptrace_helper.h"
#include "renyi_entropy.h"
#include "serialize_tuple.h"

namespace po = boost::program_options;
namespace bio = boost::iostreams;
typedef double calculation_type;

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

    po::options_description desc ( "Allowed options" );
    desc.add_options()
    ( "help,h", "produce help message" )
    ( "directory,d", po::value<std::string>()->default_value ( "." ), "Folder to export results" )
    ( "input,i", po::value<std::string>()->default_value ( "CRE.bin.zstd" ), "Raw results" )
    ( "output,o", po::value<std::string>()->default_value ( "CRE_processed.bin.zstd" ), "Processed results" );

    po::variables_map vm;
    po::store ( po::parse_command_line ( argc, argv, desc ), vm );
    po::notify ( vm );

    CPPTRACE_TRY
    {
        std::string input = vm["input"].as<std::string>();
        std::string output = vm["output"].as<std::string>();
        std::string directory = vm["directory"].as<std::string>();

        boost::filesystem::path input_filename =
        boost::filesystem::path ( directory ) / boost::filesystem::path ( input );
        boost::filesystem::ifstream input_file_handler ( input_filename );

        constexpr unsigned int maximal_neighborhood = 100;
        constexpr unsigned int runs = 11;
        constexpr unsigned int small_accumulator_cache_size = maximal_neighborhood;
        constexpr unsigned int big_accumulator_cache_size = maximal_neighborhood * runs;
        auto indices_to_use_range = std::ranges::iota_view{1U, maximal_neighborhood};
        auto indices_to_use_averaging = std::ranges::iota_view{10U, maximal_neighborhood};

        if ( input_file_handler.is_open() )
        {
            BOOST_LOG_TRIVIAL ( info ) << "Input file is opened: " << input_filename;

            auto start_processing =
                std::chrono::high_resolution_clock::now();

            bio::filtering_istream istream;
            istream.push ( bio::zstd_decompressor{} );
            istream.push ( input_file_handler );

            bool msg_pack = false;
            renyi_entropy::renyi_entropy<calculation_type>::result_RTE_t collection_result_RTEs;
            renyi_entropy::renyi_entropy<calculation_type>::type_average_result_conditional_information_transfer_type processing_RTE;
            if ( msg_pack == true )
            {
                std::string buffer ( std::istreambuf_iterator<char> {istream}, {} );
                istream >> buffer;
                msgpack::object_handle object_handle_buffer = msgpack::unpack ( buffer.data(), buffer.size() );

                // deserialized object is valid during the msgpack::object_handle instance is alive.
                msgpack::object deserialized = object_handle_buffer.get();

                deserialized.convert ( collection_result_RTEs );
            }
            else
            {
                boost::archive::binary_iarchive input_archive(istream);
                input_archive >> collection_result_RTEs;
            }

            auto start_result_processing = std::chrono::high_resolution_clock::now();
            std::set<std::tuple<bool, bool, bool, std::vector<unsigned int>, std::vector<unsigned int>, std::vector<unsigned int>>> main_RTE_titles;
            std::set<std::tuple<bool, bool, bool, std::vector<unsigned int>, std::vector<unsigned int>, std::vector<unsigned int>>> ballance_RTE_titles;
            std::set<std::tuple<bool, bool, bool, std::vector<unsigned int>, std::vector<unsigned int>, std::vector<unsigned int>>> ballance_effective_RTE_titles;
            std::set<std::tuple<bool, bool, bool, std::vector<unsigned int>, std::vector<unsigned int>, std::vector<unsigned int>>> subentropies_RTE_titles;

            {
                auto collection_key_views = std::views::keys ( collection_result_RTEs );
                std::vector<decltype ( collection_result_RTEs ) ::key_type> collection_keys{ collection_key_views.begin(), collection_key_views.end() };
                for ( auto & collection_key: collection_keys )
                {
                    const auto &[swap_datasets, shuffle_indicator, surrogate_indicator, sample, future_first, histories_first, histories_second, type_of_calculation, neighbor] = collection_key;
                    if ( type_of_calculation == renyi_entropy::conditional_renyi_entropy_label )
                    {
                        if ( ( ( shuffle_indicator == true ) && ( surrogate_indicator == false ) ) || ( ( shuffle_indicator == false ) && ( surrogate_indicator == true ) ) )
                        {
                            main_RTE_titles.insert ( std::make_tuple ( swap_datasets, shuffle_indicator, surrogate_indicator, future_first, histories_first, histories_second ) );

                        }
                        if ( ( ( shuffle_indicator == false ) && ( surrogate_indicator == false ) ) || ( ( shuffle_indicator == true ) && ( surrogate_indicator == false ) ) )
                        {
                            if ( swap_datasets == false )
                            {
                                ballance_RTE_titles.insert ( std::make_tuple ( swap_datasets, shuffle_indicator, surrogate_indicator, future_first, histories_first, histories_second ) );
                            }
                        }
                        if ( ( ( shuffle_indicator == true ) && ( surrogate_indicator == false ) ) || ( ( shuffle_indicator == false ) && ( surrogate_indicator == true ) ) )
                        {
                            if ( swap_datasets == false )
                            {
                                ballance_effective_RTE_titles.insert ( std::make_tuple ( swap_datasets, shuffle_indicator, surrogate_indicator, future_first, histories_first, histories_second ) );
                            }
                        }
                    }
                    else
                    {
                        if ( ( /* ( shuffle_indicator == false ) && */ ( surrogate_indicator == false ) ) )
                        {
                            subentropies_RTE_titles.insert ( std::make_tuple ( swap_datasets, shuffle_indicator, surrogate_indicator, future_first, histories_first, histories_second ) );
                            std::cout << std::format ( "{} {} {} {} {}", swap_datasets, shuffle_indicator, surrogate_indicator, sample, neighbor ) << std::endl;
                        }
                    }
                }
            }

            using stats_accumulators = boost::accumulators::stats<boost::accumulators::tag::count, boost::accumulators::tag::sum, boost::accumulators::tag::mean, boost::accumulators::tag::variance, boost::accumulators::tag::skewness, boost::accumulators::tag::kurtosis, boost::accumulators::tag::median, boost::accumulators::tag::tail_quantile<boost::accumulators::left>>;

            // subentropies RTE calculation, run statistics
            for ( auto & subentropies_RTE_title: subentropies_RTE_titles )
            {
                for ( const auto & type_of_calculation:
                        {
                            renyi_entropy::renyi_entropy_X_present_history_label, renyi_entropy::renyi_entropy_XY_history_label, renyi_entropy::joint_renyi_entropy_label, renyi_entropy::renyi_entropy_X_history_label
                        } )
                {
                    const auto [swap_datasets, shuffle_indicator, surrogate_indicator, future_first, histories_first, histories_second] = subentropies_RTE_title;

                    auto indices_to_use_range = std::ranges::iota_view{1U, maximal_neighborhood};
                    for ( const auto neighbor: indices_to_use_range )
                    {
                        const auto collection_key = std::make_tuple ( swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, type_of_calculation, neighbor );
                        if (!collection_result_RTEs.contains(collection_key))
                        {
                            throw collection_key;
                        }
                        auto & main_conditional_renyi_entropy_results = collection_result_RTEs[collection_key];

                        const auto CRE_key_views = std::views::keys ( main_conditional_renyi_entropy_results );
                        std::vector<renyi_entropy::renyi_entropy<double>::renyi_key_type> CRE_keys{ CRE_key_views.begin(), CRE_key_views.end() };
                        for ( auto & CRE_key: CRE_keys )
                        {
                            boost::accumulators::accumulator_set< double, stats_accumulators > accumulator_statistics ( boost::accumulators::tag::tail<boost::accumulators::left>::cache_size = small_accumulator_cache_size );

                            auto maximal_runs = ( ( shuffle_indicator||surrogate_indicator ) ? runs : 1 );
                            for ( unsigned int run = 0; run < maximal_runs ; ++ run )
                            {
                                const auto collection_key = std::make_tuple ( swap_datasets, shuffle_indicator, surrogate_indicator, run, future_first, histories_first, histories_second, type_of_calculation, neighbor );
                                auto & counterpart_conditional_renyi_entropy_results = collection_result_RTEs[collection_key];
                                if (!collection_result_RTEs.contains(collection_key))
                                {
                                    throw collection_key;
                                }
                                auto entropy_difference = counterpart_conditional_renyi_entropy_results[CRE_key];
                                accumulator_statistics ( entropy_difference );
                            }
                            auto mean = boost::accumulators::mean ( accumulator_statistics );
                            const auto final_CRE_key_mean = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "mean", renyi_entropy::average_runs );
                            processing_RTE[final_CRE_key_mean][CRE_key] = mean;

                            if (maximal_runs != 1)
                            {
                                auto variance = boost::accumulators::variance ( accumulator_statistics );
                                auto median = boost::accumulators::median ( accumulator_statistics );
                                auto quantile_01 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.1 );
                                auto quantile_02 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.2 );
                                auto quantile_03 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.3 );
                                auto quantile_04 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.4 );
                                auto quantile_06 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.6 );
                                auto quantile_07 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.7 );
                                auto quantile_08 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.8 );
                                auto quantile_09 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.9 );
                                const auto final_CRE_key_standard_deviation = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "standard_deviation", renyi_entropy::average_runs );
                                const auto final_CRE_key_median = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "median", renyi_entropy::average_runs );
                                const auto final_CRE_key_quantile_01 = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_01", renyi_entropy::average_runs );
                                const auto final_CRE_key_quantile_02 = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_02", renyi_entropy::average_runs );
                                const auto final_CRE_key_quantile_03 = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_03", renyi_entropy::average_runs );
                                const auto final_CRE_key_quantile_04 = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_04", renyi_entropy::average_runs );
                                const auto final_CRE_key_quantile_06 = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_06", renyi_entropy::average_runs );
                                const auto final_CRE_key_quantile_07 = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_07", renyi_entropy::average_runs );
                                const auto final_CRE_key_quantile_08 = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_08", renyi_entropy::average_runs );
                                const auto final_CRE_key_quantile_09 = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_09", renyi_entropy::average_runs );

                                processing_RTE[final_CRE_key_standard_deviation][CRE_key] = sqrt ( variance );
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
                }
            }

            // effective RTE calculation, run statistics
            for ( auto & main_RTE_title: main_RTE_titles )
            {
                const auto [swap_datasets, shuffle_indicator, surrogate_indicator, future_first, histories_first, histories_second] = main_RTE_title;
                const auto type_of_calculation = renyi_entropy::conditional_renyi_entropy_label;
                for ( const auto neighbor: indices_to_use_range )
                {
                    const auto collection_key = std::make_tuple ( swap_datasets, false, false, 0, future_first, histories_first, histories_second, type_of_calculation, neighbor );
                    if (!collection_result_RTEs.contains(collection_key))
                    {
                        throw collection_key;
                    }
                    auto & main_conditional_renyi_entropy_results = collection_result_RTEs[collection_key];

                    const auto CRE_key_views = std::views::keys ( main_conditional_renyi_entropy_results );
                    std::vector<renyi_entropy::renyi_entropy<double>::renyi_key_type> CRE_keys{ CRE_key_views.begin(), CRE_key_views.end() };
                    for ( auto & CRE_key: CRE_keys )
                    {
                        boost::accumulators::accumulator_set< double, stats_accumulators > accumulator_statistics ( boost::accumulators::tag::tail<boost::accumulators::left>::cache_size = small_accumulator_cache_size );

                        auto maximal_runs = ( ( shuffle_indicator||surrogate_indicator ) ? runs : 1 );
                        for ( unsigned int run = 0; run < maximal_runs ; ++ run )
                        {
                            const auto collection_key = std::make_tuple ( swap_datasets, shuffle_indicator, surrogate_indicator, run, future_first, histories_first, histories_second, type_of_calculation, neighbor );
                            if (!collection_result_RTEs.contains(collection_key))
                            {
                                throw collection_key;
                            }
                            auto & counterpart_conditional_renyi_entropy_results = collection_result_RTEs[collection_key];

                            auto entropy_difference = main_conditional_renyi_entropy_results[CRE_key] - counterpart_conditional_renyi_entropy_results[CRE_key];
                            accumulator_statistics ( entropy_difference );
                        }
                        auto mean = boost::accumulators::mean ( accumulator_statistics );
                        auto variance = boost::accumulators::variance ( accumulator_statistics );
                        auto median = boost::accumulators::median ( accumulator_statistics );
                        auto quantile_01 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.1 );
                        auto quantile_02 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.2 );
                        auto quantile_03 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.3 );
                        auto quantile_04 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.4 );
                        auto quantile_06 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.6 );
                        auto quantile_07 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.7 );
                        auto quantile_08 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.8 );
                        auto quantile_09 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.9 );
                        const auto final_CRE_key_mean = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "mean", renyi_entropy::average_runs );
                        const auto final_CRE_key_standard_deviation = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "standard_deviation", renyi_entropy::average_runs );
                        const auto final_CRE_key_median = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "median", renyi_entropy::average_runs );
                        const auto final_CRE_key_quantile_01 = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_01", renyi_entropy::average_runs );
                        const auto final_CRE_key_quantile_02 = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_02", renyi_entropy::average_runs );
                        const auto final_CRE_key_quantile_03 = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_03", renyi_entropy::average_runs );
                        const auto final_CRE_key_quantile_04 = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_04", renyi_entropy::average_runs );
                        const auto final_CRE_key_quantile_06 = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_06", renyi_entropy::average_runs );
                        const auto final_CRE_key_quantile_07 = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_07", renyi_entropy::average_runs );
                        const auto final_CRE_key_quantile_08 = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_08", renyi_entropy::average_runs );
                        const auto final_CRE_key_quantile_09 = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_09", renyi_entropy::average_runs );

                        processing_RTE[final_CRE_key_mean][CRE_key] = mean;
                        processing_RTE[final_CRE_key_standard_deviation][CRE_key] = sqrt ( variance );
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
            for ( auto & ballance_RTE_title: ballance_RTE_titles )
            {
                const auto [swap_datasets, shuffle_indicator, surrogate_indicator, future_first, histories_first, histories_second] = ballance_RTE_title;
                const auto type_of_calculation = renyi_entropy::conditional_renyi_entropy_label;

                for ( const auto neighbor: indices_to_use_range )
                {
                    const auto collection_key = std::make_tuple ( swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, type_of_calculation, neighbor );
                    if (!collection_result_RTEs.contains(collection_key))
                    {
                        throw collection_key;
                    }
                    auto & main_conditional_renyi_entropy_results = collection_result_RTEs[collection_key];

                    const auto CRE_key_views = std::views::keys ( main_conditional_renyi_entropy_results );
                    std::vector<renyi_entropy::renyi_entropy<double>::renyi_key_type> CRE_keys{ CRE_key_views.begin(), CRE_key_views.end() };
                    for ( auto & CRE_key: CRE_keys )
                    {
                        boost::accumulators::accumulator_set< double, stats_accumulators > accumulator_statistics ( boost::accumulators::tag::tail<boost::accumulators::left>::cache_size = small_accumulator_cache_size );

                        for ( unsigned int run = 0; run < ( ( shuffle_indicator||surrogate_indicator ) ? runs : 1 ) ; ++ run )
                        {
                            const auto collection_key = std::make_tuple ( swap_datasets, shuffle_indicator, surrogate_indicator, run, future_first, histories_first, histories_second, type_of_calculation, neighbor );
                            if (!collection_result_RTEs.contains(collection_key))
                            {
                                throw collection_key;
                            }
                            auto & main_conditional_renyi_entropy_results = collection_result_RTEs[collection_key];

                            const auto collection_reverse_direction_key = std::make_tuple ( !swap_datasets, shuffle_indicator, surrogate_indicator, run, future_first, histories_first, histories_second, type_of_calculation, neighbor );
                            auto & main_reverse_conditional_renyi_entropy_results = collection_result_RTEs[collection_reverse_direction_key];

                            auto entropy_difference = main_conditional_renyi_entropy_results[CRE_key] - main_reverse_conditional_renyi_entropy_results[CRE_key];
                            accumulator_statistics ( entropy_difference );
                        }

                        auto mean = boost::accumulators::mean ( accumulator_statistics );
                        auto variance = boost::accumulators::variance ( accumulator_statistics );
                        auto median = boost::accumulators::median ( accumulator_statistics );
                        auto quantile_01 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.1 );
                        auto quantile_02 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.2 );
                        auto quantile_03 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.3 );
                        auto quantile_04 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.4 );
                        auto quantile_06 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.6 );
                        auto quantile_07 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.7 );
                        auto quantile_08 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.8 );
                        auto quantile_09 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.9 );
                        const auto final_CRE_key_mean = std::make_tuple ( type_of_calculation, true, false, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "mean", renyi_entropy::average_runs );
                        const auto final_CRE_key_standard_deviation = std::make_tuple ( type_of_calculation, true, false, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "standard_deviation", renyi_entropy::average_runs );
                        const auto final_CRE_key_median = std::make_tuple ( type_of_calculation, true, false, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "median", renyi_entropy::average_runs );
                        const auto final_CRE_key_quantile_01 = std::make_tuple ( type_of_calculation, true, false, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_01", renyi_entropy::average_runs );
                        const auto final_CRE_key_quantile_02 = std::make_tuple ( type_of_calculation, true, false, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_02", renyi_entropy::average_runs );
                        const auto final_CRE_key_quantile_03 = std::make_tuple ( type_of_calculation, true, false, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_03", renyi_entropy::average_runs );
                        const auto final_CRE_key_quantile_04 = std::make_tuple ( type_of_calculation, true, false, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_04", renyi_entropy::average_runs );
                        const auto final_CRE_key_quantile_06 = std::make_tuple ( type_of_calculation, true, false, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_06", renyi_entropy::average_runs );
                        const auto final_CRE_key_quantile_07 = std::make_tuple ( type_of_calculation, true, false, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_07", renyi_entropy::average_runs );
                        const auto final_CRE_key_quantile_08 = std::make_tuple ( type_of_calculation, true, false, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_08", renyi_entropy::average_runs );
                        const auto final_CRE_key_quantile_09 = std::make_tuple ( type_of_calculation, true, false, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_09", renyi_entropy::average_runs );

                        processing_RTE[final_CRE_key_mean][CRE_key] = mean;
                        processing_RTE[final_CRE_key_standard_deviation][CRE_key] = sqrt ( variance );
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
            for ( auto & ballance_effective_RTE_title: ballance_effective_RTE_titles )
            {
                const auto [swap_datasets, shuffle_indicator, surrogate_indicator, future_first, histories_first, histories_second] = ballance_effective_RTE_title;
                const auto type_of_calculation = renyi_entropy::conditional_renyi_entropy_label;

                for ( const auto neighbor: indices_to_use_range )
                {
                    const auto collection_key = std::make_tuple ( swap_datasets, false, false, 0, future_first, histories_first, histories_second, type_of_calculation, neighbor );
                    if (!collection_result_RTEs.contains(collection_key))
                    {
                        throw collection_key;
                    }
                    auto & main_conditional_renyi_entropy_results = collection_result_RTEs[collection_key];
                    const auto swap_collection_key = std::make_tuple ( !swap_datasets, false, false, 0, future_first, histories_first, histories_second, type_of_calculation, neighbor );
                    if (!collection_result_RTEs.contains(swap_collection_key))
                    {
                        throw collection_key;
                    }
                    auto & main_reverse_conditional_renyi_entropy_results = collection_result_RTEs[swap_collection_key];

                    const auto CRE_key_views = std::views::keys ( main_conditional_renyi_entropy_results );
                    std::vector<renyi_entropy::renyi_entropy<double>::renyi_key_type> CRE_keys{ CRE_key_views.begin(), CRE_key_views.end() };
                    for ( auto & CRE_key: CRE_keys )
                    {
                        boost::accumulators::accumulator_set< double, stats_accumulators > accumulator_statistics ( boost::accumulators::tag::tail<boost::accumulators::left>::cache_size = small_accumulator_cache_size );

                        for ( unsigned int run = 0; run < runs; ++ run )
                        {
                            const auto collection_key = std::make_tuple ( swap_datasets, shuffle_indicator, surrogate_indicator, run, future_first, histories_first, histories_second, type_of_calculation, neighbor );
                            if (!collection_result_RTEs.contains(collection_key))
                            {
                                throw collection_key;
                            }
                            auto & counterpart_conditional_renyi_entropy_results = collection_result_RTEs[collection_key];

                            const auto collection_reverse_direction_key = std::make_tuple ( !swap_datasets, shuffle_indicator, surrogate_indicator, run, future_first, histories_first, histories_second, type_of_calculation, neighbor );
                            if (!collection_result_RTEs.contains(collection_reverse_direction_key))
                            {
                                throw collection_key;
                            }
                            auto & counterpart_reverse_conditional_renyi_entropy_results = collection_result_RTEs[collection_reverse_direction_key];

                            auto entropy_difference = main_conditional_renyi_entropy_results[CRE_key] - counterpart_conditional_renyi_entropy_results[CRE_key] - main_reverse_conditional_renyi_entropy_results[CRE_key] + counterpart_reverse_conditional_renyi_entropy_results[CRE_key];

                            accumulator_statistics ( entropy_difference );
                        }

                        auto mean = boost::accumulators::mean ( accumulator_statistics );
                        auto variance = boost::accumulators::variance ( accumulator_statistics );
                        auto median = boost::accumulators::median ( accumulator_statistics );
                        auto quantile_01 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.1 );
                        auto quantile_02 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.2 );
                        auto quantile_03 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.3 );
                        auto quantile_04 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.4 );
                        auto quantile_06 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.6 );
                        auto quantile_07 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.7 );
                        auto quantile_08 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.8 );
                        auto quantile_09 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.9 );
                        const auto final_CRE_key_mean = std::make_tuple ( type_of_calculation, true, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "mean", renyi_entropy::average_runs );
                        const auto final_CRE_key_standard_deviation = std::make_tuple ( type_of_calculation, true, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "standard_deviation", renyi_entropy::average_runs );
                        const auto final_CRE_key_median = std::make_tuple ( type_of_calculation, true, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "median", renyi_entropy::average_runs );
                        const auto final_CRE_key_quantile_01 = std::make_tuple ( type_of_calculation, true, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_01", renyi_entropy::average_runs );
                        const auto final_CRE_key_quantile_02 = std::make_tuple ( type_of_calculation, true, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_02", renyi_entropy::average_runs );
                        const auto final_CRE_key_quantile_03 = std::make_tuple ( type_of_calculation, true, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_03", renyi_entropy::average_runs );
                        const auto final_CRE_key_quantile_04 = std::make_tuple ( type_of_calculation, true, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_04", renyi_entropy::average_runs );
                        const auto final_CRE_key_quantile_06 = std::make_tuple ( type_of_calculation, true, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_06", renyi_entropy::average_runs );
                        const auto final_CRE_key_quantile_07 = std::make_tuple ( type_of_calculation, true, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_07", renyi_entropy::average_runs );
                        const auto final_CRE_key_quantile_08 = std::make_tuple ( type_of_calculation, true, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_08", renyi_entropy::average_runs );
                        const auto final_CRE_key_quantile_09 = std::make_tuple ( type_of_calculation, true, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_09", renyi_entropy::average_runs );

                        processing_RTE[final_CRE_key_mean][CRE_key] = mean;
                        processing_RTE[final_CRE_key_standard_deviation][CRE_key] = sqrt ( variance );
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

            // effective RTE calculation, run statistics
            for ( auto & main_RTE_title: main_RTE_titles )
            {
                const auto [swap_datasets, shuffle_indicator, surrogate_indicator, future_first, histories_first, histories_second] = main_RTE_title;
                const auto type_of_calculation = renyi_entropy::conditional_renyi_entropy_label;

                for ( const auto neighbor: indices_to_use_range )
                {
                    const auto collection_key = std::make_tuple ( swap_datasets, false, false, 0, future_first, histories_first, histories_second, type_of_calculation, indices_to_use_range[0] );
                    auto & main_conditional_renyi_entropy_results = collection_result_RTEs[collection_key];
                    const auto CRE_key_views = std::views::keys ( main_conditional_renyi_entropy_results );
                    std::vector<renyi_entropy::renyi_entropy<double>::renyi_key_type> CRE_keys{ CRE_key_views.begin(), CRE_key_views.end() };
                    for ( auto & CRE_key: CRE_keys )
                    {
                        boost::accumulators::accumulator_set< double, stats_accumulators > accumulator_statistics ( boost::accumulators::tag::tail<boost::accumulators::left>::cache_size = small_accumulator_cache_size );

                        for ( unsigned int run = 0; run < runs; ++ run )
                        {
                            const auto collection_key = std::make_tuple ( swap_datasets, shuffle_indicator, surrogate_indicator, run, future_first, histories_first, histories_second, type_of_calculation, neighbor );
                            auto & counterpart_conditional_renyi_entropy_results = collection_result_RTEs[collection_key];

                            auto entropy_difference = main_conditional_renyi_entropy_results[CRE_key] - counterpart_conditional_renyi_entropy_results[CRE_key];
                            accumulator_statistics ( entropy_difference );
                        }

                        auto mean = boost::accumulators::mean ( accumulator_statistics );
                        auto variance = boost::accumulators::variance ( accumulator_statistics );
                        auto median = boost::accumulators::median ( accumulator_statistics );
                        auto quantile_01 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.1 );
                        auto quantile_02 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.2 );
                        auto quantile_03 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.3 );
                        auto quantile_04 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.4 );
                        auto quantile_06 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.6 );
                        auto quantile_07 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.7 );
                        auto quantile_08 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.8 );
                        auto quantile_09 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.9 );
                        const auto final_CRE_key_mean = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "mean", renyi_entropy::average_runs );
                        const auto final_CRE_key_standard_deviation = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "standard_deviation", renyi_entropy::average_runs );
                        const auto final_CRE_key_median = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "median", renyi_entropy::average_runs );
                        const auto final_CRE_key_quantile_01 = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_01", renyi_entropy::average_runs );
                        const auto final_CRE_key_quantile_02 = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_02", renyi_entropy::average_runs );
                        const auto final_CRE_key_quantile_03 = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_03", renyi_entropy::average_runs );
                        const auto final_CRE_key_quantile_04 = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_04", renyi_entropy::average_runs );
                        const auto final_CRE_key_quantile_06 = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_06", renyi_entropy::average_runs );
                        const auto final_CRE_key_quantile_07 = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_07", renyi_entropy::average_runs );
                        const auto final_CRE_key_quantile_08 = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_08", renyi_entropy::average_runs );
                        const auto final_CRE_key_quantile_09 = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, neighbor, future_first, histories_first, histories_second, "quantile_09", renyi_entropy::average_runs );

                        processing_RTE[final_CRE_key_mean][CRE_key] = mean;
                        processing_RTE[final_CRE_key_standard_deviation][CRE_key] = std::sqrt ( variance );
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
            for ( auto & main_RTE_title: main_RTE_titles )
            {
                const auto [swap_datasets, shuffle_indicator, surrogate_indicator, future_first, histories_first, histories_second] = main_RTE_title;
                const auto type_of_calculation = renyi_entropy::conditional_renyi_entropy_label;

                const auto collection_key = std::make_tuple ( swap_datasets, false, false, 0, future_first, histories_first, histories_second, type_of_calculation, 1 );
                auto & main_conditional_renyi_entropy_results = collection_result_RTEs[collection_key];
                const auto CRE_key_views = std::views::keys ( main_conditional_renyi_entropy_results );
                std::vector<renyi_entropy::renyi_entropy<double>::renyi_key_type> CRE_keys{ CRE_key_views.begin(), CRE_key_views.end() };

                for ( auto & CRE_key: CRE_keys )
                {
                    boost::accumulators::accumulator_set< double, stats_accumulators >
                    // cache needs to contain complete dataset = #runs * #neighbors
                    accumulator_statistics ( boost::accumulators::tag::tail<boost::accumulators::left>::cache_size = big_accumulator_cache_size );

                    int i = 0;
                    for ( const auto neighbor: indices_to_use_averaging )
                    {
                        const auto collection_key = std::make_tuple ( swap_datasets, false, false, 0, future_first, histories_first, histories_second, type_of_calculation, neighbor );
                        auto & main_conditional_renyi_entropy_results = collection_result_RTEs[collection_key];

                        for ( unsigned int run = 0; run < runs; ++ run )
                        {
                            const auto collection_key = std::make_tuple ( swap_datasets, shuffle_indicator, surrogate_indicator, run, future_first, histories_first, histories_second, type_of_calculation, neighbor );
                            auto & counterpart_conditional_renyi_entropy_results = collection_result_RTEs[collection_key];

                            auto entropy_difference = main_conditional_renyi_entropy_results[CRE_key] - counterpart_conditional_renyi_entropy_results[CRE_key];
                            accumulator_statistics ( entropy_difference );
                            ++ i;
                        }
                    }
                    auto mean = boost::accumulators::mean ( accumulator_statistics );
                    auto variance = boost::accumulators::variance ( accumulator_statistics );
                    auto median = boost::accumulators::median ( accumulator_statistics );
                    auto quantile_01 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.1 );
                    auto quantile_02 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.2 );
                    auto quantile_03 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.3 );
                    auto quantile_04 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.4 );
                    auto quantile_06 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.6 );
                    auto quantile_07 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.7 );
                    auto quantile_08 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.8 );
                    auto quantile_09 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.9 );
                    const auto final_CRE_key_mean = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "mean", renyi_entropy::average_runs_neighbors );
                    const auto final_CRE_key_standard_deviation = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "standard_deviation", renyi_entropy::average_runs_neighbors );
                    const auto final_CRE_key_median = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "median", renyi_entropy::average_runs_neighbors );
                    const auto final_CRE_key_quantile_01 = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "quantile_01", renyi_entropy::average_runs_neighbors );
                    const auto final_CRE_key_quantile_02 = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "quantile_02", renyi_entropy::average_runs_neighbors );
                    const auto final_CRE_key_quantile_03 = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "quantile_03", renyi_entropy::average_runs_neighbors );
                    const auto final_CRE_key_quantile_04 = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "quantile_04", renyi_entropy::average_runs_neighbors );
                    const auto final_CRE_key_quantile_06 = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "quantile_06", renyi_entropy::average_runs_neighbors );
                    const auto final_CRE_key_quantile_07 = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "quantile_07", renyi_entropy::average_runs_neighbors );
                    const auto final_CRE_key_quantile_08 = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "quantile_08", renyi_entropy::average_runs_neighbors );
                    const auto final_CRE_key_quantile_09 = std::make_tuple ( type_of_calculation, false, true, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "quantile_09", renyi_entropy::average_runs_neighbors );

                    processing_RTE[final_CRE_key_mean][CRE_key] = mean;
                    processing_RTE[final_CRE_key_standard_deviation][CRE_key] = std::sqrt ( variance );
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

            // ballance effective calculation, run and neighbor statistics
            for ( auto & ballance_effective_RTE_title: ballance_effective_RTE_titles )
            {
                const auto [swap_datasets, shuffle_indicator, surrogate_indicator, future_first, histories_first, histories_second] = ballance_effective_RTE_title;
                const auto type_of_calculation = renyi_entropy::conditional_renyi_entropy_label;

                const auto collection_key = std::make_tuple ( swap_datasets, false, false, 0, future_first, histories_first, histories_second, type_of_calculation, 1 );
                auto & main_conditional_renyi_entropy_results = collection_result_RTEs[collection_key];
                const auto CRE_key_views = std::views::keys ( main_conditional_renyi_entropy_results );
                std::vector<renyi_entropy::renyi_entropy<double>::renyi_key_type> CRE_keys{ CRE_key_views.begin(), CRE_key_views.end() };

                for ( auto & CRE_key: CRE_keys )
                {
                    boost::accumulators::accumulator_set< double, stats_accumulators > accumulator_statistics ( boost::accumulators::tag::tail<boost::accumulators::left>::cache_size = big_accumulator_cache_size );

                    for ( const auto neighbor: indices_to_use_averaging )
                    {
                        const auto collection_key = std::make_tuple ( swap_datasets, false, false, 0, future_first, histories_first, histories_second, type_of_calculation, neighbor );
                        auto & main_conditional_renyi_entropy_results = collection_result_RTEs[collection_key];
                        const auto swap_collection_key = std::make_tuple ( !swap_datasets, false, false, 0, future_first, histories_first, histories_second, type_of_calculation, neighbor );
                        auto & main_reverse_conditional_renyi_entropy_results = collection_result_RTEs[swap_collection_key];

                        const auto CRE_key_views = std::views::keys ( main_conditional_renyi_entropy_results );
                        std::vector<renyi_entropy::renyi_entropy<double>::renyi_key_type> CRE_keys{ CRE_key_views.begin(), CRE_key_views.end() };

                        for ( unsigned int run = 0; run < runs; ++ run )
                        {
                            const auto collection_direction_key = std::make_tuple ( swap_datasets, shuffle_indicator, surrogate_indicator, run, future_first, histories_first, histories_second, type_of_calculation, neighbor );
                            auto & counterpart_conditional_renyi_entropy_results = collection_result_RTEs[collection_direction_key];
                            const auto collection_reverse_direction_key = std::make_tuple ( !swap_datasets, shuffle_indicator, surrogate_indicator, run, future_first, histories_first, histories_second, type_of_calculation, neighbor );
                            auto & counterpart_reverse_conditional_renyi_entropy_results = collection_result_RTEs[collection_reverse_direction_key];

                            auto entropy_difference = main_conditional_renyi_entropy_results[CRE_key] - counterpart_conditional_renyi_entropy_results[CRE_key] - main_reverse_conditional_renyi_entropy_results[CRE_key] + counterpart_reverse_conditional_renyi_entropy_results[CRE_key];
                            accumulator_statistics ( entropy_difference );
                        }
                    }

                    auto mean = boost::accumulators::mean ( accumulator_statistics );
                    auto variance = boost::accumulators::variance ( accumulator_statistics );
                    auto median = boost::accumulators::median ( accumulator_statistics );
                    auto quantile_01 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.1 );
                    auto quantile_02 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.2 );
                    auto quantile_03 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.3 );
                    auto quantile_04 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.4 );
                    auto quantile_06 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.6 );
                    auto quantile_07 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.7 );
                    auto quantile_08 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.8 );
                    auto quantile_09 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.9 );
                    const auto final_CRE_key_mean = std::make_tuple ( type_of_calculation, true, true, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "mean", renyi_entropy::average_runs_neighbors );
                    const auto final_CRE_key_standard_deviation = std::make_tuple ( type_of_calculation, true, true, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "standard_deviation", renyi_entropy::average_runs_neighbors );
                    const auto final_CRE_key_median = std::make_tuple ( type_of_calculation, true, true, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "median", renyi_entropy::average_runs_neighbors );
                    const auto final_CRE_key_quantile_01 = std::make_tuple ( type_of_calculation, true, true, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "quantile_01", renyi_entropy::average_runs_neighbors );
                    const auto final_CRE_key_quantile_02 = std::make_tuple ( type_of_calculation, true, true, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "quantile_02", renyi_entropy::average_runs_neighbors );
                    const auto final_CRE_key_quantile_03 = std::make_tuple ( type_of_calculation, true, true, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "quantile_03", renyi_entropy::average_runs_neighbors );
                    const auto final_CRE_key_quantile_04 = std::make_tuple ( type_of_calculation, true, true, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "quantile_04", renyi_entropy::average_runs_neighbors );
                    const auto final_CRE_key_quantile_06 = std::make_tuple ( type_of_calculation, true, true, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "quantile_06", renyi_entropy::average_runs_neighbors );
                    const auto final_CRE_key_quantile_07 = std::make_tuple ( type_of_calculation, true, true, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "quantile_07", renyi_entropy::average_runs_neighbors );
                    const auto final_CRE_key_quantile_08 = std::make_tuple ( type_of_calculation, true, true, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "quantile_08", renyi_entropy::average_runs_neighbors );
                    const auto final_CRE_key_quantile_09 = std::make_tuple ( type_of_calculation, true, true, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "quantile_09", renyi_entropy::average_runs_neighbors );

                    processing_RTE[final_CRE_key_mean][CRE_key] = mean;
                    processing_RTE[final_CRE_key_standard_deviation][CRE_key] = sqrt ( variance );
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

            // ballance calculation, run and neighbor statistics
            for ( auto & ballance_RTE_title: ballance_RTE_titles )
            {
                const auto [swap_datasets, shuffle_indicator, surrogate_indicator, future_first, histories_first, histories_second] = ballance_RTE_title;
                const auto type_of_calculation = renyi_entropy::conditional_renyi_entropy_label;

                const auto collection_key = std::make_tuple ( swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, type_of_calculation, 1 );
                auto & main_conditional_renyi_entropy_results = collection_result_RTEs[collection_key];

                const auto CRE_key_views = std::views::keys ( main_conditional_renyi_entropy_results );
                std::vector<renyi_entropy::renyi_entropy<double>::renyi_key_type> CRE_keys{ CRE_key_views.begin(), CRE_key_views.end() };
                for ( auto & CRE_key: CRE_keys )
                {
                    boost::accumulators::accumulator_set< double, stats_accumulators > accumulator_statistics ( boost::accumulators::tag::tail<boost::accumulators::left>::cache_size = big_accumulator_cache_size );

                    for ( const auto neighbor: indices_to_use_averaging )
                    {
                        for ( unsigned int run = 0; run < ( shuffle_indicator?runs:1 ) ; ++ run )
                        {
                            const auto collection_key = std::make_tuple ( swap_datasets, shuffle_indicator, surrogate_indicator, run, future_first, histories_first, histories_second, type_of_calculation, neighbor );
                            auto & main_conditional_renyi_entropy_results = collection_result_RTEs[collection_key];

                            const auto collection_reverse_direction_key = std::make_tuple ( !swap_datasets, shuffle_indicator, surrogate_indicator, run, future_first, histories_first, histories_second, type_of_calculation, neighbor );
                            auto & main_reverse_conditional_renyi_entropy_results = collection_result_RTEs[collection_reverse_direction_key];

                            auto entropy_difference = main_conditional_renyi_entropy_results[CRE_key] - main_reverse_conditional_renyi_entropy_results[CRE_key];
                            accumulator_statistics ( entropy_difference );
                        }
                    }

                    auto mean = boost::accumulators::mean ( accumulator_statistics );
                    auto variance = boost::accumulators::variance ( accumulator_statistics );
                    auto median = boost::accumulators::median ( accumulator_statistics );
                    auto quantile_01 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.1 );
                    auto quantile_02 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.2 );
                    auto quantile_03 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.3 );
                    auto quantile_04 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.4 );
                    auto quantile_06 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.6 );
                    auto quantile_07 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.7 );
                    auto quantile_08 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.8 );
                    auto quantile_09 = boost::accumulators::quantile ( accumulator_statistics, boost::accumulators::quantile_probability = 0.9 );
                    const auto final_CRE_key_mean = std::make_tuple ( type_of_calculation, true, false, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "mean", renyi_entropy::average_runs_neighbors );
                    const auto final_CRE_key_standard_deviation = std::make_tuple ( type_of_calculation, true, false, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "standard_deviation", renyi_entropy::average_runs_neighbors );
                    const auto final_CRE_key_median = std::make_tuple ( type_of_calculation, true, false, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "median", renyi_entropy::average_runs_neighbors );
                    const auto final_CRE_key_quantile_01 = std::make_tuple ( type_of_calculation, true, false, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "quantile_01", renyi_entropy::average_runs_neighbors );
                    const auto final_CRE_key_quantile_02 = std::make_tuple ( type_of_calculation, true, false, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "quantile_02", renyi_entropy::average_runs_neighbors );
                    const auto final_CRE_key_quantile_03 = std::make_tuple ( type_of_calculation, true, false, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "quantile_03", renyi_entropy::average_runs_neighbors );
                    const auto final_CRE_key_quantile_04 = std::make_tuple ( type_of_calculation, true, false, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "quantile_04", renyi_entropy::average_runs_neighbors );
                    const auto final_CRE_key_quantile_06 = std::make_tuple ( type_of_calculation, true, false, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "quantile_06", renyi_entropy::average_runs_neighbors );
                    const auto final_CRE_key_quantile_07 = std::make_tuple ( type_of_calculation, true, false, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "quantile_07", renyi_entropy::average_runs_neighbors );
                    const auto final_CRE_key_quantile_08 = std::make_tuple ( type_of_calculation, true, false, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "quantile_08", renyi_entropy::average_runs_neighbors );
                    const auto final_CRE_key_quantile_09 = std::make_tuple ( type_of_calculation, true, false, swap_datasets, shuffle_indicator, surrogate_indicator, 0, future_first, histories_first, histories_second, "quantile_09", renyi_entropy::average_runs_neighbors );

                    processing_RTE[final_CRE_key_mean][CRE_key] = mean;
                    processing_RTE[final_CRE_key_standard_deviation][CRE_key] = sqrt ( variance );
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

            // save results
            std::stringstream ss;
            msgpack::pack ( ss, processing_RTE ); //renyi_entropy::renyi_entropy<calculation_type>::storage_RTE(collection_result_RTE, processing_RTE)
            std::cout << "Result size: " << ss.str().size() << std::endl;
            boost::filesystem::path output_file =
                boost::filesystem::path ( directory ) / boost::filesystem::path ( output );
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

            auto end_processing =
                std::chrono::high_resolution_clock::now();
            auto elapsed_result_processing =
                end_processing - start_processing;
            auto microseconds_result_processing =
                std::chrono::duration_cast<std::chrono::microseconds> ( elapsed_result_processing ).count();
            BOOST_LOG_TRIVIAL ( trace ) << "Total time processing " << microseconds_result_processing / microseconds_in_second << " seconds";
        }
        else
        {
            BOOST_LOG_TRIVIAL ( error ) << "Input file cannot be opened: " << input_filename;
        }
    }
    CPPTRACE_CATCH ( renyi_entropy::renyi_entropy<calculation_type>::result_RTE_key& exc )
    {
        std::cerr << std::endl;
    }
//    CPPTRACE_CATCH ( std::exception& exc )
//    {
//        std::cerr << "Exception: "<< exc.what() << std::endl;
//        cpptrace::from_current_exception().print();
//    }
}
// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
