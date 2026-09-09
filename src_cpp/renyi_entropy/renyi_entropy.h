//
// Created by hynek on 5.3.25.
//
#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <format>
#include <iostream>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <numeric>
#include <sstream>
#include <vector>
#include <ranges>
#include <any>
#include <thread>
#include <random>
#include <memory>
#include <type_traits>
#include <typeindex>
#include <typeinfo>
#include <generator>
#include <latch>
#include <cmath>
#include <numbers>

#include <real2DFFT.H>
#include <complexFFT.H>

#include <boost/log/trivial.hpp>

#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/gamma.hpp>

#include <boost/serialization/map.hpp>

#include <boost/tuple/tuple.hpp>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/zstd.hpp>
#include <boost/iostreams/filter/lzma.hpp>
#include <boost/iostreams/filter/bzip2.hpp>

#include <eigen3/Eigen/Core>

#include "KDTree.cpp"

#include "msgpack.hpp"

namespace bio = boost::iostreams;

struct hash_tuple {

    template <class T1, class T2>

    size_t operator()(
        const std::tuple<T1, T2>& x)
        const
    {
        return std::get<0>(x)
               ^ std::get<1>(x);
    }
};

namespace renyi_entropy
{

const auto microseconds_in_second =
        static_cast<double> (std::chrono::duration_cast<std::chrono::microseconds> (1s).count());

struct zstd_ostream : boost::iostreams::filtering_ostream
{
    zstd_ostream(std::ostream& os)
    {
        bio::zstd_params zstd_params ( 13 );
        bio::filtering_ostream::push( bio::zstd_compressor{zstd_params} );
        bio::filtering_ostream::push( os );
    }
};

inline extern const std::string conditional_renyi_entropy_label{"CRE"};
inline extern const std::string renyi_entropy_X_present_history_label{"RE_ph_X"};
inline extern const std::string renyi_entropy_XY_history_label{"RE_h_XY"};
inline extern const std::string joint_renyi_entropy_label{"RE_ph_XY"};
inline extern const std::string renyi_entropy_X_history_label{"RE_h_X"};
inline extern const std::string average_runs{"average_runs"};
inline extern const std::string average_runs_neighbors{"average_runs_neighbors"};
inline extern const std::vector<const std::string*> RTE_types {&conditional_renyi_entropy_label, &renyi_entropy_X_present_history_label, &renyi_entropy_XY_history_label, &joint_renyi_entropy_label, &renyi_entropy_X_history_label};

template<typename Scalar, typename Matrix>
inline std::vector< std::vector<Scalar> > fromEigenMatrix ( const Matrix & M )
{
    std::vector< std::vector<Scalar> > m;
    m.resize ( M.cols(), std::vector<Scalar> ( M.rows(), 0 ) );
    for ( size_t i = 0; i < m.size(); i++ )
        for ( size_t j = 0; j < m.front().size(); j++ ) {
            m[i][j] = M ( j, i );
        }
    return m;
}

template <typename TYPE>
TYPE volume_of_hypersphere ( TYPE radius, TYPE metric, unsigned int dimension )
{
    const TYPE static_part = pow ( 2 / metric, dimension - 1 );
    const TYPE radius_part = pow ( radius, dimension ) / dimension;
    const TYPE complete_angular_part = 2 * boost::math::tgamma ( 1 / metric ) *
                                       boost::math::tgamma ( 1 / metric ) /
                                       boost::math::tgamma ( 2 / metric );

    TYPE noncomplete_angular_part =
        pow ( std::numbers::pi_v<double>, ( dimension - 2 ) / 2. );
    for ( unsigned int angular_momentum_index = 1;
            angular_momentum_index < dimension - 1; ++angular_momentum_index ) {
        auto nominator = boost::math::lgamma ( 0.5 + ( dimension - ( angular_momentum_index + 1 ) ) / metric );
        auto denominator= boost::math::lgamma ( ( dimension + metric - ( angular_momentum_index + 1 ) ) / metric );
        noncomplete_angular_part *=
            exp ( nominator - denominator );
    }

    return static_part * radius_part * complete_angular_part *
           noncomplete_angular_part;
}

template <typename TYPE>
class renyi_entropy
{
public:
    typedef TYPE renyi_key_type;
    typedef std::map<renyi_key_type, TYPE> renyi_entropy_storage;
    typedef std::vector<renyi_entropy_storage> renyi_entropy_storage_collection;
    typedef std::unordered_map<std::string, std::vector<renyi_entropy_storage>> conditional_renyi_entropy_strorage;
    typedef std::tuple<std::vector<unsigned int>, std::vector<unsigned int>, std::vector<unsigned int>> collection_conditional_information_transfer_key_type;
    typedef std::tuple<bool, bool, bool, unsigned int> conditional_information_transfer_key_type;
    typedef std::tuple<bool, bool, bool, bool, unsigned int> result_conditional_information_transfer_key_type;
    typedef std::tuple<bool, bool, bool, unsigned int, std::vector<unsigned int>, std::vector<unsigned int>, std::vector<unsigned int>, std::string, unsigned int> result_RTE_key;
    typedef std::tuple<bool, bool, bool, bool, unsigned int, std::vector<unsigned int>, std::vector<unsigned int>, std::vector<unsigned int>, std::string, std::string> collection_RTE_results_key;
    typedef std::map<conditional_information_transfer_key_type, conditional_renyi_entropy_strorage> result_conditional_information_transfer_type;
    typedef std::map<result_conditional_information_transfer_key_type, renyi_entropy_storage> processed_conditional_information_transfer_type;
    typedef std::map<collection_conditional_information_transfer_key_type, result_conditional_information_transfer_type> collection_result_conditional_information_transfer_type;
    typedef std::map<result_RTE_key, renyi_entropy_storage> result_RTE_t;
    typedef std::map<collection_conditional_information_transfer_key_type, processed_conditional_information_transfer_type> collection_processed_conditional_information_transfer_type;
    typedef std::map<std::string, collection_processed_conditional_information_transfer_type> average_result_conditional_information_transfer_type;
    typedef std::map<collection_RTE_results_key, renyi_entropy_storage> type_average_result_conditional_information_transfer_type;
    typedef std::tuple<collection_result_conditional_information_transfer_type, type_average_result_conditional_information_transfer_type> storage_RTE;

    renyi_entropy()
        : _multithreading ( false ), _number_of_threads(std::thread::hardware_concurrency())
    {}

    ~renyi_entropy()
    {}

    KDTree<TYPE> create_KDtree ( std::vector<std::vector<TYPE>> &dataset )
    {
        KDTree<TYPE> tree ( dataset );
        return tree;
    }

    KDTree<TYPE> create_KDtree ( const Eigen::MatrixXd &dataset )
    {
        KDTree<TYPE> tree ( dataset );
        return tree;
    }

    void calculation()
    {
        std::vector<double> dataset_x;
        std::vector<std::vector<double>> distances;
        for ( auto const alpha : _alphas ) {
            if ( alpha == 1.0 ) {
                auto result =
                    entropy_sum_Shannon_LeonenkoProzanto ( dataset_x, distances );
            } else {
                auto entropy_sum =
                    entropy_sum_generic_LeonenkoProzanto ( dataset_x, distances, alpha );
            }
        }
    }

    void entropy_sum_Renyi_LeonenkoProzanto (
        renyi_entropy_storage_collection &results,
        int dimension_of_data, const std::vector<std::vector<TYPE>> &distances,
        const TYPE alpha, bool log_calculation )
    {
        std::map<unsigned int, TYPE> entropy;
        TYPE one_minus_alpha = TYPE ( 1.0 ) - alpha;

        for ( const unsigned int use_index : GetIndices() ) {
            auto const number_of_data =
                std::accumulate ( distances.begin(), distances.end(), 0,
            [=] ( unsigned int count, auto const &element ) {
                if ( element[use_index] > 0.0 ) {
                    return ++count;
                } else {
                    return count;
                }
            } );

            auto maximum_distance_element =
                std::max_element ( distances.cbegin(), distances.cend(),
            [use_index] ( const auto &first, const auto &second ) {
                auto first_data = first[use_index];
                auto second_data = second[use_index];
                return first_data > second_data ? false : true;
            } );
            auto maximum_distance = ( *maximum_distance_element ) [use_index];
            TYPE exponent = dimension_of_data * one_minus_alpha;

            auto sum_of_power_of_distances = std::accumulate (
                                                 distances.cbegin(), distances.cend(), static_cast<TYPE> ( 0 ),
                                                 [use_index, maximum_distance, exponent, this] ( TYPE count,
            auto &element ) {
                auto power_base = element[use_index];
                if ( power_base > 0.0 ) {
                    return count + _power ( power_base / maximum_distance, exponent );
                } else {
                    return count;
                }
            } );

            auto index_exponent = one_minus_alpha + TYPE ( use_index );
            if ( log_calculation ) {
                auto multiplicator =
                    boost::math::lgamma ( use_index ) -
                    boost::math::lgamma ( index_exponent ) -
                    one_minus_alpha * boost::math::lgamma ( dimension_of_data / 2. + 1. ) +
                    exponent / 2. * _logarithm ( TYPE ( std::numbers::pi_v<double> ) ) +
                    exponent * _logarithm ( maximum_distance ) -
                    _logarithm ( number_of_data ) +
                    one_minus_alpha * _logarithm ( number_of_data - 1 );
                {
                    std::lock_guard<std::recursive_mutex> lock ( _result_mutex );
                    results[use_index-1][alpha] =
                        _exp ( multiplicator + _logarithm ( sum_of_power_of_distances ) );
                }
            } else {
                {
                    std::lock_guard<std::recursive_mutex> lock ( _result_mutex );
                    results[use_index-1][alpha] =
                        sum_of_power_of_distances / number_of_data *
                        pow ( number_of_data - 1, one_minus_alpha ) *
                        boost::math::tgamma ( use_index ) /
                        boost::math::tgamma ( index_exponent ) *
                        _power ( TYPE ( std::numbers::pi_v<double> ), exponent / 2 ) /
                        _power ( boost::math::tgamma ( dimension_of_data / 2. + 1. ),
                                 one_minus_alpha );
                }
            }
        }
    }

    void entropy_sum_Shannon_LeonenkoProzanto (
        renyi_entropy_storage_collection &results,
        int dimension_of_data, std::vector<std::vector<TYPE>> &distances,
        bool log_calculation )
    {
        for ( const unsigned int use_index : GetIndices() ) {
            auto const number_of_data =
                std::accumulate ( distances.cbegin(), distances.cend(), 0,
            [=] ( unsigned int count, auto const &element ) {
                if ( element[use_index] > 0.0 ) {
                    return ++count;
                } else {
                    return count;
                }
            } );

            if ( log_calculation ) {
                const auto V_m =
                    _power ( std::numbers::pi_v<double>, dimension_of_data / 2.0 ) /
                    boost::math::tgamma ( dimension_of_data / 2.0 + 1 );
                const auto digamma = boost::math::digamma ( use_index );
                const auto multiplicator = ( number_of_data - 1 ) * _exp ( -digamma );
                auto const sum_of_weighted_distances = std::accumulate (
                        distances.begin(), distances.end(), static_cast<TYPE> ( 0 ),
                [&] ( TYPE count, auto const &element ) {
                    const auto base = element[use_index];
                    if ( base > 0.0 ) {
                        return count + _logarithm ( multiplicator * V_m *
                                                    +_power ( base, dimension_of_data ) );
                    } else {
                        return count;
                    }
                } );
                {
                    std::lock_guard<std::recursive_mutex> lock ( _result_mutex );
                    results[use_index-1][static_cast<TYPE> ( 1.0 )] =
                        sum_of_weighted_distances / number_of_data;
                }
            } else {
                auto const sum_of_log_distances = std::accumulate (
                                                      distances.begin(), distances.end(), static_cast<TYPE> ( 0 ),
                [use_index, this] ( TYPE count, auto const &element ) {
                    const auto log_base = element[use_index];
                    if ( log_base > 0.0 ) {
                        return count + _logarithm ( log_base );
                    } else {
                        return count;
                    }
                } );

                auto const addition_to_entropy =
                    dimension_of_data * sum_of_log_distances / number_of_data;
                auto digamma = boost::math::digamma ( use_index );
                auto V_m = _power ( std::numbers::pi_v<double>, dimension_of_data / 2.0 ) /
                           boost::math::tgamma ( dimension_of_data / 2.0 + 1 );
                auto argument_log = V_m * ( number_of_data - 1 );
                auto log_volume = _logarithm ( argument_log );
                {
                    std::lock_guard<std::recursive_mutex> lock ( _result_mutex );
                    results[use_index-1][static_cast<TYPE> ( 1.0 )] =
                        addition_to_entropy + log_volume - digamma;
                }
            }
        }
    }

    renyi_entropy_storage_collection
    renyi_entropy_LeonenkoProzanto ( std::vector<std::vector<TYPE>> &dataset,
                                     double metric = 2 )
    {
        auto nearest_iterator =
            std::max_element ( GetIndices().begin(), GetIndices().end() );
        // calculation how large of index of distance needs to be calculated
        // + 1 for skipping 0-neighbor
        auto nearest = ( *nearest_iterator ) + 1;
        renyi_entropy_storage_collection results(GetIndices().back());
        int dimension_of_data = dataset[0].size();
        auto start_calculation = std::chrono::high_resolution_clock::now();
        auto kdtree = create_KDtree ( dataset );
        auto end_tree_construction_calculation =
            std::chrono::high_resolution_clock::now();

        // calculate distances
        std::vector<std::vector<TYPE>> distances;
        calculate_distances_from_each_point ( kdtree, dataset, nearest, distances,
                                              metric );

        // remove kdtree to same memory
        auto end_distance_calculation = std::chrono::high_resolution_clock::now();

        // calculate renyi entropy
        auto log_calculation = true;
        auto one = static_cast<TYPE> ( 1 );
        if ( ! GetMultithreading () ) {
            for ( const auto alpha : GetAlphas() ) {
                if ( alpha == one ) {
                    entropy_sum_Shannon_LeonenkoProzanto ( results, dimension_of_data,
                                                        distances, log_calculation );
                } else {
                    entropy_sum_Renyi_LeonenkoProzanto ( results, dimension_of_data,
                                                        distances, alpha, log_calculation );
                }
            }
        } else {
            const auto processor_count = GetNumberOfThreads ();
            const auto maximal_number_of_jobs = GetAlphas().size();
            const auto real_thread_count = ( processor_count >= maximal_number_of_jobs ? maximal_number_of_jobs : processor_count );
            const auto jobs_to_process { maximal_number_of_jobs / real_thread_count };
            const auto excess_jobs_to_process { maximal_number_of_jobs % real_thread_count };
            std::vector<std::thread> threads ( real_thread_count );
            std::latch work_done{static_cast<long int>(real_thread_count)};
            for ( int thread_count : std::ranges::iota_view{0U, real_thread_count} ) {
                std::thread thread_worker (
                [&] ( int thread_count ) {
                    const int start_job = jobs_to_process * thread_count + ( ( 0 <= thread_count ) && ( thread_count < excess_jobs_to_process ) ? thread_count : excess_jobs_to_process );
                    const int end_job = ( jobs_to_process * ( thread_count + 1 ) ) + ( ( 0 <= thread_count ) && ( thread_count < excess_jobs_to_process ) ? thread_count + 1 : excess_jobs_to_process );
                    //std::cout << thread_count << " " << start_job << " " << end_job << std::endl;
                    for ( int alpha_index : std::ranges::iota_view{start_job, end_job} ) {
                        const auto alpha = GetAlphas() [alpha_index];
                        if ( alpha == one ) {
                            entropy_sum_Shannon_LeonenkoProzanto ( results, dimension_of_data, distances,
                                                         log_calculation );
                        } else {
                            entropy_sum_Renyi_LeonenkoProzanto ( results, dimension_of_data, distances, alpha,
                                                       log_calculation );
                        }
                    }
                    work_done.count_down();
                }, thread_count
                );
                threads[thread_count] = std::move ( thread_worker );
            }

            work_done.wait();

            for ( int thread_count : std::ranges::iota_view{0U, real_thread_count} ) {
                threads[thread_count].join();
            }
            threads.clear();
        }

        if ( ! GetMultithreading () ) {
            // calculation of entropy for Renyi case
            // Shannon case already holds entropy
            for ( auto &item_neighbor : results ) {
                for (auto & item: item_neighbor)
                {
                    auto alpha = item.first;
                    if ( alpha != one ) {
                        auto entropy = _logarithm ( item.second ) / ( 1 - alpha );
                        item.second = entropy;
                    }
                }
            }
        } else {
            auto ks = std::views::keys(results.back());
            std::vector<renyi_key_type> keys{ ks.begin(), ks.end() };

            const auto processor_count = GetNumberOfThreads ();
            const auto maximal_number_of_jobs = GetAlphas().size() * GetIndices().size();
            const auto real_thread_count = ( processor_count >= maximal_number_of_jobs ? maximal_number_of_jobs : processor_count );
            const auto jobs_to_process { maximal_number_of_jobs / real_thread_count };
            const auto excess_jobs_to_process { maximal_number_of_jobs % real_thread_count };
            std::vector<std::thread> threads ( processor_count );
            std::latch work_done{processor_count};
            for ( int thread_count : std::ranges::iota_view{0U, processor_count} ) {
                std::thread thread_worker (
                [&] ( int thread_count ) {
                    const int start_job = jobs_to_process * thread_count + ( ( 0 <= thread_count ) && ( thread_count < excess_jobs_to_process ) ? thread_count : excess_jobs_to_process );
                    const int end_job = ( jobs_to_process * ( thread_count + 1 ) ) + ( ( 0 <= thread_count ) && ( thread_count < excess_jobs_to_process ) ? thread_count + 1 : excess_jobs_to_process );
                    //std::cout << thread_count << " " << start_job << " " << end_job << std::endl;
                    for ( int key_index : std::ranges::iota_view{start_job, end_job} ) {
                        const auto index_neighbor = key_index % GetAlphas().size();
                        const auto index_alpha = key_index / GetAlphas().size();

                        const auto alpha = GetAlphas()[index_alpha];
                        const auto neighbor = GetIndices()[index_neighbor] - 1;
                        if ( alpha != one ) {
                            auto &RTE = results[neighbor][alpha];
                            if (RTE>0)
                            {
                                auto entropy = _logarithm ( results[neighbor][alpha] ) / ( 1 - alpha );
                                RTE = entropy;
                            }
                        }
                    }
                    work_done.count_down();
                }, thread_count
                );
                threads[thread_count] = std::move ( thread_worker );
            }

            work_done.wait();

            for ( int thread_count : std::ranges::iota_view{0U, jobs_to_process} ) {
                threads[thread_count].join();
            }
            threads.clear();
        }

        auto end_entropy_calculation = std::chrono::high_resolution_clock::now();

        auto elapsed_tree_construction =
            end_tree_construction_calculation - start_calculation;
        auto elapsed_distance_calculation =
            end_distance_calculation - end_tree_construction_calculation;
        auto elapsed_entropy_calculation =
            end_entropy_calculation - end_distance_calculation;

        auto microseconds_tree_construction =
            std::chrono::duration_cast<std::chrono::microseconds> (
                elapsed_tree_construction )
            .count();
        auto microseconds_distance_calculation =
            std::chrono::duration_cast<std::chrono::microseconds> (
                elapsed_distance_calculation )
            .count();
        auto microseconds_entropy_calculation =
            std::chrono::duration_cast<std::chrono::microseconds> (
                elapsed_entropy_calculation )
            .count();

        BOOST_LOG_TRIVIAL ( trace )
                << "Tree constructed in " << microseconds_tree_construction / microseconds_in_second << " seconds";
        BOOST_LOG_TRIVIAL ( trace )
                << "Distances calculated in " << microseconds_distance_calculation / microseconds_in_second << " seconds";
        BOOST_LOG_TRIVIAL ( trace )
                << "Entropy calculated in " << microseconds_entropy_calculation / microseconds_in_second << " seconds";
        return results;
    }

    void entropy_sum_Renyi_metric (
        renyi_entropy_storage_collection &results,
        int dimension_of_data, const std::vector<std::vector<TYPE>> &distances,
        const TYPE alpha, bool log_calculation, TYPE metric )
    {
        std::map<unsigned int, TYPE> entropy;
        TYPE one_minus_alpha = TYPE ( 1.0 ) - alpha;

        for ( const auto [index_used_index, use_index] : std::views::enumerate(GetIndices()) ) {
            auto const number_of_data =
                std::accumulate ( distances.begin(), distances.end(), 0,
            [=] ( unsigned int count, auto const &element ) {
                if ( element[use_index] > 0.0 ) {
                    return ++count;
                } else {
                    return count;
                }
            } );

            const auto maximum_distance_element =
                std::max_element ( distances.cbegin(), distances.cend(),
            [use_index] ( const auto &first, const auto &second ) {
                auto first_data = first[use_index];
                auto second_data = second[use_index];
                return first_data > second_data ? false : true;
            } );
            const auto maximum_distance = ( *maximum_distance_element ) [use_index];
            const TYPE exponent = dimension_of_data * one_minus_alpha;

            const auto sum_of_power_of_distances = std::accumulate (
                    distances.cbegin(), distances.cend(), static_cast<TYPE> ( 0 ),
                    [use_index, maximum_distance, exponent, this] ( TYPE count,
            auto &element ) {
                auto power_base = element[use_index];
                if ( power_base > 0.0 ) {
                    return count + _power ( power_base / maximum_distance, exponent );
                } else {
                    return count;
                }
            } );

            const auto dimensional_multiplicator =
                volume_of_hypersphere ( TYPE ( 1 ), metric, dimension_of_data );
            const auto index_exponent = one_minus_alpha + TYPE ( use_index );
            auto floor_index_exponent = floor(index_exponent);
            if ( log_calculation ) {
                if ( (index_exponent <= 0) && (index_exponent == floor_index_exponent) )
                {
                    results[index_used_index][alpha] = NAN;
                }
                else
                {
                    const auto multiplicator = one_minus_alpha * _logarithm ( dimensional_multiplicator ) +
                                            boost::math::lgamma ( use_index ) -
                                            boost::math::lgamma ( index_exponent ) +
                                            exponent * _logarithm ( maximum_distance ) -
                                            _logarithm ( number_of_data ) +
                                            one_minus_alpha * _logarithm ( number_of_data - 1 );
                    std::lock_guard<std::recursive_mutex> lock ( _result_mutex );
                    results[index_used_index][alpha] = _exp ( multiplicator + _logarithm ( sum_of_power_of_distances ) );
                }
            } else {
                std::lock_guard<std::recursive_mutex> lock ( _result_mutex );
                results[index_used_index][alpha] =
                    sum_of_power_of_distances / number_of_data *
                    pow ( number_of_data - 1, one_minus_alpha ) *
                    boost::math::tgamma ( use_index ) /
                    boost::math::tgamma ( index_exponent ) * _power ( dimensional_multiplicator, one_minus_alpha );
            }
        }
    }

    void entropy_sum_Shannon_metric (
        renyi_entropy_storage_collection &results,
        int dimension_of_data, std::vector<std::vector<TYPE>> &distances,
        bool log_calculation, TYPE metric )
    {
        for ( const auto [index_used_index, use_index] : std::views::enumerate(GetIndices()) ) {
            auto const number_of_data =
                std::accumulate ( distances.cbegin(), distances.cend(), 0,
            [=] ( unsigned int count, auto const &element ) {
                if ( element[use_index] > 0.0 ) {
                    return ++count;
                } else {
                    return count;
                }
            } );

            const auto V_m =
                volume_of_hypersphere ( TYPE ( 1 ), metric, dimension_of_data );
            if ( log_calculation ) {
                const auto digamma = boost::math::digamma ( use_index );
                const auto multiplicator = ( number_of_data - 1 ) * _exp ( -digamma );
                auto const sum_of_weighted_distances = std::accumulate (
                        distances.begin(), distances.end(), static_cast<TYPE> ( 0 ),
                [&] ( TYPE count, auto const &element ) {
                    const auto base = element[use_index];
                    if ( base > 0.0 ) {
                        return count + _logarithm ( multiplicator * V_m *
                                                    +_power ( base, dimension_of_data ) );
                    } else {
                        return count;
                    }
                } );
                {
                    std::lock_guard<std::recursive_mutex> lock ( _result_mutex );
                    results[index_used_index] [static_cast<TYPE> ( 1.0 ) ] = sum_of_weighted_distances / number_of_data;
                }
            } else {
                auto const sum_of_log_distances = std::accumulate (
                                                      distances.begin(), distances.end(), static_cast<TYPE> ( 0 ),
                [use_index, this] ( TYPE count, auto const &element ) {
                    const auto log_base = element[use_index];
                    if ( log_base > 0.0 ) {
                        return count + _logarithm ( log_base );
                    } else {
                        return count;
                    }
                } );

                const auto addition_to_entropy =
                    dimension_of_data * sum_of_log_distances / number_of_data;
                const auto digamma = boost::math::digamma ( use_index );
                auto argument_log = V_m * ( number_of_data - 1 );
                auto log_volume = _logarithm ( argument_log );
                {
                    std::lock_guard<std::recursive_mutex> lock ( _result_mutex );
                    results[index_used_index][ static_cast<TYPE> ( 1.0 )] = addition_to_entropy + log_volume - digamma;
                }
            }
        }
    }


    renyi_entropy_storage_collection
    renyi_entropy_metric ( const Eigen::MatrixXd &dataset, double metric = 2 )
    {
        auto nearest_iterator =
            std::max_element ( GetIndices().begin(), GetIndices().end() );
        // calculation how large of index of distance needs to be calculated
        // + 1 for skipping 0-neighbor
        auto nearest = ( *nearest_iterator ) + 1;
        renyi_entropy_storage_collection results(GetIndices().size());
        const int dimension_of_data = dataset.rows();
        auto start_calculation = std::chrono::high_resolution_clock::now();
        auto kdtree = create_KDtree ( dataset );
        auto end_tree_construction_calculation =
            std::chrono::high_resolution_clock::now();

        // calculate distances
        std::vector<std::vector<TYPE>> distances;
        calculate_distances_from_each_point ( kdtree, dataset, nearest, distances,
                                              metric );

        // remove kdtree to same memory
        auto end_distance_calculation = std::chrono::high_resolution_clock::now();

        // calculate renyi entropy
        auto log_calculation = true;
        auto one = static_cast<TYPE> ( 1 );
        if ( ! GetMultithreading () ) {
            for ( const auto alpha : GetAlphas() ) {
                if ( alpha == one ) {
                    entropy_sum_Shannon_metric ( results, dimension_of_data, distances,
                                                log_calculation, metric );
                } else {
                    entropy_sum_Renyi_metric ( results, dimension_of_data, distances, alpha,
                                            log_calculation, metric );
                }
            }
        } else {
            const auto processor_count = GetNumberOfThreads ();
            const auto maximal_number_of_jobs = GetAlphas().size();
            const auto real_thread_count = ( processor_count >= maximal_number_of_jobs ? maximal_number_of_jobs : processor_count );
            const auto jobs_to_process { maximal_number_of_jobs / real_thread_count };
            const auto excess_jobs_to_process { maximal_number_of_jobs % real_thread_count };
            std::vector<std::thread> threads ( real_thread_count );
            std::latch work_done{static_cast<long>(real_thread_count)};
            for ( int thread_count : std::ranges::iota_view{0U, real_thread_count} ) {
                std::thread thread_worker (
                [&] ( int thread_count ) {
                    const int start_job = jobs_to_process * thread_count + ( ( 0 <= thread_count ) && ( thread_count < excess_jobs_to_process ) ? thread_count : excess_jobs_to_process );
                    const int end_job = ( jobs_to_process * ( thread_count + 1 ) ) + ( ( 0 <= thread_count ) && ( thread_count < excess_jobs_to_process ) ? thread_count + 1 : excess_jobs_to_process );
                    //std::cout << thread_count << " " << start_job << " " << end_job << std::endl;
                    for ( int alpha_index : std::ranges::iota_view{start_job, end_job} ) {
                        const auto alpha = GetAlphas() [alpha_index];
                        if ( alpha == one ) {
                            entropy_sum_Shannon_metric ( results, dimension_of_data, distances,
                                                         log_calculation, metric );
                        } else {
                            entropy_sum_Renyi_metric ( results, dimension_of_data, distances, alpha,
                                                       log_calculation, metric );
                        }
                    }
                    work_done.count_down();
                }, thread_count
                );
                threads[thread_count] = std::move ( thread_worker );
            }

            work_done.wait();

            for ( int thread_count : std::ranges::iota_view{0U, real_thread_count} ) {
                threads[thread_count].join();
            }
            threads.clear();
        }

        // fitted compensation to large and small alphas of Renyi entropy calculation
        renyi_entropy_LeonenkoProzanto_compensation ( results, dimension_of_data );
        auto end_entropy_calculation = std::chrono::high_resolution_clock::now();

        auto elapsed_tree_construction =
            end_tree_construction_calculation - start_calculation;
        auto elapsed_distance_calculation =
            end_distance_calculation - end_tree_construction_calculation;
        auto elapsed_entropy_calculation =
            end_entropy_calculation - end_distance_calculation;

        auto microseconds_tree_construction =
            std::chrono::duration_cast<std::chrono::microseconds> (
                elapsed_tree_construction )
            .count();
        auto microseconds_distance_calculation =
            std::chrono::duration_cast<std::chrono::microseconds> (
                elapsed_distance_calculation )
            .count();
        auto microseconds_entropy_calculation =
            std::chrono::duration_cast<std::chrono::microseconds> (
                elapsed_entropy_calculation )
            .count();

        BOOST_LOG_TRIVIAL ( trace )
                << "Tree constructed in " << microseconds_tree_construction / microseconds_in_second << " seconds";
        BOOST_LOG_TRIVIAL ( trace )
                << "Distances calculated in " << microseconds_distance_calculation / microseconds_in_second << " seconds";
        BOOST_LOG_TRIVIAL ( trace )
                << "Entropy calculated in " << microseconds_entropy_calculation / microseconds_in_second << " seconds";
        return results;
    }

    renyi_entropy_storage_collection renyi_entropy_metric ( std::vector<std::vector<TYPE>> &dataset, double metric = 2 )
    {
        auto nearest_iterator =
            std::max_element ( GetIndices().begin(), GetIndices().end() );
        // calculation how large of index of distance needs to be calculated
        // + 1 for skipping 0-neighbor
        auto nearest = ( *nearest_iterator ) + 1;
        renyi_entropy_storage_collection results;
        const int dimension_of_data = dataset[0].size();
        auto start_calculation = std::chrono::high_resolution_clock::now();
        auto kdtree = create_KDtree ( dataset );
        auto end_tree_construction_calculation =
            std::chrono::high_resolution_clock::now();

        // calculate distances
        std::vector<std::vector<TYPE>> distances;
        calculate_distances_from_each_point ( kdtree, dataset, nearest, distances,
                                              metric );

        // remove kdtree to same memory
        auto end_distance_calculation = std::chrono::high_resolution_clock::now();

        // calculate renyi entropy
        auto log_calculation = true;
        auto one = static_cast<TYPE> ( 1 );
        if ( ! GetMultithreading () ) {
            for ( const auto alpha : GetAlphas() ) {
                if ( alpha == one ) {
                    entropy_sum_Shannon_metric ( results, dimension_of_data, distances,
                                                 log_calculation, metric );
                } else {
                    entropy_sum_Renyi_metric ( results, dimension_of_data, distances, alpha,
                                               log_calculation, metric );
                }
            }
        } else {
            const auto maximal_number_of_jobs = GetAlphas().size();
            const auto processor_count = GetNumberOfThreads ();
            const auto real_thread_count = ( processor_count >= maximal_number_of_jobs ? maximal_number_of_jobs : processor_count );
            const auto jobs_to_process { maximal_number_of_jobs / real_thread_count };
            const auto excess_jobs_to_process { maximal_number_of_jobs % real_thread_count };
            std::vector<std::thread> threads ( real_thread_count );
            std::latch work_done{static_cast<long>(real_thread_count)};
            for ( int thread_count : std::ranges::iota_view{0U, real_thread_count} ) {
                std::thread thread_worker (
                [&] ( int thread_count ) {
                    const int start_job = jobs_to_process * thread_count + ( ( 0 <= thread_count ) && ( thread_count < excess_jobs_to_process ) ? thread_count : excess_jobs_to_process );
                    const int end_job = ( jobs_to_process * ( thread_count + 1 ) ) + ( ( 0 <= thread_count ) && ( thread_count < excess_jobs_to_process ) ? thread_count + 1 : excess_jobs_to_process );
                    //std::cout << thread_count << " " << start_job << " " << end_job << std::endl;
                    for ( int alpha_index : std::ranges::iota_view{start_job, end_job} ) {
                        const auto alpha = GetAlphas() [alpha_index];
                        if ( alpha == one ) {
                            entropy_sum_Shannon_metric ( results, dimension_of_data, distances,
                                                         log_calculation, metric );
                        } else {
                            entropy_sum_Renyi_metric ( results, dimension_of_data, distances, alpha,
                                                       log_calculation, metric );
                        }
                    }
                    work_done.count_down();
                }, thread_count
                );
                threads[thread_count] = std::move ( thread_worker );
            }

            work_done.wait();

            for ( int thread_count : std::ranges::iota_view{0U, real_thread_count} ) {
                threads[thread_count].join();
            }
            threads.clear();
        }

        renyi_entropy_LeonenkoProzanto_compensation ( results, dimension_of_data );
        auto end_entropy_calculation = std::chrono::high_resolution_clock::now();

        auto elapsed_tree_construction =
            end_tree_construction_calculation - start_calculation;
        auto elapsed_distance_calculation =
            end_distance_calculation - end_tree_construction_calculation;
        auto elapsed_entropy_calculation =
            end_entropy_calculation - end_distance_calculation;

        auto microseconds_tree_construction =
            std::chrono::duration_cast<std::chrono::microseconds> (
                elapsed_tree_construction )
            .count();
        auto microseconds_distance_calculation =
            std::chrono::duration_cast<std::chrono::microseconds> (
                elapsed_distance_calculation )
            .count();
        auto microseconds_entropy_calculation =
            std::chrono::duration_cast<std::chrono::microseconds> (
                elapsed_entropy_calculation )
            .count();

        BOOST_LOG_TRIVIAL ( trace )
                << "Tree constructed in " << microseconds_tree_construction / microseconds_in_second << " seconds";
        BOOST_LOG_TRIVIAL ( trace )
                << "Distances calculated in " << microseconds_distance_calculation / microseconds_in_second << " seconds";
        BOOST_LOG_TRIVIAL ( trace )
                << "Entropy calculated in " << microseconds_entropy_calculation / microseconds_in_second << " seconds";

        return results;
    }

    inline void renyi_entropy_LeonenkoProzanto_compensation ( renyi_entropy_storage_collection &results, const int dimension_of_data )
    {
        auto one = static_cast<TYPE> ( 1 );
        // calculation of entropy for Renyi case
        // Shannon case already holds entropy
        const auto fixed_compesation = -0.001921 * dimension_of_data + 1.02726;

        const auto exponent1 = ( 0.167 / dimension_of_data - 0.181 );
        const auto multiplication_factor1 = 0.95;
        const auto threshold1 = exp ( log ( fixed_compesation ) - log ( multiplication_factor1 ) / exponent1 );

        const auto exponent2 = ( -0.0015 * dimension_of_data + 0.00067 );
        const auto multiplication_factor2 = 0.99;
        const auto threshold2 = exp ( log ( fixed_compesation ) - log ( multiplication_factor2 ) / exponent2 );

        for ( auto &item_neighbor : results ) {
            for( auto &item : item_neighbor )
            {
                auto alpha = item.first;
                const auto compensation1 = ( alpha < threshold1 ) ? ( multiplication_factor1 * pow ( alpha, exponent1 ) ) : 1;
                const auto compensation2 = ( alpha > threshold2 ) ? ( multiplication_factor2 * pow ( alpha, exponent2 ) ) : 1;

                if ( alpha != one ) {
                    auto entropy = _logarithm ( item.second ) / ( 1 - alpha ) * compensation1 * compensation2 * fixed_compesation;
                    item.second = entropy;
                } else {
                    auto entropy = item.second * compensation1 * compensation2 * fixed_compesation;
                    item.second = entropy;
                }
            }
        }
    }

    inline void calculate_distances_from_each_point (
        KDTree<TYPE> &kdtree, std::vector<std::vector<TYPE>> &dataset,
        unsigned int nearest, std::vector<std::vector<TYPE>> &distances,
        const TYPE metric )
    {
        // calculate distances
        distances.resize ( dataset.size() );
        if ( ! GetMultithreading () ) {
            for ( const auto &[index, point] : std::views::enumerate ( dataset ) ) {
                distances[index].resize ( nearest );
                auto points = kdtree.nearest_points ( point, nearest );
                for ( int nearest_index = 0; nearest_index < nearest; ++nearest_index ) {
                    distances[index][nearest_index] = ( distance ( point, points[nearest_index], metric ) );
                }
            }
        } else {
            const auto maximal_number_of_jobs = dataset.size();
            const auto processor_count = GetNumberOfThreads ();
            const auto real_thread_count = ( processor_count >= maximal_number_of_jobs ? maximal_number_of_jobs : processor_count );
            const auto jobs_to_process { maximal_number_of_jobs / real_thread_count };
            const auto excess_jobs_to_process { maximal_number_of_jobs % real_thread_count };
            std::vector<std::thread> threads ( real_thread_count );
            std::latch work_done{static_cast<long>(real_thread_count)};
            for ( int thread_count : std::ranges::iota_view{0U, real_thread_count} ) {
                std::thread thread_worker (
                [&real_thread_count, &jobs_to_process, &dataset, &kdtree, &nearest, &distances, &metric, &work_done, &maximal_number_of_jobs, &excess_jobs_to_process] ( int thread_count ) {
                    const int start_job = jobs_to_process * thread_count + ( ( 0 <= thread_count ) && ( thread_count < excess_jobs_to_process ) ? thread_count : excess_jobs_to_process );
                    const int end_job = ( jobs_to_process * ( thread_count + 1 ) ) + ( ( 0 <= thread_count ) && ( thread_count < excess_jobs_to_process ) ? thread_count + 1 : excess_jobs_to_process );
                    //std::cout << thread_count << " " << start_job << " " << end_job << std::endl;
                    for ( int point_index : std::ranges::iota_view{start_job, end_job} ) {
                        distances[point_index].resize ( nearest );
                        const auto point = dataset[point_index];
                        auto points = kdtree.nearest_points ( point, nearest );
                        for ( int nearest_index = 0; nearest_index < nearest; ++nearest_index ) {
                            distances[point_index][nearest_index] = ( distance ( point, points[nearest_index], metric ) );
                        }
                    }
                    work_done.count_down();
                }, thread_count
                );
                threads[thread_count] = std::move ( thread_worker );
            }

            work_done.wait();

            for ( int thread_count : std::ranges::iota_view{0U, real_thread_count} ) {
                threads[thread_count].join();
            }
            threads.clear();
        }
    }

    inline void calculate_distances_from_each_point (
        KDTree<TYPE> &kdtree, const Eigen::MatrixXd &dataset,
        unsigned int nearest, std::vector<std::vector<TYPE>> &distances,
        const TYPE metric )
    {
        // calculate distances
        const auto number_data = dataset.cols();
        distances.resize ( number_data );
        if ( ! GetMultithreading () ) {
            for ( int index = 0; index < number_data; ++index ) {
                auto item = dataset.col ( index );
                const auto start = item.data();
                auto end =  item.data() + item.rows();
                std::vector<double> point ( start, end );
                auto points = kdtree.nearest_points ( point, nearest );
                distances[index].resize(nearest);
                for ( int j = 0; j < nearest; ++j ) {
                    distances[index][j] = ( distance ( point, points[j], metric ) );
                }
            }
        }
        else
        {
            const auto maximal_number_of_jobs = dataset.cols();
            const auto processor_count = GetNumberOfThreads ();
            const auto real_thread_count = static_cast<unsigned int>( processor_count >= maximal_number_of_jobs ? maximal_number_of_jobs : processor_count );
            const auto jobs_to_process { maximal_number_of_jobs / real_thread_count };
            const auto excess_jobs_to_process { maximal_number_of_jobs - jobs_to_process * real_thread_count };
            std::vector<std::thread> threads ( real_thread_count );
            std::latch work_done{static_cast<long>(real_thread_count)};
            for ( int thread_count : std::ranges::iota_view{0U, real_thread_count} ) {
                std::thread thread_worker (
                [&real_thread_count, &jobs_to_process, &dataset, &kdtree, &nearest, &distances, &metric, &work_done, &maximal_number_of_jobs, &excess_jobs_to_process, &processor_count] ( int thread_count ) {
                    const int start_job = jobs_to_process * thread_count + ( ( 0 <= thread_count ) && ( thread_count < excess_jobs_to_process ) ? thread_count : excess_jobs_to_process );
                    const int end_job = ( jobs_to_process * ( thread_count + 1 ) ) + ( ( 0 <= thread_count ) && ( thread_count < excess_jobs_to_process ) ? thread_count + 1 : excess_jobs_to_process );
                    //std::cout << thread_count << " " << start_job << " " << end_job << std::endl;
                    for ( int point_index : std::ranges::iota_view{start_job, end_job}  ) {
                        distances[point_index].resize ( nearest );
                        auto item = dataset.col ( point_index );
                        const auto start = item.data();
                        auto end =  item.data() + item.rows();
                        std::vector<double> point ( start, end );
                        auto points = kdtree.nearest_points ( point, nearest );
                        for ( int nearest_index = 0; nearest_index < nearest; ++nearest_index ) {
                            distances[point_index][nearest_index] = ( distance ( point, points[nearest_index], metric ) );
                        }
                    }
                    work_done.count_down();
                }, thread_count
                );
                threads[thread_count] = std::move ( thread_worker );
            }

            work_done.wait();

            for ( int thread_count : std::ranges::iota_view{0U, real_thread_count} ) {
                threads[thread_count].join();
            }
            threads.clear();
        }
    }

    std::vector<unsigned int> &GetIndices()
    {
        return _indices;
    }

    const std::vector<TYPE> &GetAlphas()
    {
        return _alphas;
    }

    static TYPE distance ( const std::vector<TYPE> &a, const std::vector<TYPE> &b,
                           TYPE power = 2 )
    {
        if ( a.size() == b.size() ) {
            TYPE distance = 0;
            for ( int i = 0; i < a.size(); ++i ) {
                distance += pow ( abs ( a[i] - b[i] ), power );
            }
            return pow ( distance, 1 / power );
        }
        return NAN;
    }

    void SetAlpha ( const std::vector<TYPE> &alpha )
    {
        _alphas = alpha;
    }

    void SetIndices ( const std::vector<unsigned int> &indices )
    {
        _indices = indices;
    }

    void SetPower ( std::function<TYPE ( TYPE, TYPE ) > power )
    {
        _power = power;
    }

    void SetExp ( std::function<TYPE ( TYPE ) > exp )
    {
        _exp = exp;
    }

    void SetLog ( std::function<TYPE ( TYPE ) > log )
    {
        _logarithm = log;
    }

    static void SaveRenyiEntropy ( renyi_entropy_storage_collection &result, const std::string &directory, const std::string &file )
    {
        std::stringstream ss;
        msgpack::pack(ss, result); //renyi_entropy::renyi_entropy<calculation_type>::storage_RTE(collection_result_RTE, processing_RTE)
        std::cout << ss.str().size() << std::endl;
        boost::filesystem::path output_file =
        boost::filesystem::path(directory) / boost::filesystem::path(file);
        boost::filesystem::ofstream output_file_handler ( output_file );
        zstd_ostream zstd_compression_stream{output_file_handler};

        zstd_compression_stream << ss.str();
    }

    static void LoadRenyiEntropy ( renyi_entropy_storage_collection &result, const std::string &file )
    {
        boost::filesystem::path myFile =
            boost::filesystem::current_path() / file;
        boost::filesystem::ifstream ifs ( myFile );
        boost::archive::binary_iarchive iarch ( ifs );
        iarch >> result;
    }

    static void renyi_conditional_information_transfer ( result_RTE_t &RTE_result, const Eigen::MatrixXd &y_future, const Eigen::MatrixXd &y_history, const Eigen::MatrixXd &z_history, std::map<std::string, std::any> parameters, std::function< result_RTE_key ( std::string key, unsigned int neighbor) > key_generator )
    {
        bool enhanced_calculation {true};
        if ( parameters.contains ( "enhanced_calculation" ) ) {
            enhanced_calculation = std::any_cast<bool> ( parameters["enhanced_calculation"] );
        }

        renyi_entropy<double> calculator;
        if ( parameters.contains ( "indices" ) ) {
            auto indices = std::any_cast<std::vector<unsigned int>> ( parameters["indices"] );
            calculator.SetIndices ( indices );
        }

        if ( parameters.contains ( "alphas" ) ) {
            auto alphas = std::any_cast<std::vector<double>> ( parameters["alphas"] );
            calculator.SetAlpha ( alphas );
        }

        if ( parameters.contains ( "multithreading" ) ) {
            auto multithreading = std::any_cast<bool> ( parameters["multithreading"] );
            calculator.SetMultithreading ( multithreading );
        }
        calculator.SetExp ( [&] (double x) { return exp(x);} );
        calculator.SetLog ( [&] (double x) { return log(x);} );
        calculator.SetPower ( [&] (double x, double y) { return pow(x, y);} );

        if ( enhanced_calculation ) {
            Eigen::MatrixXd joint_dataset ( y_future.rows() + y_history.rows(), y_future.cols() );
            joint_dataset << y_future, y_history;
            auto entropy_present_X_history_X = calculator.renyi_entropy_metric ( joint_dataset, 2 );

            joint_dataset = Eigen::MatrixXd ( y_history.rows() + z_history.rows(), y_history.cols() );
            joint_dataset << y_history, z_history;
            auto entropy_history_X_history_Y = calculator.renyi_entropy_metric ( joint_dataset, 2 );

            joint_dataset = Eigen::MatrixXd ( y_future.rows() + y_history.rows() + z_history.rows(), y_future.cols() );
            joint_dataset << y_future, y_history, z_history;
            auto entropy_joint = calculator.renyi_entropy_metric ( joint_dataset, 2 );

            auto entropy_history_X = calculator.renyi_entropy_metric ( y_history, 2 );

            decltype ( entropy_present_X_history_X ) conditional_information_transfer;

            auto conditional_information_transfer_calculator = [] ( auto a, auto b, auto c, auto d ) {
                const auto alpha_a = a.first;
                const auto alpha_b = b.first;
                const auto alpha_c = c.first;
                const auto alpha_d = d.first;
                assert(alpha_a == alpha_b);
                assert(alpha_b == alpha_c);
                assert(alpha_c == alpha_d);
                return std::tuple<TYPE, TYPE>(alpha_a, a.second + b.second - c.second - d.second );
            };

#ifndef NDEBUG
            auto keys_view = std::views::keys(entropy_present_X_history_X.back());
            std::set<renyi_key_type> keys{ keys_view.begin(), keys_view.end() };
            keys_view = std::views::keys(entropy_history_X_history_Y.back());
            keys.insert( keys_view.begin(), keys_view.end() );
            keys_view = std::views::keys(entropy_joint.back());
            keys.insert( keys_view.begin(), keys_view.end() );
            keys_view = std::views::keys(entropy_history_X.back());
            keys.insert( keys_view.begin(), keys_view.end() );

            if (keys.size() != std::views::keys(entropy_present_X_history_X.back()).size())
            {
                BOOST_LOG_TRIVIAL ( debug ) << std::format("entropy_present_X_history_X {} {}", keys.size(), std::views::keys(entropy_present_X_history_X.back()).size() );
            }
            if (keys.size() != std::views::keys(entropy_history_X_history_Y.back()).size())
            {
                BOOST_LOG_TRIVIAL ( debug ) << std::format("entropy_history_X_history_Y {} {}", keys.size(), std::views::keys(entropy_history_X_history_Y.back()).size() );
            }
            if (keys.size() != std::views::keys(entropy_joint.back()).size())
            {
                BOOST_LOG_TRIVIAL ( debug ) << std::format("entropy_joint {} {}", keys.size(), std::views::keys(entropy_joint.back()).size() );
            }
            if (keys.size() != std::views::keys(entropy_history_X.back()).size())
            {
                BOOST_LOG_TRIVIAL ( debug ) << std::format("entropy_history_X {} {}", keys.size(), std::views::keys(entropy_history_X.back()).size() );
            }

            auto keys2(keys);
            for ( auto item: std::views::keys(entropy_present_X_history_X.back()))
            {
                keys2.erase(item);
            }
            auto keys3(keys);
            for ( auto item: std::views::keys(entropy_history_X_history_Y.back()))
            {
                keys3.erase(item);
            }
            auto keys4(keys);
            for ( auto item: std::views::keys(entropy_joint.back()))
            {
                keys4.erase(item);
            }
            auto keys5(keys);
            for ( auto item: std::views::keys(entropy_history_X.back()))
            {
                keys5.erase(item);
            }

            for ( const auto & key : keys )
            {
                if (!entropy_present_X_history_X.back().contains(key))
                {
                    BOOST_LOG_TRIVIAL ( debug ) << std::format("key({}) is missing entropy_present_X_history_X", (key));
                }
                if (!entropy_history_X_history_Y.back().contains(key))
                {
                    BOOST_LOG_TRIVIAL ( debug ) << std::format("key({}) is missing entropy_history_X_history_Y", (key));
                }
                if (!entropy_joint.back().contains(key))
                {
                    BOOST_LOG_TRIVIAL ( debug ) << std::format("key({}) is missing entropy_joint", (key));
                }
                if (!entropy_history_X.back().contains(key))
                {
                    BOOST_LOG_TRIVIAL ( debug ) << std::format("key({}) is missing entropy_history_X", (key));
                }
            }
#endif

            // conditional Renyi entropy
            for (int neighbor_index = 0; neighbor_index < calculator.GetIndices().size(); ++ neighbor_index)
            {
                const auto & neighbor = calculator.GetIndices()[neighbor_index];
                auto sum = std::views::zip_transform ( conditional_information_transfer_calculator, entropy_present_X_history_X[neighbor_index], entropy_history_X_history_Y[neighbor_index], entropy_joint[neighbor_index], entropy_history_X[neighbor_index] );
                decltype(entropy_present_X_history_X) conditional_information_transfer_result(calculator.GetIndices().size());
                for (auto [key, item]: sum)
                {
                    conditional_information_transfer_result[neighbor_index][key] = item;
                }
                RTE_result[ key_generator(conditional_renyi_entropy_label, neighbor) ] = std::move(conditional_information_transfer_result[neighbor_index]);
                RTE_result[ key_generator(renyi_entropy_X_present_history_label, neighbor) ] = std::move(entropy_present_X_history_X[neighbor_index]);
                RTE_result[ key_generator(renyi_entropy_XY_history_label, neighbor) ] = std::move(entropy_history_X_history_Y[neighbor_index]);
                RTE_result[ key_generator(joint_renyi_entropy_label, neighbor) ] = std::move(entropy_joint[neighbor_index]);
                RTE_result[ key_generator(renyi_entropy_X_history_label, neighbor) ] = std::move(entropy_history_X[neighbor_index]);
            }
        } else {

        }
        /*
                    results = {}
            if enhanced_calculation:
                joint_dataset = np.concatenate((data_x_fut, data_x_hist), axis=axis_to_join)
                entropy_present_X_history_X = renyi_entropy(joint_dataset, **kwargs)

                joint_dataset = np.concatenate((data_x_hist, data_y), axis=axis_to_join)
                entropy_history_X_history_Y = renyi_entropy(joint_dataset, **kwargs)

                joint_dataset = np.concatenate(
                    (data_x_fut, data_x_hist, data_y), axis=axis_to_join
                )
                entropy_joint = renyi_entropy(joint_dataset, **kwargs)

                entropy_history_X = renyi_entropy(data_x_hist, **kwargs)

                for alpha in kwargs["alphas"]:
                    try:
                        result = [
                            entropy_present_X_history_X[alpha][index]
                            + entropy_history_X_history_Y[alpha][index]
                            - entropy_joint[alpha][index]
                            - entropy_history_X[alpha][index]
                            for index in range(len(kwargs["indices_to_use"]))
                        ]
                        results[alpha] = result
                    except KeyError as exc:
                        tb = sys.exception().__traceback__
                        print(f"Key {alpha} is missing: {exc.with_traceback(tb)}")
                return results
            else:
                joint_dataset = np.concatenate((data_x_hist, data_y), axis=axis_to_join)

                joint_part = renyi_mutual_entropy(data_x_fut, joint_dataset, **kwargs)
                marginal_part = renyi_mutual_entropy(data_x_fut, data_x_hist, **kwargs)

                results = {}
                for alpha in kwargs["alphas"]:
                    result = [
                        joint_part[alpha][index] - marginal_part[alpha][index]
                        for index in range(len(kwargs["indices_to_use"]))
                    ]
                    results[alpha] = result

                return results*/
    }

    static std::tuple<const Eigen::MatrixXd, const Eigen::MatrixXd, const Eigen::MatrixXd> PrepareDatasetForTransferEntropy ( Eigen::MatrixXd & marginal_solution_1, Eigen::MatrixXd & marginal_solution_2, std::map<std::string, std::any> parameters )
    {
        unsigned int history_first, history_second, skip_first, skip_last, history_x, history_y, time_shift_between_X_Y;
        std::vector<unsigned int> history_index_x, history_index_y, future_index_x;
        unsigned int postselection_y_fut, postselection_z_hist, postselection_y_hist, random_source;
        Eigen::MatrixXd marginal_solution_1_selected, marginal_solution_2_selected;
        if ( parameters.contains ( "history_first" ) ) {
            history_first = std::any_cast<unsigned int> ( parameters["history_first"] );
        } else {
            history_first = 5;
        }

        if ( parameters.contains ( "history_second" ) ) {
            history_second = std::any_cast<unsigned int> ( parameters["history_second"] );
        } else {
            history_second = 5;
        }

        if ( parameters.contains ( "skip_first" ) ) {
            skip_first = std::any_cast<unsigned int> ( parameters["skip_first"] );
        } else {
            skip_first = 0;
        }

        if ( parameters.contains ( "skip_last" ) ) {
            skip_last = std::any_cast<unsigned int> ( parameters["skip_last"] );
        } else {
            skip_last = 0;
        }

        if ( parameters.contains ( "history_index_x" ) ) {
            history_index_x = std::any_cast<std::vector<unsigned int>> ( parameters["history_index_x"] );
            history_x = *max_element ( history_index_x.begin(), history_index_x.end() );
        }

        if ( parameters.contains ( "history_index_y" ) ) {
            history_index_y = std::any_cast<std::vector<unsigned int>> ( parameters["history_index_y"] );
            history_y = *max_element ( history_index_y.begin(), history_index_y.end() );
        }

        if ( parameters.contains ( "time_shift_between_X_Y" ) ) {
            time_shift_between_X_Y = std::any_cast<unsigned long> ( parameters["time_shift_between_X_Y"] );
        } else {
            time_shift_between_X_Y = 1;
        }

        if ( parameters.contains ( "future_index_x" ) ) {
            future_index_x = std::any_cast<std::vector<unsigned int>> ( parameters["future_index_x"] );
        } else {
            future_index_x = std::vector<unsigned int> ( {0U} );
        }

        if ( parameters.contains ( "transpose" ) ) {
            marginal_solution_1 = marginal_solution_1.transpose();
            marginal_solution_2 = marginal_solution_2.transpose();
        }
        unsigned int dimension_timeseries_1 = marginal_solution_1.cols();
        unsigned int dimension_timeseries_2 = marginal_solution_2.cols();

        if ( parameters.contains ( "postselection_y_fut" ) ) {
            postselection_y_fut = std::any_cast<unsigned int> ( parameters["postselection_y_fut"] );
        } else {
            postselection_y_fut = 0;
        }

        if ( parameters.contains ( "postselection_z_hist" ) ) {
            postselection_z_hist = std::any_cast<unsigned int> ( parameters["postselection_z_hist"] );
        } else {
            postselection_z_hist = 0;
        }

        if ( parameters.contains ( "postselection_y_hist" ) ) {
            postselection_y_hist = std::any_cast<unsigned int> ( parameters["postselection_y_hist"] );
        } else {
            postselection_y_hist = 0;
        }

        if ( parameters.contains ( "random_source" ) ) {
            random_source = std::any_cast<bool> ( parameters["random_source"] );
        }

        //std::cout << marginal_solution_1 << std::endl;
        marginal_solution_1_selected = marginal_solution_1.block ( skip_first, 0, marginal_solution_1.rows() - skip_last, dimension_timeseries_1 );
        marginal_solution_2_selected = marginal_solution_2.block ( skip_first, 0, marginal_solution_2.rows() - skip_last, dimension_timeseries_2 );

        const Eigen::MatrixXd matrix;

        parameters["transpose"]= false;
        // additional move in history is there because then actual timeserie is separated
        parameters["history"] = history_x;
        parameters["skip_first"] = ( history_x > history_y ) ? time_shift_between_X_Y : history_y - history_x + time_shift_between_X_Y;
        parameters["skip_last"] = 0U;
        std::vector<unsigned int> join_histories_x_y ( history_index_x.size() + history_index_y.size() );
        join_histories_x_y.insert ( join_histories_x_y.end(), history_index_x.begin(), history_index_x.end() );
        join_histories_x_y.insert ( join_histories_x_y.end(), history_index_y.begin(), history_index_y.end() );
        auto max_join_histories_x_y = std::max_element ( join_histories_x_y.begin(), join_histories_x_y.end() );
        auto max_future_index_x = std::max_element ( future_index_x.begin(), future_index_x.end() );
        auto max_history_index_x =  std::max_element ( history_index_x.begin(), history_index_x.end() );
        auto max_history_index_y =  std::max_element ( history_index_y.begin(), history_index_y.end() );

        if ( parameters.contains ( "history_index_x" ) ) {
            std::vector<unsigned int> indices;
            if ( parameters.contains ( "future_index_x" ) ) {
                for ( unsigned int item : future_index_x ) {
                    indices.push_back ( ( *max_future_index_x ) - item );
                }
                for ( unsigned int item : history_index_x ) {
                    indices.push_back ( ( *max_future_index_x ) + item );
                }
            } else {
                for ( unsigned int item : history_index_x ) {
                    indices.push_back ( 1 + item );
                }
            }
            parameters["select_indices"] = indices;
            parameters["skip_first"] = static_cast<unsigned int> ( ( * max_join_histories_x_y ) + ( * max_future_index_x ) );
        }
        auto samples_marginal_1 = samples_from_arrays ( marginal_solution_1_selected, parameters );

        //std::cout << samples_marginal_1->rows() << " " << samples_marginal_1->cols() << *samples_marginal_1 << std::endl;
        auto separator_row_x = dimension_timeseries_1 * history_index_x.size(); //bug number of original dimension of data
        Eigen::MatrixXd y_fut, y_history, x_hist;
        y_history = ( *samples_marginal_1 ) ( Eigen::seq ( 0, separator_row_x ), Eigen::all );
        y_fut = ( *samples_marginal_1 ) ( Eigen::seq ( separator_row_x, samples_marginal_1->rows() - 1 ), Eigen::all );

        parameters["history"] = history_y;
        parameters["skip_first"] = static_cast<unsigned int> ( ( * max_join_histories_x_y ) + ( * max_future_index_x ) );
        parameters["skip_last"] = 0U;
        if ( parameters.contains ( "history_index_y" ) ) {
            std::vector<unsigned int> indices;
            for ( unsigned int item : history_index_y ) {
                indices.push_back ( ( * max_history_index_y ) - item );
            }
            parameters["select_indices"] = indices;
        }

        auto samples_marginal_2 = samples_from_arrays ( marginal_solution_2_selected, parameters );
        x_hist = *(samples_marginal_2.get());

        if ( parameters.contains ( "transpose" ) ) {
            marginal_solution_1 = marginal_solution_1.transpose();
            marginal_solution_2 = marginal_solution_2.transpose();
        }
        return std::tuple ( y_fut, y_history, x_hist );
    }

    static Eigen::MatrixXd shuffle_sample ( const Eigen::MatrixXd &dataset )
    {
        std::random_device random_device_instance;
        std::seed_seq rng_seed{random_device_instance(), random_device_instance(), random_device_instance(), random_device_instance(), random_device_instance(), random_device_instance(), random_device_instance(), random_device_instance()};
        std::mt19937 random_generator (rng_seed);
        // https://stackoverflow.com/questions/15858569/randomly-permute-rows-columns-of-a-matrix-with-eigen
        //Eigen::MatrixXd final_dataset ( dataset.cols(), dataset.rows() );

        Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm ( dataset.cols() );
        perm.setIdentity();
        std::shuffle ( perm.indices().data(), perm.indices().data() + perm.indices().size(), random_generator );
        // permutate columns
        //dataset.transpose();
        auto result = dataset * perm;
        //std::cout << dataset << "\nres:\n" << result << std::endl;
        return result;
    }

    static Eigen::MatrixXd surrogate_sample ( const Eigen::MatrixXd &dataset )
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(-std::numbers::pi, std::numbers::pi);

        Eigen::MatrixXd result(dataset.rows(), dataset.cols());

        int count = dataset.cols();
        std::vector<double> phases(count / 2);
        for (int i = 1; i < count / 2; ++i)
        {
            double phase = dis(gen);
        }

        for (unsigned int row = 0; row < dataset.rows(); ++ row )
        {
            complexFFTData fftData(count);
            complexFFT fft(&fftData);

            for ( unsigned int column = 0 ; column < count; ++ column )
            {
                c_re(fftData.in[column]) = static_cast<double>(dataset(row, column) );
                c_im(fftData.in[column]) = static_cast<double>(0);
            }

            // forward transform
            fft.fwdTransform();

            // power spectrum
            fftData.compPowerSpec();

            for (int i = 1; i < count / 2; ++ i)
            {
                double phase = phases[i];
                double real_shift = cos(phase);
                double imag_shift = sin(phase);
                c_re(fftData.out[i]) = sqrt( fftData.power_spectrum[i] ) * real_shift;
                c_im(fftData.out[i]) = sqrt( fftData.power_spectrum[i] ) * imag_shift;

                c_re(fftData.out[count - i]) = c_re(fftData.out[i]);
                c_im(fftData.out[count - i]) = - c_im(fftData.out[i]);
            }

            // inverse transform
            fft.invTransform();

            // normalization
            for (int i = 0; i < count; ++ i)
            {
                c_re(fftData.in[i]) /= count;
                c_im(fftData.in[i]) /= count;
            }

            for ( unsigned int column = 0 ; column < count; ++ column )
            {
                result(row, column) = c_re(fftData.in[column]);
                assert(std::fabs(c_im(fftData.in[column])) < 1e-10);

                //std::cout << c_re(fftData.in[column]) << " + i* " << c_im(fftData.in[column]) << std::endl;
            }
        }
        return result;
    }

    static std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> prepare_dataset ( const Eigen::MatrixXd &joint_dataset, bool swap_datasets, bool shuffle_dataset, bool surrogate_dataset, unsigned int selection_1, unsigned int selection_2 )
    {
        auto marginal_solution_1 = joint_dataset ( Eigen::seqN ( 0, selection_1 ), Eigen::all );
        auto marginal_solution_2 = joint_dataset ( Eigen::seqN ( selection_1, selection_2 ), Eigen::all );

        auto * swaped_marginal_solution_1 = &marginal_solution_1;
        auto * swaped_marginal_solution_2 = &marginal_solution_2;

        if ( swap_datasets ) {
            std::swap ( swaped_marginal_solution_1, swaped_marginal_solution_2 );
        }

        Eigen::MatrixXd modified_dataset;
        if ( shuffle_dataset ) {
            modified_dataset = shuffle_sample ( *swaped_marginal_solution_1 );
        }
        else if ( surrogate_dataset ) {
            modified_dataset = surrogate_sample ( *swaped_marginal_solution_1 );
        }
        else
        {
            modified_dataset = * swaped_marginal_solution_1;
        }

        return std::tuple<const Eigen::MatrixXd, const Eigen::MatrixXd> ( modified_dataset, *swaped_marginal_solution_2 );
    }

    static std::shared_ptr<const Eigen::MatrixXd> samples_from_arrays ( Eigen::MatrixXd &data, std::map<std::string, std::any> parameters )
    {
        const unsigned int dimension_of_data = data.cols();
        unsigned int history, allocated_space, skip_first, skip_last;
        std::vector<unsigned int> select_indices;
        std::vector<unsigned int> &range_of_history ( select_indices );
        if ( parameters.contains ( "history" ) ) {
            history = std::any_cast<unsigned int> ( parameters["history"] );
            select_indices = std::views::iota ( 1U, history ) | std::ranges::to<std::vector<unsigned int>>();
        } else {
            allocated_space = 5;
            history = 5;
        }

        if ( parameters.contains ( "select_indices" ) ) {
            select_indices = std::any_cast<std::vector<unsigned int>& > ( parameters["select_indices"] );
            allocated_space = select_indices.size();
            history = *max_element ( select_indices.begin(), select_indices.end() );
            range_of_history = select_indices;
        }

        if ( parameters.contains ( "skip_first" ) ) {
            skip_first = std::any_cast<unsigned int> ( parameters["skip_first"] );
        } else {
            skip_first = history;
        }

        if ( parameters.contains ( "skip_last" ) ) {
            skip_last = std::any_cast<unsigned int> ( parameters["skip_last"] );
        } else {
            skip_last = 0;
        }

        const auto length_of_timeserie = data.rows() - skip_first - skip_last;
        if ( length_of_timeserie <= 0 ) {
            throw std::logic_error ( "Calculation of allocated matrix size gives zero or negative" );
        }
        std::shared_ptr<Eigen::MatrixXd> sampled_dataset ( new Eigen::MatrixXd ( data.cols() * allocated_space, length_of_timeserie ) );
        for ( unsigned int row{skip_first}, count_rows{0} ; row < skip_first + length_of_timeserie; ++row, ++count_rows ) {
            for ( auto [count_history, item_history]: std::views::enumerate ( range_of_history ) ) {
                for ( auto dimension_of_original_dataset: std::views::iota ( 0U, dimension_of_data ) ) {
                    const auto original_dataset_row = row - item_history;
                    const auto inserted_data = data ( original_dataset_row, dimension_of_original_dataset );
                    sampled_dataset->operator() ( count_history * dimension_of_data + dimension_of_original_dataset, count_rows ) = inserted_data;
                }
            }
        }

        return sampled_dataset;
    }

    static std::generator<std::shared_ptr<const Eigen::MatrixXd>> generator_samples_from_arrays ( Eigen::MatrixXd &data, std::map<std::string, std::any> parameters )
    {
        const unsigned int dimension_of_data = data.cols();
        unsigned int history, allocated_space, skip_first, skip_last;
        std::vector<unsigned int> select_indices;
        std::vector<unsigned int> &range_of_history ( select_indices );
        if ( parameters.contains ( "history" ) ) {
            history = std::any_cast<unsigned int> ( parameters["history"] );
            select_indices = std::views::iota ( 1U, history ) | std::ranges::to<std::vector<unsigned int>>();
        } else {
            allocated_space = 5;
            history = 5;
        }

        if ( parameters.contains ( "select_indices" ) ) {
            select_indices = std::any_cast<std::vector<unsigned int>& > ( parameters["select_indices"] );
            allocated_space = select_indices.size();
            history = *max_element ( select_indices.begin(), select_indices.end() );
            range_of_history = select_indices;
        }

        if ( parameters.contains ( "skip_first" ) ) {
            skip_first = std::any_cast<unsigned int> ( parameters["skip_first"] );
        } else {
            skip_first = history;
        }

        if ( parameters.contains ( "skip_last" ) ) {
            skip_last = std::any_cast<unsigned int> ( parameters["skip_last"] );
        } else {
            skip_last = 0;
        }

        const auto length_of_timeserie = data.rows() - skip_first - skip_last;
        std::shared_ptr<Eigen::MatrixXd> sampled_dataset ( new Eigen::MatrixXd ( data.cols() * allocated_space, 1 ) );
        for ( unsigned int row{skip_first}, count_rows{0} ; row < skip_first + length_of_timeserie; ++row, ++count_rows ) {
            for ( auto [count_history, item_history]: std::views::enumerate ( range_of_history ) ) {
                for ( auto dimension_of_original_dataset: std::views::iota ( 0U, dimension_of_data ) ) {
                    const auto original_dataset_row = row - item_history;
                    const auto inserted_data = data ( original_dataset_row, dimension_of_original_dataset );
                    sampled_dataset->operator() ( count_history * dimension_of_data + dimension_of_original_dataset, 0 ) = inserted_data;
                }
            }
            co_yield sampled_dataset;
        }
    }

    bool GetMultithreading ()
    {
        return _multithreading;
    }

    void SetMultithreading ( bool multithreading )
    {
        _multithreading = multithreading;
    }

    const unsigned int GetNumberOfThreads ()
    {
        return _number_of_threads;
    }

    void SetNumberOfThreads ( unsigned int number_of_threads )
    {
        _number_of_threads = number_of_threads;
    }

protected:
    bool _multithreading;

    unsigned int _number_of_threads;

    std::recursive_mutex _result_mutex;

    std::vector<TYPE> _alphas;

    std::vector<unsigned int> _indices;

    std::function<TYPE ( TYPE ) > _logarithm;

    std::function<TYPE ( TYPE ) > _exp;

    std::function<TYPE ( TYPE, TYPE ) > _power;
};

} // namespace renyi_entropy
