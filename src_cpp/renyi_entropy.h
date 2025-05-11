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
#include <numeric>
#include <sstream>
#include <vector>
#include <ranges>
#include <any>
#include <random>
#include <memory>
#include <type_traits>
#include <typeindex>
#include <typeinfo>
#include <generator>

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

#include <eigen3/Eigen/Core>

#include "KDTree.cpp"

#include "ndarray.h"

namespace renyi_entropy {


template<typename Scalar, typename Matrix>
inline std::vector< std::vector<Scalar> > fromEigenMatrix( const Matrix & M ){
	std::vector< std::vector<Scalar> > m;
	m.resize(M.cols(), std::vector<Scalar>(M.rows(), 0));
	for(size_t i = 0; i < m.size(); i++)
		for(size_t j = 0; j < m.front().size(); j++)
			m[i][j] = M(j, i);
	return m;
}
  
template <typename TYPE>
TYPE volume_of_hypersphere(TYPE radius, TYPE metric, unsigned int dimension) {
  const TYPE static_part = pow(2 / metric, dimension - 1);
  const TYPE radius_part = pow(radius, dimension) / dimension;
  const TYPE complete_angular_part = 2 * boost::math::tgamma(1 / metric) *
                                     boost::math::tgamma(1 / metric) /
                                     boost::math::tgamma(2 / metric);

  TYPE noncomplete_angular_part =
      pow(std::numbers::pi_v<double>, (dimension - 2) / 2.);
  for (unsigned int angular_momentum_index = 1;
       angular_momentum_index < dimension - 1; ++angular_momentum_index) {
    noncomplete_angular_part *=
        boost::math::tgamma(0.5 + (dimension - (angular_momentum_index + 1)) /
                                      metric) /
        boost::math::tgamma(
            (dimension + metric - (angular_momentum_index + 1)) / metric);
  }

  return static_part * radius_part * complete_angular_part *
         noncomplete_angular_part;
}

template <typename TYPE> 
class renyi_entropy {

public:
  renyi_entropy() {}

  ~renyi_entropy() {}

  KDTree<TYPE> create_KDtree(std::vector<std::vector<TYPE>> &dataset) {
    KDTree<TYPE> tree(dataset);
    return tree;
  }

  void calculation() {
    std::vector<double> dataset_x;
    std::vector<std::vector<double>> distances;
    for (auto const alpha : _alphas) {
      if (alpha == 1.0) {
        auto result =
            entropy_sum_Shannon_LeonenkoProzanto(dataset_x, distances);
      } else {
        auto entropy_sum =
            entropy_sum_generic_LeonenkoProzanto(dataset_x, distances, alpha);
      }
    }
  }

  void entropy_sum_Renyi_LeonenkoProzanto(
      std::map<std::tuple<unsigned int, TYPE>, TYPE> &results,
      int dimension_of_data, const std::vector<std::vector<TYPE>> &distances,
      const TYPE alpha, bool log_calculation) {
    std::map<unsigned int, TYPE> entropy;
    TYPE one_minus_alpha = TYPE(1.0) - alpha;

    for (const unsigned int use_index : GetIndices()) {
      auto const number_of_data =
          std::accumulate(distances.begin(), distances.end(), 0,
                          [=](unsigned int count, auto const &element) {
                            if (element[use_index] > 0.0) {
                              return ++count;
                            } else {
                              return count;
                            }
                          });

      auto maximum_distance_element =
          std::max_element(distances.cbegin(), distances.cend(),
                           [use_index](const auto &first, const auto &second) {
                             auto first_data = first[use_index];
                             auto second_data = second[use_index];
                             return first_data > second_data ? false : true;
                           });
      auto maximum_distance = (*maximum_distance_element)[use_index];
      TYPE exponent = dimension_of_data * one_minus_alpha;

      auto sum_of_power_of_distances = std::accumulate(
          distances.cbegin(), distances.cend(), static_cast<TYPE>(0),
          [use_index, maximum_distance, exponent, this](TYPE count,
                                                        auto &element) {
            auto power_base = element[use_index];
            if (power_base > 0.0) {
              return count + _power(power_base / maximum_distance, exponent);
            } else {
              return count;
            }
          });

      auto index_exponent = one_minus_alpha + TYPE(use_index);
      if (log_calculation) {
        auto multiplicator =
            boost::math::lgamma(use_index) -
            boost::math::lgamma(index_exponent) -
            one_minus_alpha * boost::math::lgamma(dimension_of_data / 2. + 1.) +
            exponent / 2. * _logarithm(TYPE(std::numbers::pi_v<double>)) +
            exponent * _logarithm(maximum_distance) -
            _logarithm(number_of_data) +
            one_minus_alpha * _logarithm(number_of_data - 1);

        results[std::make_tuple(use_index, alpha)] =
            _exp(multiplicator + _logarithm(sum_of_power_of_distances));
      } else {
        results[std::make_tuple(use_index, alpha)] =
            sum_of_power_of_distances / number_of_data *
            pow(number_of_data - 1, one_minus_alpha) *
            boost::math::tgamma(use_index) /
            boost::math::tgamma(index_exponent) *
            _power(TYPE(std::numbers::pi_v<double>), exponent / 2) /
            _power(boost::math::tgamma(dimension_of_data / 2. + 1.),
                   one_minus_alpha);
      }
    }
  }

  void entropy_sum_Shannon_LeonenkoProzanto(
      std::map<std::tuple<unsigned int, TYPE>, TYPE> &results,
      int dimension_of_data, std::vector<std::vector<TYPE>> &distances,
      bool log_calculation) {
    for (const unsigned int use_index : GetIndices()) {
      auto const number_of_data =
          std::accumulate(distances.cbegin(), distances.cend(), 0,
                          [=](unsigned int count, auto const &element) {
                            if (element[use_index] > 0.0) {
                              return ++count;
                            } else {
                              return count;
                            }
                          });

      if (log_calculation) {
        const auto V_m =
            _power(std::numbers::pi_v<double>, dimension_of_data / 2.0) /
            boost::math::tgamma(dimension_of_data / 2.0 + 1);
        const auto digamma = boost::math::digamma(use_index);
        const auto multiplicator = (number_of_data - 1) * _exp(-digamma);
        auto const sum_of_weighted_distances = std::accumulate(
            distances.begin(), distances.end(), static_cast<TYPE>(0),
            [&](TYPE count, auto const &element) {
              const auto base = element[use_index];
              if (base > 0.0) {
                return count + _logarithm(multiplicator * V_m *
                                          +_power(base, dimension_of_data));
              } else {
                return count;
              }
            });
        results[std::make_tuple(use_index, static_cast<TYPE>(1.0))] =
            sum_of_weighted_distances / number_of_data;
      } else {
        auto const sum_of_log_distances = std::accumulate(
            distances.begin(), distances.end(), static_cast<TYPE>(0),
            [use_index, this](TYPE count, auto const &element) {
              const auto log_base = element[use_index];
              if (log_base > 0.0) {
                return count + _logarithm(log_base);
              } else {
                return count;
              }
            });

        auto const addition_to_entropy =
            dimension_of_data * sum_of_log_distances / number_of_data;
        auto digamma = boost::math::digamma(use_index);
        auto V_m = _power(std::numbers::pi_v<double>, dimension_of_data / 2.0) /
                   boost::math::tgamma(dimension_of_data / 2.0 + 1);
        auto argument_log = V_m * (number_of_data - 1);
        auto log_volume = _logarithm(argument_log);
        results[std::make_tuple(use_index, static_cast<TYPE>(1.0))] =
            addition_to_entropy + log_volume - digamma;
      }
    }
  }

  std::map<std::tuple<unsigned int, TYPE>, TYPE>
  renyi_entropy_LeonenkoProzanto(std::vector<std::vector<TYPE>> &dataset,
                                 double metric = 2) {
    auto nearest_iterator =
        std::max_element(GetIndices().begin(), GetIndices().end());
    // calculation how large of index of distance needs to be calculated
    // + 1 for skipping 0-neighbor
    auto nearest = (*nearest_iterator) + 1;
    std::map<std::tuple<unsigned int, TYPE>, TYPE> results;
    int dimension_of_data = dataset[0].size();
    auto start_calculation = std::chrono::high_resolution_clock::now();
    auto kdtree = create_KDtree(dataset);
    auto end_tree_construction_calculation =
        std::chrono::high_resolution_clock::now();

    // calculate distances
    std::vector<std::vector<TYPE>> distances;
    calculate_distances_from_each_point(kdtree, dataset, nearest, distances,
                                        metric);

    // remove kdtree to same memory
    auto end_distance_calculation = std::chrono::high_resolution_clock::now();

    // calculate renyi entropy
    auto log_calculation = true;
    auto one = static_cast<TYPE>(1);
    for (const auto alpha : GetAlphas()) {
      if (alpha == one) {
        entropy_sum_Shannon_LeonenkoProzanto(results, dimension_of_data,
                                             distances, log_calculation);
      } else {
        entropy_sum_Renyi_LeonenkoProzanto(results, dimension_of_data,
                                           distances, alpha, log_calculation);
      }
    }

    // calculation of entropy for Renyi case
    // Shannon case already holds entropy
    for (auto &item : results) {
      auto alpha = std::get<1>(item.first);
      if (alpha != one) {
        auto entropy = _logarithm(item.second) / (1 - alpha);
        item.second = entropy;
      }
    }

    auto end_entropy_calculation = std::chrono::high_resolution_clock::now();

    auto elapsed_tree_construction =
        end_tree_construction_calculation - start_calculation;
    auto elapsed_distance_calculation =
        end_distance_calculation - end_tree_construction_calculation;
    auto elapsed_entropy_calculation =
        end_entropy_calculation - end_distance_calculation;

    auto milliseconds_tree_construction =
        std::chrono::duration_cast<std::chrono::milliseconds>(
            elapsed_tree_construction)
            .count();
    auto milliseconds_distance_calculation =
        std::chrono::duration_cast<std::chrono::milliseconds>(
            elapsed_distance_calculation)
            .count();
    auto milliseconds_entropy_calculation =
        std::chrono::duration_cast<std::chrono::milliseconds>(
            elapsed_entropy_calculation)
            .count();

    BOOST_LOG_TRIVIAL(trace)
        << std::format("Tree constructed in {:.5f} seconds",
                       milliseconds_tree_construction / 1000.0);
    BOOST_LOG_TRIVIAL(trace)
        << std::format("Distances calculated in {:.5f} seconds",
                       milliseconds_distance_calculation / 1000.0);
    BOOST_LOG_TRIVIAL(trace)
        << std::format("Entropy calculated in {:.5f} seconds",
                       milliseconds_entropy_calculation / 1000.0);
    return results;
  }

  void entropy_sum_Renyi_metric(
      std::map<std::tuple<unsigned int, TYPE>, TYPE> &results,
      int dimension_of_data, const std::vector<std::vector<TYPE>> &distances,
      const TYPE alpha, bool log_calculation, TYPE metric) {
    std::map<unsigned int, TYPE> entropy;
    TYPE one_minus_alpha = TYPE(1.0) - alpha;

    for (const unsigned int use_index : GetIndices()) {
      auto const number_of_data =
          std::accumulate(distances.begin(), distances.end(), 0,
                          [=](unsigned int count, auto const &element) {
                            if (element[use_index] > 0.0) {
                              return ++count;
                            } else {
                              return count;
                            }
                          });

      const auto maximum_distance_element =
          std::max_element(distances.cbegin(), distances.cend(),
                           [use_index](const auto &first, const auto &second) {
                             auto first_data = first[use_index];
                             auto second_data = second[use_index];
                             return first_data > second_data ? false : true;
                           });
      const auto maximum_distance = (*maximum_distance_element)[use_index];
      const TYPE exponent = dimension_of_data * one_minus_alpha;

      const auto sum_of_power_of_distances = std::accumulate(
          distances.cbegin(), distances.cend(), static_cast<TYPE>(0),
          [use_index, maximum_distance, exponent, this](TYPE count,
                                                        auto &element) {
            auto power_base = element[use_index];
            if (power_base > 0.0) {
              return count + _power(power_base / maximum_distance, exponent);
            } else {
              return count;
            }
          });

      const auto dimensional_multiplicator =
          volume_of_hypersphere(TYPE(1), metric, dimension_of_data);
      const auto index_exponent = one_minus_alpha + TYPE(use_index);
      if (log_calculation) {
        const auto multiplicator = one_minus_alpha * _logarithm(dimensional_multiplicator) +
            boost::math::lgamma(use_index) -
            boost::math::lgamma(index_exponent) +
            exponent * _logarithm(maximum_distance) -
            _logarithm(number_of_data) +
            one_minus_alpha * _logarithm(number_of_data - 1);
        results[std::make_tuple(use_index, alpha)] = _exp(multiplicator + _logarithm(sum_of_power_of_distances));
      } else {
        results[std::make_tuple(use_index, alpha)] =
            sum_of_power_of_distances / number_of_data *
            pow(number_of_data - 1, one_minus_alpha) *
            boost::math::tgamma(use_index) /
            boost::math::tgamma(index_exponent) * _power( dimensional_multiplicator, one_minus_alpha );
      }
    }
  }

  void entropy_sum_Shannon_metric(
      std::map<std::tuple<unsigned int, TYPE>, TYPE> &results,
      int dimension_of_data, std::vector<std::vector<TYPE>> &distances,
      bool log_calculation, TYPE metric) {
    for (const unsigned int use_index : GetIndices()) {
      auto const number_of_data =
          std::accumulate(distances.cbegin(), distances.cend(), 0,
                          [=](unsigned int count, auto const &element) {
                            if (element[use_index] > 0.0) {
                              return ++count;
                            } else {
                              return count;
                            }
                          });

      const auto V_m =
          volume_of_hypersphere(TYPE(1), metric, dimension_of_data);
      if (log_calculation) {
        const auto digamma = boost::math::digamma(use_index);
        const auto multiplicator = (number_of_data - 1) * _exp(-digamma);
        auto const sum_of_weighted_distances = std::accumulate(
            distances.begin(), distances.end(), static_cast<TYPE>(0),
            [&](TYPE count, auto const &element) {
              const auto base = element[use_index];
              if (base > 0.0) {
                return count + _logarithm(multiplicator * V_m *
                                          +_power(base, dimension_of_data));
              } else {
                return count;
              }
            });
        results[std::make_tuple(use_index, static_cast<TYPE>(1.0))] =
            sum_of_weighted_distances / number_of_data;
      } else {
        auto const sum_of_log_distances = std::accumulate(
            distances.begin(), distances.end(), static_cast<TYPE>(0),
            [use_index, this](TYPE count, auto const &element) {
              const auto log_base = element[use_index];
              if (log_base > 0.0) {
                return count + _logarithm(log_base);
              } else {
                return count;
              }
            });

        const auto addition_to_entropy =
            dimension_of_data * sum_of_log_distances / number_of_data;
        const auto digamma = boost::math::digamma(use_index);
        auto argument_log = V_m * (number_of_data - 1);
        auto log_volume = _logarithm(argument_log);
        results[std::make_tuple(use_index, static_cast<TYPE>(1.0))] =
            addition_to_entropy + log_volume - digamma;
      }
    }
  }

  std::map<std::tuple<unsigned int, TYPE>, TYPE>
  renyi_entropy_metric(std::vector<std::vector<TYPE>> &dataset,
                       double metric = 2) {
    auto nearest_iterator =
        std::max_element(GetIndices().begin(), GetIndices().end());
    // calculation how large of index of distance needs to be calculated
    // + 1 for skipping 0-neighbor
    auto nearest = (*nearest_iterator) + 1;
    std::map<std::tuple<unsigned int, TYPE>, TYPE> results;
    const int dimension_of_data = dataset[0].size();
    auto start_calculation = std::chrono::high_resolution_clock::now();
    auto kdtree = create_KDtree(dataset);
    auto end_tree_construction_calculation =
        std::chrono::high_resolution_clock::now();

    // calculate distances
    std::vector<std::vector<TYPE>> distances;
    calculate_distances_from_each_point(kdtree, dataset, nearest, distances,
                                        metric);

    // remove kdtree to same memory
    auto end_distance_calculation = std::chrono::high_resolution_clock::now();

    // calculate renyi entropy
    auto log_calculation = true;
    auto one = static_cast<TYPE>(1);
    for (const auto alpha : GetAlphas()) {
      if (alpha == one) {
        entropy_sum_Shannon_metric(results, dimension_of_data, distances,
                                   log_calculation, metric);
      } else {
        entropy_sum_Renyi_metric(results, dimension_of_data, distances, alpha,
                                 log_calculation, metric);
      }
    }

    // calculation of entropy for Renyi case
    // Shannon case already holds entropy
    const auto fixed_compesation = -0.001921 * dimension_of_data + 1.02726;
    
    const auto exponent1 = (0.167 / dimension_of_data - 0.181);
    const auto multiplication_factor1 = 0.95;
    const auto threshold1 = exp( log(fixed_compesation) - log(multiplication_factor1) / exponent1 );

    const auto exponent2 = (-0.0015 * dimension_of_data + 0.00067);
    const auto multiplication_factor2 = 0.99;
    const auto threshold2 = exp( log(fixed_compesation) - log(multiplication_factor2) / exponent2 );
    
    for (auto &item : results) {
      auto alpha = std::get<1>(item.first);
        const auto compensation1 = ( alpha < threshold1 ) ? (multiplication_factor1 * pow(alpha, exponent1)) : 1;        
        const auto compensation2 = ( alpha > threshold2 ) ? (multiplication_factor2 * pow(alpha, exponent2)) : 1;

      if (alpha != one) {        
        auto entropy = _logarithm(item.second) / (1 - alpha) * compensation1 * compensation2 * fixed_compesation;
        item.second = entropy;
      }
      else
      {
        auto entropy = item.second * compensation1 * compensation2 * fixed_compesation;
        item.second = entropy;
      }
    }

    auto end_entropy_calculation = std::chrono::high_resolution_clock::now();

    auto elapsed_tree_construction =
        end_tree_construction_calculation - start_calculation;
    auto elapsed_distance_calculation =
        end_distance_calculation - end_tree_construction_calculation;
    auto elapsed_entropy_calculation =
        end_entropy_calculation - end_distance_calculation;

    auto milliseconds_tree_construction =
        std::chrono::duration_cast<std::chrono::milliseconds>(
            elapsed_tree_construction)
            .count();
    auto milliseconds_distance_calculation =
        std::chrono::duration_cast<std::chrono::milliseconds>(
            elapsed_distance_calculation)
            .count();
    auto milliseconds_entropy_calculation =
        std::chrono::duration_cast<std::chrono::milliseconds>(
            elapsed_entropy_calculation)
            .count();

    BOOST_LOG_TRIVIAL(trace)
        << std::format("Tree constructed in {:.5f} seconds",
                       milliseconds_tree_construction / 1000.0);
    BOOST_LOG_TRIVIAL(trace)
        << std::format("Distances calculated in {:.5f} seconds",
                       milliseconds_distance_calculation / 1000.0);
    BOOST_LOG_TRIVIAL(trace)
        << std::format("Entropy calculated in {:.5f} seconds",
                       milliseconds_entropy_calculation / 1000.0);
    return results;
  }

  inline void calculate_distances_from_each_point(
      KDTree<TYPE> &kdtree, std::vector<std::vector<TYPE>> &dataset,
      unsigned int nearest, std::vector<std::vector<TYPE>> &distances,
      const TYPE metric) {
    // calculate distances
    for (const auto &point : dataset) {
      auto points = kdtree.nearest_points(point, nearest);
      std::vector<TYPE> distance_vector;
      for (int i = 0; i < nearest; ++i) {
        distance_vector.push_back(distance(point, points[i], metric));
      }
      distances.push_back(distance_vector);
    }
  }

  std::vector<unsigned int> &GetIndices() { return _indices; }

  const std::vector<TYPE> &GetAlphas() { return _alphas; }

  static TYPE distance(const std::vector<TYPE> &a, const std::vector<TYPE> &b,
                       TYPE power = 2) {
    if (a.size() == b.size()) {
      TYPE distance = 0;
      for (int i = 0; i < a.size(); ++i) {
        distance += pow(abs(a[i] - b[i]), power);
      }
      return pow(distance, 1 / power);
    }
    return NAN;
  }

  void SetAlpha(const std::vector<TYPE> &alpha) { _alphas = alpha; }

  void SetIndices(const std::vector<unsigned int> &indices) {
    _indices = indices;
  }

  void SetPower(std::function<TYPE(TYPE, TYPE)> power) { _power = power; }

  void SetExp(std::function<TYPE(TYPE)> exp) { _exp = exp; }

  void SetLog(std::function<TYPE(TYPE)> log) { _logarithm = log; }

  static void SaveRenyiEntropy(std::map<std::tuple<unsigned int, TYPE>, TYPE> &result, const std::string &file)
  {
    std::stringstream ss;
    boost::filesystem::path myFile =
        boost::filesystem::current_path() / file;
    boost::filesystem::ofstream ofs ( myFile );
    boost::archive::binary_oarchive oarch ( ofs );
    oarch << result;
    std::cout << ss.str();
  }
  
  static void LoadRenyiEntropy(std::map<std::tuple<unsigned int, TYPE>, TYPE> &result, const std::string &file)
  {
    boost::filesystem::path myFile =
        boost::filesystem::current_path() / file;
    boost::filesystem::ifstream ifs(myFile);
    boost::archive::binary_iarchive iarch(ifs);
    iarch >> result;
  }
  
  static int renyi_conditional_information_transfer(const Eigen::MatrixXd &y_future, const Eigen::MatrixXd &y_history, const Eigen::MatrixXd &z_history, std::map<std::string, std::any> parameters)
  {
    
  }
  
  static std::tuple<const Eigen::MatrixXd, const Eigen::MatrixXd, const Eigen::MatrixXd> PrepareDatasetForTransferEntropy ( Eigen::MatrixXd & marginal_solution_1, Eigen::MatrixXd & marginal_solution_2, std::map<std::string, std::any> parameters)
  {
      unsigned int history_first, history_second, skip_first, skip_last, history_x, history_y, time_shift_between_X_Y;
      std::vector<unsigned int> history_index_x, history_index_y, future_index_x;
      bool postselection_y_fut, postselection_z_hist, postselection_y_hist, random_source;
      Eigen::MatrixXd marginal_solution_1_selected, marginal_solution_2_selected;
      if ( parameters.contains("history_first") )
      {
          history_first = std::any_cast<unsigned int>(parameters["history_first"]);
      }
      else
      {
          history_first = 5;
      }

      if ( parameters.contains("history_second") )
      {
          history_second = std::any_cast<unsigned int>(parameters["history_second"]);
      }
      else
      {
          history_second = 5;
      }

      if (parameters.contains("skip_first"))
      {
        skip_first = std::any_cast<unsigned int>(parameters["skip_first"]);
      }
      else
      {
          skip_first = 0;
      }

      if (parameters.contains("skip_last"))
      {
        skip_last = std::any_cast<unsigned int>(parameters["skip_last"]);
      }
      else
      {
          skip_last = 0;
      }

      if (parameters.contains("history_index_x"))
      {
        history_index_x = std::any_cast<std::vector<unsigned int>>(parameters["history_index_x"]);
        history_x = *max_element(history_index_x.begin(), history_index_x.end());
      }

      if (parameters.contains("history_index_y"))
      {
        history_index_y = std::any_cast<std::vector<unsigned int>>(parameters["history_index_y"]);
        history_y = *max_element(history_index_y.begin(), history_index_y.end());
      }

      if (parameters.contains("time_shift_between_X_Y"))
      {
        time_shift_between_X_Y = std::any_cast<unsigned int>(parameters["time_shift_between_X_Y"]);
      }
      else
      {
          time_shift_between_X_Y = 1;
      }

      if (parameters.contains("future_index_x"))
      {
        future_index_x = std::any_cast<std::vector<unsigned int>>(parameters["future_index_x"]);
      }
      else
      {
          future_index_x = std::vector<unsigned int>({0U});
      }

      if (parameters.contains("transpose"))
      {
          marginal_solution_1 = marginal_solution_1.transpose();
          marginal_solution_2 = marginal_solution_2.transpose();
      }

      if (parameters.contains("postselection_y_fut"))
      {
          postselection_y_fut = std::any_cast<bool>(parameters["postselection_y_fut"]);
      }
      else
      {
          postselection_y_fut = false;
      }

      if (parameters.contains("postselection_z_hist"))
      {
          postselection_z_hist = std::any_cast<bool>(parameters["postselection_z_hist"]);
      }
      else
      {
          postselection_z_hist = false;
      }

      if (parameters.contains("postselection_y_hist"))
      {
        postselection_y_hist = std::any_cast<bool>(parameters["postselection_y_hist"]);
      }
      else
      {
          postselection_y_hist = false;
      }

      if (parameters.contains("random_source"))
      {
          random_source = std::any_cast<bool>(parameters["random_source"]);
      }

      marginal_solution_1_selected = marginal_solution_1.block(0, skip_first, marginal_solution_1.cols(), marginal_solution_1.rows() - skip_last);
      marginal_solution_2_selected = marginal_solution_2.block(0, skip_first, marginal_solution_2.cols(), marginal_solution_2.rows() - skip_last);

      const Eigen::MatrixXd matrix;
      return std::tuple(matrix, matrix, matrix);
  }

    static Eigen::MatrixXd shuffle_sample ( const Eigen::MatrixXd &dataset )
    {
      // https://stackoverflow.com/questions/15858569/randomly-permute-rows-columns-of-a-matrix-with-eigen
        Eigen::MatrixXd final_dataset ( dataset.cols(), dataset.rows() );
        
        Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(dataset.rows());
        perm.setIdentity();
        std::random_shuffle(perm.indices().data(), perm.indices().data() + perm.indices().size());
        auto result = perm * dataset;
        return result;
    }

    static std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> prepare_dataset( const Eigen::MatrixXd &joint_dataset, bool swap_datasets, bool shuffle_dataset, unsigned int selection_1, unsigned int selection_2)
    {
        auto marginal_solution_1 = joint_dataset( Eigen::all, Eigen::seqN(0, selection_1) );
        auto marginal_solution_2 = joint_dataset( Eigen::all, Eigen::seqN(selection_1, selection_2) );
        
        auto * swaped_marginal_solution_1 = &marginal_solution_1;
        auto * swaped_marginal_solution_2 = &marginal_solution_2;
        
        if (swap_datasets)
        {
            std::swap( swaped_marginal_solution_1, swaped_marginal_solution_2);
        }
        
        Eigen::MatrixXd shuffled_dataset;
        if (shuffle_dataset)
        {
            shuffled_dataset = shuffle_sample( *swaped_marginal_solution_1 );
        }
        else
        {
            shuffled_dataset = * swaped_marginal_solution_1;
        }
        
        return std::tuple<const Eigen::MatrixXd, const Eigen::MatrixXd>(shuffled_dataset, *swaped_marginal_solution_2);
    }
  
    static std::shared_ptr<const Eigen::MatrixXd> samples_from_arrays( Eigen::MatrixXd &data, std::map<std::string, std::any> parameters)
    {
        const unsigned int dimension_of_data = data.cols();
        unsigned int history, allocated_space, skip_first, skip_last;
        std::vector<unsigned int> select_indices;
        std::vector<unsigned int> &range_of_history (select_indices);
        if (parameters.contains("history"))
        {
            history = std::any_cast<unsigned int>(parameters["history"]);
            select_indices = std::views::iota(1U, history) | std::ranges::to<std::vector<unsigned int>>();
        }
        else
        {
            allocated_space = 5;
            history = 5;
        }

        if (parameters.contains("select_indices"))
        {
            select_indices = std::any_cast<std::vector<unsigned int>& >( parameters["select_indices"] );
            allocated_space = select_indices.size();
            history = *max_element(select_indices.begin(), select_indices.end());
            range_of_history = select_indices;
        }

        if (parameters.contains("skip_first"))
        {
          skip_first = std::any_cast<unsigned int>(parameters["skip_first"]);
        }
        else
        {
            skip_first = history;
        }

        if (parameters.contains("skip_last"))
        {
          skip_last = std::any_cast<unsigned int>(parameters["skip_last"]);
        }
        else
        {
            skip_last = 0;
        }

        const auto length_of_timeserie = data.rows() - skip_first - skip_last;
        std::shared_ptr<Eigen::MatrixXd> sampled_dataset ( new Eigen::MatrixXd( data.cols() * allocated_space, length_of_timeserie) );
        for ( unsigned int row{skip_first}, count_rows{0} ; row < skip_first + length_of_timeserie; ++row, ++count_rows )
        {
            for( auto [count_history , item_history]: std::views::enumerate(range_of_history))
            {
                for (auto dimension_of_original_dataset: std::views::iota(0U, dimension_of_data))
                {
                    const auto original_dataset_row = row - item_history;
                    const auto inserted_data = data(original_dataset_row, dimension_of_original_dataset);
                    sampled_dataset->operator()(count_history * dimension_of_data + dimension_of_original_dataset, count_rows) = inserted_data;
                }
            }
        }
      
        return sampled_dataset;
    }

    static std::generator<std::shared_ptr<const Eigen::MatrixXd>> generator_samples_from_arrays( Eigen::MatrixXd &data, std::map<std::string, std::any> parameters)
    {
        const unsigned int dimension_of_data = data.cols();
        unsigned int history, allocated_space, skip_first, skip_last;
        std::vector<unsigned int> select_indices;
        std::vector<unsigned int> &range_of_history (select_indices);
        if (parameters.contains("history"))
        {
            history = std::any_cast<unsigned int>(parameters["history"]);
            select_indices = std::views::iota(1U, history) | std::ranges::to<std::vector<unsigned int>>();
        }
        else
        {
            allocated_space = 5;
            history = 5;
        }

        if (parameters.contains("select_indices"))
        {
            select_indices = std::any_cast<std::vector<unsigned int>& >( parameters["select_indices"] );
            allocated_space = select_indices.size();
            history = *max_element(select_indices.begin(), select_indices.end());
            range_of_history = select_indices;
        }

        if (parameters.contains("skip_first"))
        {
          skip_first = std::any_cast<unsigned int>(parameters["skip_first"]);
        }
        else
        {
            skip_first = history;
        }

        if (parameters.contains("skip_last"))
        {
          skip_last = std::any_cast<unsigned int>(parameters["skip_last"]);
        }
        else
        {
            skip_last = 0;
        }

        const auto length_of_timeserie = data.rows() - skip_first - skip_last;
        std::shared_ptr<Eigen::MatrixXd> sampled_dataset ( new Eigen::MatrixXd( data.cols() * allocated_space, 1) );
        for ( unsigned int row{skip_first}, count_rows{0} ; row < skip_first + length_of_timeserie; ++row, ++count_rows )
        {
            for( auto [count_history , item_history]: std::views::enumerate(range_of_history))
            {
                for (auto dimension_of_original_dataset: std::views::iota(0U, dimension_of_data))
                {
                    const auto original_dataset_row = row - item_history;
                    const auto inserted_data = data(original_dataset_row, dimension_of_original_dataset);
                    sampled_dataset->operator()(count_history * dimension_of_data + dimension_of_original_dataset, 0) = inserted_data;
                    
                }
            }
            co_yield sampled_dataset;
        }
    }
    
protected:
    std::vector<TYPE> _alphas;

    std::vector<unsigned int> _indices;

    std::function<TYPE(TYPE)> _logarithm;

    std::function<TYPE(TYPE)> _exp;

    std::function<TYPE(TYPE, TYPE)> _power;
};

} // namespace renyi_entropy
