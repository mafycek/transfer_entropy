//
// Created by hynek on 9.4.25.
//

#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <ranges>
#include <string>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/multiprecision/number.hpp>
#include <boost/program_options.hpp>

#include <gtest/gtest.h>

#include "renyi_entropy.h"

int main(int argc, char *argv[]) {
  std::vector<std::string> methods{"LeonenkoProzanto",
                                   "LeonenkoWithGeneralizedMetric"};
  // Declare the supported options.
  boost::program_options::options_description desc("Allowed options");
  std::string version("1.0.0");
  desc.add_options()("help", "Gaussian distribution analyzed by Renyi entropy")(
      "version", "Version of the program")(
      "dimension",
      boost::program_options::value<unsigned int>()->default_value(2),
      "Dimension")(
      "metric", boost::program_options::value<double>()->default_value(2),
      "Metric")("method",
                boost::program_options::value<std::string>()->default_value(
                    "LeonenkoProzanto"),
                "Method");

  boost::program_options::variables_map vm;
  boost::program_options::store(
      boost::program_options::parse_command_line(argc, argv, desc), vm);
  boost::program_options::notify(vm);

  // show help
  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 1;
  }

  // show version of the program
  if (vm.count("version")) {
    std::cout << version << std::endl;
    return 1;
  }

  const double mean_gaussion_distribution = 0;
  const double sigma_gaussion_distribution = 1;
  const unsigned int dimension = vm["dimension"].as<unsigned int>();
  const double metric = vm["metric"].as<double>();
  const std::string method = vm["method"].as<std::string>();

  const double delta_alpha = 0.01;
  auto number_samples = 50000;
  std::random_device random{};
  std::mt19937 generator{random()};
  std::normal_distribution distribution{mean_gaussion_distribution,
                                        sigma_gaussion_distribution};

  auto random_double = [&distribution, &generator] {
    return distribution(generator);
  };

  std::vector<std::vector<double>> dataset;
  for (auto n{number_samples}; n; --n) {
    std::vector<double> sample;
    for (int i :
         std::ranges::iota_view<unsigned int, unsigned int>(0, dimension)) {
      sample.push_back(random_double());
    }
    dataset.push_back(sample);
  }

  std::vector<double> alphas;
  for (double alpha = 0.05; alpha < 1; alpha += delta_alpha) {
    alphas.push_back(alpha);
  }
  alphas.push_back(0.999);
  alphas.push_back(1);
  alphas.push_back(1.001);
  for (double alpha = 1 + delta_alpha; alpha < 5; alpha += delta_alpha) {
    alphas.push_back(alpha);
  }

  renyi_entropy::renyi_entropy<double> calculator;
  std::vector<unsigned int> indices{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

  calculator.SetAlpha(alphas);
  calculator.SetIndices(indices);
  calculator.SetExp(exp);
  calculator.SetLog(log);
  calculator.SetPower(pow);

  std::map<std::tuple<unsigned int, double>, double> result;
  if (method == methods[0]) {
    result = calculator.renyi_entropy_LeonenkoProzanto(dataset, 2);
  } else if (method == methods[1]) {
    result = calculator.renyi_entropy_metric(dataset, metric);
  }

  std::stringstream ss;
  boost::filesystem::path myFile =
      boost::filesystem::current_path() / "myfile.dat";
  boost::filesystem::ofstream ofs(myFile);
  boost::archive::binary_oarchive oarch(ofs);
  oarch << result;
  std::cout << ss.str();

  Eigen::MatrixXd sigma = Eigen::MatrixXd ::Identity(dimension, dimension) *
                          sigma_gaussion_distribution *
                          sigma_gaussion_distribution;
  // Eigen::MatrixXd
  // sigma{{sigma_gaussion_distribution*sigma_gaussion_distribution, 0}, {0,
  // sigma_gaussion_distribution*sigma_gaussion_distribution}};
  for (const auto &item_alpha : calculator.GetAlphas()) {
    std::cout << item_alpha << " ";
    for (const auto &item_index : calculator.GetIndices()) {
      auto result_entropy = result[std::make_tuple(item_index, item_alpha)];
      std::cout << result_entropy << " ";
    }
    std::cout << renyi_entropy::Renyi_entropy_normal_distribution(
                     static_cast<double>(item_alpha), sigma)
              << " ";
    for (const auto &item_index : calculator.GetIndices()) {
      auto result_relative_entropy =
          renyi_entropy::Renyi_entropy_normal_distribution(
              static_cast<double>(item_alpha), sigma) /
          result[std::make_tuple(item_index, item_alpha)];
      std::cout << result_relative_entropy << " ";
    }
    std::cout << std::endl;
  }
}
