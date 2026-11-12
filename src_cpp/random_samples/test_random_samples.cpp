
#include <algorithm>
#include <format>
#include <iostream>
#include <utils/enumerate.hpp>

#include <gtest/gtest.h>

#include "random_samples.h"

TEST(RandomSamples, SimpleNormalDistribution) {
  Eigen::MatrixXd dataset;

  random_samples::samples_normal_distribution_uncorrelated(dataset, 0, 1, 20,
                                                           5);
  std::cout << dataset.transpose() << std::endl;
}

TEST(RandomSamples, CorrelatedNormalDistribution) {
  Eigen::MatrixXd dataset2;
  Eigen::MatrixXd Sigma(5, 5);
  Eigen::VectorXd mean_gaussion_distribution =
      Eigen::VectorXd::Zero(Sigma.cols());
  Sigma << 1, 0.5, 0, 0, 0.5, 0.5, 1, 0.5, 0, 0, 0, 0.5, 1, 0.5, 0, 0, 0, 0.5,
      1, 0.5, 0.5, 0, 0, 0.5, 1;

  random_samples::samples_normal_distribution_correlated(
      dataset2, mean_gaussion_distribution, Sigma, 1000);
  std::cout << dataset2.transpose() << std::endl;
}

TEST(RandomSamples, SimpleStudentTlDistribution) {
  double degrees_of_freedom = 2;
  double sigma_gaussion_distribution = 1;
  unsigned int dimension = 5;
  Eigen::MatrixXd dataset;
  Eigen::MatrixXd sigma = Eigen::MatrixXd::Identity(dimension, dimension) *
                          sigma_gaussion_distribution *
                          sigma_gaussion_distribution;
  Eigen::VectorXd means = Eigen::VectorXd::Zero(dimension);

  random_samples::sample_student_t_distribution(dataset, sigma, means,
                                                degrees_of_freedom, 1000);
  // std::cout << dataset.transpose()  << std::endl;

  auto random_dataset = dataset.row(0);
  std::sort(std::begin(random_dataset), std::end(random_dataset),
            std::less<double>());
  std::cout << random_dataset.transpose() << std::endl;
}

TEST(RandomSamples, AlphaStableDistributionSamples) {
  unsigned int dimension = 5;
  Eigen::MatrixXd dataset;
  Eigen::MatrixXd Sigma = Eigen::MatrixXd::Identity(dimension, dimension);
  random_samples::sample_alpha_stable_distribution(dataset, Sigma, 1.9, 0, 0, 1,
                                                   1000);
  std::cout << dataset.transpose() << std::endl;
}

TEST(RandomSamples, AB_process) {
  std::random_device random{};
  std::mt19937 generator{random()};
  std::normal_distribution distributionA{0., 1.};
  std::normal_distribution distributionB{0., 1.};

  random_samples::noise_struct_tuple_type noise;
  noise.push_back(random_samples::noise_struct_type{
      [&]() { return distributionA(generator); }});
  noise.push_back(random_samples::noise_struct_type{
      [&]() { return distributionB(generator); }});

  Eigen::MatrixXd dataset;
  random_samples::AB_process(dataset, 0.1, 0.1, 0.1, 100, noise,
                             std::make_tuple(0., 0.));
  std::cout << dataset << std::endl;
}

TEST(RandomSamples, CB_process) {
  std::random_device random{};
  std::mt19937 generator{random()};
  std::normal_distribution distributionC{0., 1.};
  std::normal_distribution distributionB{0., 1.};
  random_samples::noise_struct_tuple_type noise;
  noise.push_back(random_samples::noise_struct_type{
      [&]() { return distributionC(generator); }});
  noise.push_back(random_samples::noise_struct_type{
      [&]() { return distributionB(generator); }});

  Eigen::MatrixXd dataset;
  random_samples::CB_process(dataset, 0.1, 0.1, 0.1, 100, noise,
                             std::make_tuple(0., 0.));
  std::cout << dataset << std::endl;
}

TEST(RandomSamples, ADB_process) {
  std::random_device random{};
  std::mt19937 generator{random()};
  std::normal_distribution distributionA{0., 1.};
  std::normal_distribution distributionD{0., 1.};
  std::normal_distribution distributionB{0., 1.};

  random_samples::noise_struct_tuple_type noise;
  noise.push_back(random_samples::noise_struct_type{
      [&]() { return distributionA(generator); }});
  noise.push_back(random_samples::noise_struct_type{
      [&]() { return distributionD(generator); }});
  noise.push_back(random_samples::noise_struct_type{
      [&]() { return distributionB(generator); }});

  Eigen::MatrixXd dataset;
  random_samples::ADB_process(dataset, 0.1, 0.1, 0.1, 0.1, 0.1, 100, noise,
                              std::make_tuple(0., 0., 0.));
  std::cout << dataset << std::endl;
}

TEST(RandomSamples, ACB_process) {
  std::random_device random{};
  std::mt19937 generator{random()};
  std::normal_distribution distributionA{0., 1.};
  std::normal_distribution distributionC{0., 1.};
  std::normal_distribution distributionB{0., 1.};

  random_samples::noise_struct_tuple_type noise;
  noise.push_back(random_samples::noise_struct_type{
      [&]() { return distributionA(generator); }});
  noise.push_back(random_samples::noise_struct_type{
      [&]() { return distributionC(generator); }});
  noise.push_back(random_samples::noise_struct_type{
      [&]() { return distributionB(generator); }});

  Eigen::MatrixXd dataset;
  random_samples::ACB_process(dataset, 0.1, 0.1, 0.1, 0.1, 0.1, 100, noise,
                              std::make_tuple(0., 0., 0.));
  std::cout << dataset << std::endl;
}

TEST(RandomSamples, A_process) {
  const unsigned int length = 1000;
  const double rho1 = 0.0;
  const double rho2 = 0.5;
  const double variance = 1.0;
  const double alpha = 2.0;
  const double skewness_parameter = 0.0;
  const unsigned long int cut_off_length = 1000;
  const unsigned int skip = 0U;
  int index = 0;
  for (auto &&item : random_samples::A_process<double>(
           length, rho1, rho2, variance, alpha, skewness_parameter,
           cut_off_length, skip)) {
    std::cout << std::format("{} {}\n", index, item);
    ++index;
  }
}

TEST(RandomSamples, B_process) {
  const unsigned int length = 1000;
  const double rho1 = 0.0;
  const double lambda = 0.0;
  const double rho2 = 0.5;
  const double rho3 = 0.0;
  const double variance = 1.0;
  const double alpha = 2.0;
  const double skewness_parameter = 0.0;
  const unsigned long int cut_off_length = 1000;
  const unsigned int skip = 0U;
  int index = 0;
  for (auto &&item : random_samples::generalized_B_process<double>(
           length, rho1, lambda, rho2, rho3, variance, alpha,
           skewness_parameter, cut_off_length, skip)) {
    std::cout << std::format("{} {}\n", index, item);
    ++index;
  }
}
