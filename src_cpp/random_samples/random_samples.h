#pragma once

#include <algorithm>
#include <cmath>
#include <coroutine>
#include <ctime>
#include <functional>
#include <generator>
#include <iostream>
#include <list>
#include <map>
#include <random>
#include <ranges>
#include <tuple>
#include <vector>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/mpfr.hpp>

#include <boost/math/special_functions.hpp>
#include <boost/math/special_functions/gamma.hpp>

namespace random_samples {

using noise_struct_type =
    std::variant<std::function<double()>, std::generator<double>>;
using noise_struct_tuple_type = std::vector<random_samples::noise_struct_type>;

double full_correlation_matrix(int dimension, double q);

double tridiagonal_matrix_determinant(int dimension, double q);

double Renyi_entropy_normal_distribution(const double q,
                                         const Eigen::MatrixXd &Sigma);

double Renyi_student_t_distribution(double degrees_of_freedom, const double q,
                                    const Eigen::MatrixXd &sigma);

double normal_distribution_sample(const double mean_gaussion_distribution,
                                  const double sigma_gaussion_distribution);

std::vector<double>
normal_distribution_uncorrelated(double mean_gaussion_distribution,
                                 double sigma_gaussion_distribution,
                                 unsigned int dimension);

void samples_normal_distribution_uncorrelated(
    std::vector<std::vector<double>> &dataset,
    double mean_gaussion_distribution, double sigma_gaussion_distribution,
    unsigned int number_samples, unsigned int dimension);

void samples_normal_distribution_uncorrelated(
    Eigen::MatrixXd &dataset, double mean_gaussion_distribution,
    double sigma_gaussion_distribution, unsigned int number_samples,
    unsigned int dimension);

void samples_normal_distribution_correlated(
    Eigen::MatrixXd &dataset, const Eigen::VectorXd &mean_gaussion_distribution,
    const Eigen::MatrixXd &Sigma, unsigned int number_samples);

void sample_student_t_distribution(Eigen::MatrixXd &dataset,
                                   const Eigen::MatrixXd &Sigma,
                                   const Eigen::VectorXd &mean,
                                   const double degrees_of_freedom,
                                   unsigned int number_samples);

void sample_alpha_stable_distribution(Eigen::MatrixXd &dataset,
                                      const Eigen::MatrixXd &Sigma,
                                      const double alpha, const double beta,
                                      const double mu, const double param,
                                      unsigned int number_samples);

void AB_process(Eigen::MatrixXd &dataset, const double alpha, const double beta,
                const double eta, const unsigned int number_samples,
                noise_struct_tuple_type &random_generator,
                const std::tuple<double, double> initial_confition =
                    std::tuple<double, double>({0., 0.}));

void CB_process(Eigen::MatrixXd &dataset, const double alpha, const double beta,
                const double eta, const unsigned int number_samples,
                noise_struct_tuple_type &random_generator,
                const std::tuple<double, double> initial_confition =
                    std::tuple<double, double>({0., 0.}),
                const double lambda = 3);

void ADB_process(Eigen::MatrixXd &dataset, const double alpha,
                 const double beta, const double gamma, const double eta1,
                 const double eta2, const unsigned int number_samples,
                 noise_struct_tuple_type &random_generator,
                 const std::tuple<double, double, double> initial_confition =
                     std::tuple<double, double, double>({
                         0.,
                         0.,
                         0,
                     }));

void ACB_process(Eigen::MatrixXd &dataset, const double alpha,
                 const double beta, const double gamma, const double eta1,
                 const double eta2, const unsigned int number_samples,
                 noise_struct_tuple_type &random_generator,
                 const std::tuple<double, double, double> initial_confition =
                     std::tuple<double, double, double>({
                         0.,
                         0.,
                         0,
                     }),
                 const double lambda = 3);

template <class TYPE>
void make_cached_coefficients(std::vector<TYPE> &precalculated,
                              const double rho, const unsigned int number) {
  precalculated.resize(number);
  for (int i = 1; i <= number; i++) {
    precalculated[i - 1] =
        static_cast<TYPE>(rho) *
        (boost::math::expm1(boost::math::lgamma(static_cast<TYPE>(i - rho)) -
                            boost::math::lgamma(static_cast<TYPE>(1 - rho)) -
                            boost::math::lgamma(static_cast<TYPE>(1 + i))) +
         static_cast<TYPE>(1));
  }
}

template <class TYPE>
TYPE mean_absolute_value(const std::list<TYPE> &serie,
                         std::function<TYPE(TYPE)> absolute_value) {
  TYPE inter_result = static_cast<TYPE>(0.0);
  for (auto &item : serie) {
    inter_result += absolute_value(item);
  }

  return (inter_result / static_cast<TYPE>(serie.size()));
}

template <class TYPE>
std::generator<TYPE> A_process(
    const unsigned int length = 1000, const double rho1 = 0.0,
    const double rho2 = 0.5, const double variance = 1.0,
    const double alpha = 2.0, const double skewness_parameter = 0.0,
    const unsigned int cut_off_length = 1000, const unsigned int skip = 0U,
    const unsigned long int seed = 0,
    std::function<TYPE(TYPE)> absolute_value = [](double value) -> double {
      return abs(value);
    }) {
  std::vector<TYPE> precalculated_coef_autocorrelation;
  std::vector<TYPE> precalculated_coef_distribution;

  make_cached_coefficients(precalculated_coef_autocorrelation, rho1,
                           cut_off_length);
  make_cached_coefficients(precalculated_coef_distribution, rho2,
                           cut_off_length);

  std::list<TYPE> differences;
  std::list<TYPE> noise;
  std::list<TYPE> noise_volatility;

  std::random_device random_device;
  const gsl_rng_type *rng_type_instance;
  gsl_rng *rnd_generator;

  gsl_rng_env_setup();

  rng_type_instance = gsl_rng_default;
  rnd_generator = gsl_rng_alloc(rng_type_instance);
  gsl_rng_set(rnd_generator, seed);

  // loop over generated data
  for (unsigned int length_of_data :
       std::ranges::iota_view{1U, length + skip + 1}) {
    TYPE noise_gen = static_cast<TYPE>(
        gsl_ran_levy_skew(rnd_generator, variance, alpha, skewness_parameter));
    noise.push_back(noise_gen);
    TYPE sigma;
    TYPE sigma_autocorrelation = TYPE(0);

    if (length_of_data == 1) {
      sigma = static_cast<TYPE>(1);
    } else {
      TYPE inter_sigma = TYPE(0);
      typename std::list<TYPE>::reverse_iterator iter = differences.rbegin();
      for (unsigned int i = 1;
           i <= std::min(cut_off_length,
                         static_cast<unsigned int>(differences.size()));
           i++) {
        TYPE value = (*iter);
        inter_sigma +=
            precalculated_coef_distribution[i - 1] * absolute_value(value);
        sigma_autocorrelation +=
            precalculated_coef_autocorrelation[i - 1] * value;
        iter++;
      }
      sigma = inter_sigma / mean_absolute_value(differences, absolute_value);

      if (rho2 <= 0.0) {
        sigma = TYPE(1.0);
      }
    }

    TYPE next_result = sigma_autocorrelation + sigma * noise_gen;
    differences.push_back(next_result);
    if (length_of_data >= skip) {
      co_yield static_cast<TYPE>(next_result);
    }
  }
}

template <class TYPE>
TYPE Scalar_Product(std::list<TYPE> noise, std::list<TYPE> noise_volatility,
                    unsigned int max_count = 1000) {
  TYPE result = static_cast<TYPE>(0.0);
  typename std::list<TYPE>::reverse_iterator iter_noise = noise.rbegin();
  typename std::list<TYPE>::reverse_iterator iter_noise_volatility =
      noise_volatility.rbegin();
  for (int count = 0; ((iter_noise != noise.rend()) &&
                       (iter_noise_volatility != noise_volatility.rend()) &&
                       (count < max_count));
       ++count) {
    result += (*iter_noise) * (*iter_noise_volatility);
    ++iter_noise;
    ++iter_noise_volatility;
  }

  return (result);
}

// generalization of the stochastic B process for non-symmetric distributions
template <class TYPE>
std::generator<TYPE> generalized_B_process(
    const unsigned int length = 1000, const double rho1 = 0.0,
    const double lambda = 0.0, const double rho2 = 0.5, const double rho3 = 0.0,
    const double variance = 1.0, const double alpha = 2.0,
    const double skewness_parameter = 0.0,
    const unsigned long int cut_off_length = 1000U,
    const unsigned int skip = 0U, const unsigned long int seed = 0U,
    std::function<TYPE(TYPE)> absolute_value = [](double value) -> double {
      return abs(value);
    }) {
  std::vector<TYPE> precalculated_coef_autocorrelation(
      static_cast<size_t>(cut_off_length), 0.);
  std::vector<TYPE> precalculated_coef_distribution{
      static_cast<size_t>(cut_off_length), 0.};
  std::vector<TYPE> precalculated_coef_noise_correlation{
      static_cast<size_t>(cut_off_length), 0.};

  make_cached_coefficients(precalculated_coef_autocorrelation, rho1,
                           cut_off_length);
  make_cached_coefficients(precalculated_coef_distribution, rho2,
                           cut_off_length);
  make_cached_coefficients(precalculated_coef_noise_correlation, rho3,
                           cut_off_length);

  std::list<TYPE> differences{0};
  std::list<TYPE> noise;
  std::list<TYPE> noise_volatility;

  std::random_device random_device;
  const gsl_rng_type *rng_type_instance;
  gsl_rng *rnd_generator;

  gsl_rng_env_setup();

  rng_type_instance = gsl_rng_default;
  rnd_generator = gsl_rng_alloc(rng_type_instance);
  gsl_rng_set(rnd_generator, seed);

  for (unsigned int length_of_data :
       std::ranges::iota_view{1U, length + skip + 1}) {
    const TYPE noise_gen = static_cast<TYPE>(
        gsl_ran_levy_skew(rnd_generator, variance, alpha, skewness_parameter));

    noise.push_back(noise_gen);
    TYPE sigma;
    TYPE sigma_autocorrelation = static_cast<TYPE>(0);

    if (length_of_data == 1) {
      sigma = static_cast<TYPE>(1);
    } else {
      TYPE inter_sigma = static_cast<TYPE>(0);
      typename std::list<TYPE>::reverse_iterator iter = differences.rbegin();
      for (unsigned int i = 1;
           i <= std::min(cut_off_length,
                         static_cast<unsigned long int>(differences.size()));
           i++) {
        TYPE value = (*iter);
        typename std::list<TYPE>::reverse_iterator iter2 = iter;
        TYPE previous_value = (*(iter2++));
        inter_sigma +=
            precalculated_coef_distribution[i - 1] * absolute_value(value);
        sigma_autocorrelation +=
            precalculated_coef_autocorrelation[i - 1] *
            (value - lambda * absolute_value(previous_value));
        ++iter;
      }
      const auto abs_value = mean_absolute_value(differences, absolute_value);
      sigma = inter_sigma / abs_value;

      if (rho2 <= 0.0) {
        sigma = static_cast<TYPE>(1.0);
      }
    }

    noise_volatility.push_back(sigma);

    TYPE next_result = lambda * abs(differences.back()) +
                       sigma_autocorrelation +
                       Scalar_Product(noise_volatility, noise, cut_off_length);
    differences.push_back(next_result);
    if (length_of_data >= skip) {
      co_yield static_cast<TYPE>(next_result);
    }
  }
  gsl_rng_free(rnd_generator);
}

} // namespace random_samples
// kate: indent-mode cstyle; indent-width 2; replace-tabs on; ;
