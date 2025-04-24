
#include "random_samples.h"

#include <list>
#include <random>
#include <exception>
#include <cmath>

#include <boost/math/special_functions/gamma.hpp>

#include <eigen3/Eigen/Eigenvalues>

namespace random_samples {

double full_correlation_matrix(int dimension, double q)
{
    return pow(1 - q, dimension - 1) * (1 + (dimension - 1) * q);
}


double  tridiagonal_matrix_determinant(int dimension, double q)
{
  if (dimension == 1)
  {
        return 1;
  } else if (dimension == 2)
  {
        return 1 - pow (q, 2);
  }
  else
  {
        std::vector<double> sample {1, 1 - pow (q, 2)};
        for (int n{2}; n < dimension; ++n)
        {
            auto result = sample[n - 1] - pow(q, 2) * sample[n - 2];
            sample.push_back(result);
        }
        return sample[dimension - 1];
  }
}

  
double Renyi_entropy_normal_distribution(const double q, const Eigen::MatrixXd &Sigma) {
  if (Sigma.rows() == Sigma.cols()) {
    auto m = Sigma.cols();
    if (q != 1) {
      return log(2 * std::numbers::pi_v<double>) * m / 2 +
             log(Sigma.eval().determinant()) / 2 - m * log(q) / (1 - q) / 2;
    } else {
      return log(2 * std::numbers::pi_v<double> * std::exp(1)) * m / 2 +
             log(Sigma.eval().determinant()) / 2;
    }
  } else {
    std::cerr << "Nonsymmetric matrix" << std::endl;
    return 0;
  }
}

void samples_normal_distribution_uncorrelated(std::vector<std::vector<double>> &dataset, double mean_gaussion_distribution, double sigma_gaussion_distribution, unsigned int number_samples, unsigned int dimension)
{
  std::random_device random{};
  std::mt19937 generator{random()};
  std::normal_distribution distribution{mean_gaussion_distribution, sigma_gaussion_distribution};

  auto random_double = [&distribution, &generator] {
    return distribution(generator);
  };

  for (auto n{number_samples}; n; --n) 
  {
    dataset.push_back(std::vector<double>());
    for (int i : std::ranges::iota_view<unsigned int, unsigned int>(0, dimension)) 
    {
      dataset.back().push_back(random_double());
    }
  }
}

void samples_normal_distribution_uncorrelated(Eigen::MatrixXd &dataset, double mean_gaussion_distribution, double sigma_gaussion_distribution, unsigned int number_samples, unsigned int dimension)
{
  std::random_device random{};
  std::mt19937 generator{random()};
  std::normal_distribution distribution{mean_gaussion_distribution, sigma_gaussion_distribution};

  auto random_double = [&distribution, &generator] {
    return distribution(generator);
  };

  dataset.conservativeResize(dimension, number_samples);
  for (auto n{0}; n<number_samples; ++n) 
  {
    Eigen::VectorXd random_vector(dimension);
    for (int i : std::ranges::iota_view<unsigned int, unsigned int>(0, dimension)) 
    {
      random_vector[i] = random_double();
    }
    dataset.col(n) = random_vector;
  }
}

void samples_normal_distribution_correlated(Eigen::MatrixXd &dataset, const Eigen::VectorXd & mean_gaussion_distribution, const Eigen::MatrixXd &Sigma, unsigned int number_samples)
{
  if (Sigma.cols() == Sigma.rows())
  {
    auto dimension = Sigma.cols();

    Eigen::MatrixXd normal_dataset;

    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> eigenValueSolver(Sigma, Eigen::MatrixXd::Identity(Sigma.cols(), Sigma.cols()));
    const auto eigenvalueMatrix = eigenValueSolver.eigenvalues();
    const auto eigenvectorMatrix = eigenValueSolver.eigenvectors();
    const auto diagonal_matrix = eigenvalueMatrix.array().sqrt().matrix().asDiagonal();
    const auto transform = eigenvectorMatrix * diagonal_matrix;
    //std::cout << eigenvalueMatrix << std::endl;

    samples_normal_distribution_uncorrelated(normal_dataset, 0., 1., number_samples, dimension);
    //std::cout << transform << " " << eigenvalueMatrix << " " << eigenvectorMatrix << std::endl;
    dataset = (transform * normal_dataset).colwise() + mean_gaussion_distribution;
    //dataset = dataset.colwise() + mean_gaussion_distribution;
  }
  else
  {
    //throw Error("sigma parameter has wrong type");
  }
}

void sample_student_t_distribution(Eigen::MatrixXd &dataset, const Eigen::MatrixXd &Sigma, const Eigen::VectorXd & mean, const double degrees_of_freedom, unsigned int number_samples)
{
    std::random_device random{};
    std::mt19937 generator{random()};
    std::chi_squared_distribution<double> chi2_generator{degrees_of_freedom};
    
    auto random_double = [&chi2_generator, &generator] {
      return chi2_generator(generator);
    };
    
    Eigen::MatrixXd dataset_normal;
    samples_normal_distribution_correlated(dataset_normal, mean, Sigma, number_samples);

    auto dimension = Sigma.cols();
    dataset.conservativeResize(dimension, number_samples);
    for (unsigned int i{0}; i < number_samples; ++i)
    {
        auto random_chi2_number = random_double();
        auto random_vector = mean + dataset_normal.col(i) * ( sqrt( degrees_of_freedom / random_chi2_number ) );
        dataset.col(i) = random_vector;
    }
    //std::cout << dataset.transpose() << std::endl;
}

double Renyi_student_t_distribution(double degrees_of_freedom, const double q, const Eigen::MatrixXd & sigma )
{
    auto dimension = sigma.cols();
    if (q * (degrees_of_freedom + dimension) / 2 - dimension / 2 <= 0 )
    {
        return std::nan("");
    }

    auto arg = pow(sigma.eval().determinant(), (1.0 - q) / 2);

    auto arg2 = pow(degrees_of_freedom * std::numbers::pi_v<double>, dimension / 2 * (1.0 - q));
    auto arg31 = pow(boost::math::tgamma((degrees_of_freedom + dimension) / 2), q);
    auto arg32 = boost::math::tgamma(q * (degrees_of_freedom + dimension) / 2 - dimension / 2);
    auto arg33 = boost::math::tgamma(q * (degrees_of_freedom + dimension) / 2);
    auto arg34 = pow(boost::math::tgamma(degrees_of_freedom / 2), q);
    auto arg3 = (arg31 * arg32) / (arg33 * arg34);

    return 1 / (1 - q) * log2(arg3 * arg * arg2);
}

}
