#pragma once

#include <iostream>
#include <random>
#include <ranges>
#include <vector>
#include <functional>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

namespace random_samples
{

double full_correlation_matrix ( int dimension, double q );

double tridiagonal_matrix_determinant ( int dimension, double q );

double Renyi_entropy_normal_distribution ( const double q,
        const Eigen::MatrixXd &Sigma );

double Renyi_student_t_distribution ( double degrees_of_freedom, const double q,
                                      const Eigen::MatrixXd &sigma );

void samples_normal_distribution_uncorrelated (
    std::vector<std::vector<double>> &dataset,
    double mean_gaussion_distribution, double sigma_gaussion_distribution,
    unsigned int number_samples, unsigned int dimension );

void samples_normal_distribution_uncorrelated (
    Eigen::MatrixXd &dataset, double mean_gaussion_distribution,
    double sigma_gaussion_distribution, unsigned int number_samples,
    unsigned int dimension );

void samples_normal_distribution_correlated (
    Eigen::MatrixXd &dataset, const Eigen::VectorXd &mean_gaussion_distribution,
    const Eigen::MatrixXd &Sigma, unsigned int number_samples );

void sample_student_t_distribution ( Eigen::MatrixXd &dataset,
                                     const Eigen::MatrixXd &Sigma,
                                     const Eigen::VectorXd &mean,
                                     const double degrees_of_freedom,
                                     unsigned int number_samples );

void sample_alpha_stable_distribution ( Eigen::MatrixXd &dataset,
                                        const Eigen::MatrixXd &Sigma,
                                        const double alpha, const double beta,
                                        const double mu, const double param,
                                        unsigned int number_samples );

void AB_process(Eigen::MatrixXd &dataset, double alpha, double beta, double eta, unsigned int number_samples, std::tuple<std::function<double ()>, std::function<double ()>> random_generator, std::tuple<double, double> initial_confition = std::tuple<double, double>({0.,0.}));

void CB_process(Eigen::MatrixXd &dataset, double alpha, double beta, double eta, unsigned int number_samples, std::tuple<std::function<double ()>, std::function<double ()>> random_generator, std::tuple<double, double> initial_confition = std::tuple<double, double>({0.,0.}), double lambda=3);

void ADB_process(Eigen::MatrixXd &dataset, double alpha, double beta, double gamma, double eta1, double eta2, unsigned int number_samples, std::tuple<std::function<double ()>, std::function<double ()>, std::function<double ()>> random_generator, std::tuple<double, double, double> initial_confition = std::tuple<double, double, double>({0.,0.,0,}));

void ACB_process(Eigen::MatrixXd &dataset, double alpha, double beta, double gamma, double eta1, double eta2, unsigned int number_samples, std::tuple<std::function<double ()>, std::function<double ()>, std::function<double ()>> random_generator, std::tuple<double, double, double> initial_confition = std::tuple<double, double, double>({0.,0.,0,}), double lambda=3);

} // namespace random_samples
