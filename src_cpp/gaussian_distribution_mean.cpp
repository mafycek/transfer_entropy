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

#include "random_samples.h"
#include "renyi_entropy.h"

int main ( int argc, char *argv[] )
{
    std::vector<std::string> methods{"LeonenkoProzanto",
                                     "LeonenkoWithGeneralizedMetric"};

    std::vector<std::string> random_noise_types
    {
        "uncorrelated_normal",    "correlated_normal",
        "full_correlated_normal", "uncorrelated_student-t",
        "correlated_student-t",   "full_correlated_student-t"};

    // Declare the supported options.
    boost::program_options::options_description desc ( "Allowed options" );
    std::string version ( "1.0.0" );
    desc.add_options() 
    ( "help", "Gaussian distribution analyzed by Renyi entropy" ) 
    ( "version", "Version of the program" ) 
    ( "dimension", boost::program_options::value<unsigned int>()->default_value ( 2 ), "Dimension" ) 
    ( "neighborhood", boost::program_options::value<unsigned int>()->default_value ( 21 ), "Maximal neighborhood" ) 
    ( "metric", boost::program_options::value<double>()->default_value ( 2 ), "Metric" ) 
    ( "method", boost::program_options::value<std::string>()->default_value ( methods[0] ), "Method - LeonenkoProzanto, LeonenkoWithGeneralizedMetric" ) 
    ( "random", boost::program_options::value<std::string>()->default_value ( random_noise_types[0] ), "Random noise - uncorrelated_normal, correlated_normal, full_correlated_normal, uncorrelated_student-t, correlated_student-t, full_correlated_student-t" ) 
    ( "sample", boost::program_options::value<unsigned int>()->default_value ( 1000 ), "Sample size" ) 
    ( "multithreading", "Multithreading" )
    ;

    boost::program_options::variables_map vm;
    boost::program_options::store (
        boost::program_options::parse_command_line ( argc, argv, desc ), vm );
    boost::program_options::notify ( vm );

    // show help
    if ( vm.count ( "help" ) )
    {
        std::cout << desc << std::endl;
        return 1;
    }

    bool multithreading = false;
    if ( vm.count ( "multithreading" ) )
    {
        multithreading = true;
    }

    // show version of the program
    if ( vm.count ( "version" ) )
    {
        std::cout << version << std::endl;
        return 1;
    }

    const double alpha_max = 5;
    const double mean_gaussion_distribution = 0;
    const double sigma_gaussion_distribution = 1;
    const unsigned int dimension = vm["dimension"].as<unsigned int>();
    const unsigned int neighborhood = vm["neighborhood"].as<unsigned int>();
    const double metric = vm["metric"].as<double>();
    const std::string method = vm["method"].as<std::string>();
    const std::string random_noise_type = vm["random"].as<std::string>();
    const unsigned int number_samples = vm["sample"].as<unsigned int>();

    std::cout << "N= " << number_samples << ", D= " << dimension
              << ", k_max= " << neighborhood << ", metric= " << metric
              << ", random_noise_type= " << random_noise_type
              << ", method=" << method << std::endl;

    const double delta_alpha = 0.005;
    Eigen::MatrixXd sigma;
    std::vector<std::vector<double>> dataset;
    Eigen::MatrixXd dataset2;
    std::vector<double> means ( dimension );
    std::function<double ( double, Eigen::MatrixXd & ) > renyi_entropy_function;

    auto start_dataset_generation_calculation =
        std::chrono::high_resolution_clock::now();
    if ( random_noise_type == random_noise_types[0] )
    {
        sigma = Eigen::MatrixXd::Identity ( dimension, dimension ) *
                sigma_gaussion_distribution * sigma_gaussion_distribution;
        random_samples::samples_normal_distribution_uncorrelated (
            dataset, 0, 1, number_samples, dimension );
        renyi_entropy_function = [] ( double q, Eigen::MatrixXd &sigma )
        {
            return random_samples::Renyi_entropy_normal_distribution ( q, sigma );
        };
    }
    else if ( random_noise_type == random_noise_types[1] )
    {
        double correlation = 0.4;
        double offdiagonal_element =
            sigma_gaussion_distribution * sigma_gaussion_distribution * correlation;
        Eigen::VectorXd means = Eigen::VectorXd::Zero ( dimension );
        sigma = Eigen::MatrixXd::Identity ( dimension, dimension ) *
                sigma_gaussion_distribution * sigma_gaussion_distribution;
        sigma.topRightCorner ( dimension - 1, dimension - 1 ) +=
            offdiagonal_element *
            Eigen::MatrixXd::Identity ( dimension - 1, dimension - 1 );
        sigma.bottomLeftCorner ( dimension - 1, dimension - 1 ) +=
            offdiagonal_element *
            Eigen::MatrixXd::Identity ( dimension - 1, dimension - 1 );
        sigma ( dimension - 1, 0 ) = offdiagonal_element;
        sigma ( 0, dimension - 1 ) = offdiagonal_element;

        random_samples::samples_normal_distribution_correlated (
            dataset2, means, sigma, number_samples );

        dataset = renyi_entropy::fromEigenMatrix<double, Eigen::MatrixXd> ( dataset2 );
        renyi_entropy_function = [] ( double q, Eigen::MatrixXd &sigma )
        {
            return random_samples::Renyi_entropy_normal_distribution ( q, sigma );
        };
    }
    else if ( random_noise_type == random_noise_types[2] )
    {
        double correlation = 0.4;
        double offdiagonal_element =
            sigma_gaussion_distribution * sigma_gaussion_distribution * correlation;
        Eigen::VectorXd means = Eigen::VectorXd::Zero ( dimension );
        sigma =
            Eigen::MatrixXd::Constant ( dimension, dimension, offdiagonal_element );
        sigma += Eigen::MatrixXd::Identity ( dimension, dimension ) *
                 ( 1 - correlation ) * sigma_gaussion_distribution *
                 sigma_gaussion_distribution;

        random_samples::samples_normal_distribution_correlated (
            dataset2, means, sigma, number_samples );
        dataset = renyi_entropy::fromEigenMatrix<double, Eigen::MatrixXd> ( dataset2 );
        renyi_entropy_function = [] ( double q, Eigen::MatrixXd &sigma )
        {
            return random_samples::Renyi_entropy_normal_distribution ( q, sigma );
        };
    }
    else if ( random_noise_type == random_noise_types[3] )
    {
        double degrees_of_freedom = 10;
        sigma = Eigen::MatrixXd::Identity ( dimension, dimension ) *
                sigma_gaussion_distribution * sigma_gaussion_distribution;
        Eigen::VectorXd means = Eigen::VectorXd::Zero ( dimension );
        random_samples::sample_student_t_distribution (
            dataset2, sigma, means, degrees_of_freedom, number_samples );
        dataset = renyi_entropy::fromEigenMatrix<double, Eigen::MatrixXd> ( dataset2 );
        renyi_entropy_function = [degrees_of_freedom] ( double q,
                                 Eigen::MatrixXd &sigma )
        {
            return random_samples::Renyi_student_t_distribution ( degrees_of_freedom, q,
                    sigma );
        };
    }
    else if ( random_noise_type == random_noise_types[4] )
    {
        double degrees_of_freedom = 2;
        double correlation = 0.4;
        double offdiagonal_element =
            sigma_gaussion_distribution * sigma_gaussion_distribution * correlation;
        sigma = Eigen::MatrixXd::Identity ( dimension, dimension ) *
                sigma_gaussion_distribution * sigma_gaussion_distribution;
        sigma.topRightCorner ( dimension - 1, dimension - 1 ) +=
            offdiagonal_element *
            Eigen::MatrixXd::Identity ( dimension - 1, dimension - 1 );
        sigma.bottomLeftCorner ( dimension - 1, dimension - 1 ) +=
            offdiagonal_element *
            Eigen::MatrixXd::Identity ( dimension - 1, dimension - 1 );
        sigma ( dimension - 1, 0 ) = offdiagonal_element;
        sigma ( 0, dimension - 1 ) = offdiagonal_element;
        Eigen::VectorXd means = Eigen::VectorXd::Zero ( dimension );
        random_samples::sample_student_t_distribution (
            dataset2, sigma, means, degrees_of_freedom, number_samples );
        dataset = renyi_entropy::fromEigenMatrix<double, Eigen::MatrixXd> ( dataset2 );
        renyi_entropy_function = [degrees_of_freedom] ( double q,
                                 Eigen::MatrixXd &sigma )
        {
            return random_samples::Renyi_student_t_distribution ( degrees_of_freedom, q,
                    sigma );
        };
    }
    else if ( random_noise_type == random_noise_types[5] )
    {
        double degrees_of_freedom = 2;
        double correlation = 0.4;
        double offdiagonal_element =
            sigma_gaussion_distribution * sigma_gaussion_distribution * correlation;
        sigma =
            Eigen::MatrixXd::Constant ( dimension, dimension, offdiagonal_element );
        sigma += Eigen::MatrixXd::Identity ( dimension, dimension ) *
                 ( 1 - correlation ) * sigma_gaussion_distribution *
                 sigma_gaussion_distribution;
        Eigen::VectorXd means = Eigen::VectorXd::Zero ( dimension );
        random_samples::sample_student_t_distribution (
            dataset2, sigma, means, degrees_of_freedom, number_samples );
        dataset = renyi_entropy::fromEigenMatrix<double, Eigen::MatrixXd> ( dataset2 );
        renyi_entropy_function = [degrees_of_freedom] ( double q,
                                 Eigen::MatrixXd &sigma )
        {
            return random_samples::Renyi_student_t_distribution ( degrees_of_freedom, q,
                    sigma );
        };
    }
    auto end_dataset_generation_calculation =
        std::chrono::high_resolution_clock::now();
    auto elapsed_edata_generation =
        end_dataset_generation_calculation - start_dataset_generation_calculation;
    auto milliseconds_elapsed_edata_generation =
        std::chrono::duration_cast<std::chrono::milliseconds> (
            elapsed_edata_generation )
        .count();
    BOOST_LOG_TRIVIAL ( trace ) << std::format (
                                    "Dataset generated in {:.5f} seconds",
                                    milliseconds_elapsed_edata_generation / 1000.0 );

    std::vector<double> alphas;
    for ( double alpha = delta_alpha; alpha < 1; alpha += delta_alpha )
    {
        alphas.push_back ( alpha );
    }
    alphas.push_back ( 0.999 );
    alphas.push_back ( 1 );
    alphas.push_back ( 1.001 );
    for ( double alpha = 1 + delta_alpha; alpha < alpha_max;
            alpha += delta_alpha )
    {
        alphas.push_back ( alpha );
    }

    renyi_entropy::renyi_entropy<double> calculator;
    auto neighbors = std::ranges::iota_view{1U, neighborhood};
    std::vector<unsigned int> indices (
        neighbors.begin(),
        neighbors.end() ); //{1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13,
    //14, 15, 16, 17, 18, 19, 20});

    calculator.SetAlpha ( alphas );
    calculator.SetIndices ( indices );
    calculator.SetExp ( [&] (double x) { return exp(x);} );
    calculator.SetLog ( [&] (double x) { return log(x);} );
    calculator.SetPower ( [&] (double x, double y) { return pow(x, y);} );
    calculator.SetMultithreading ( multithreading );

    std::map<std::tuple<unsigned int, double>, double> result;
    if ( method == methods[0] )
    {
        result = calculator.renyi_entropy_LeonenkoProzanto ( dataset, 2 );
    }
    else if ( method == methods[1] )
    {
        result = calculator.renyi_entropy_metric ( dataset, metric );
    }

    renyi_entropy::renyi_entropy<double>::SaveRenyiEntropy ( result, "myfile.dat" );

    for ( const auto &item_alpha : calculator.GetAlphas() )
    {
        std::cout << item_alpha << " ";
        for ( const auto &item_index : calculator.GetIndices() )
        {
            auto result_entropy = result[std::make_tuple ( item_index, item_alpha )];
            std::cout << result_entropy << " ";
        }
        auto theoretical_renyi_entropy =
            renyi_entropy_function ( static_cast<double> ( item_alpha ), sigma );
        std::cout << theoretical_renyi_entropy << " ";
        for ( const auto &item_index : calculator.GetIndices() )
        {
            auto result_relative_entropy =
                theoretical_renyi_entropy /
                result[std::make_tuple ( item_index, item_alpha )];
            std::cout << result_relative_entropy << " ";
        }
        std::cout << std::endl;
    }
}
