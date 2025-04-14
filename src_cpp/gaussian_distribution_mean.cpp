//
// Created by hynek on 9.4.25.
//

#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <string>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/multiprecision/number.hpp>
#include <gtest/gtest.h>

#include "renyi_entropy.h"

int main(int argc, char *argv[])
{
        double mean_gaussion_distribution = 0;
        double sigma_gaussion_distribution = 1;
        auto number_samples = 100000;
        std::random_device random{};
        std::mt19937 generator{random() };
        std::normal_distribution distribution{mean_gaussion_distribution, sigma_gaussion_distribution};

        auto random_double = [&distribution, &generator] {
                return distribution ( generator );
        };

        std::vector<std::vector<double>> dataset;
        for ( auto n{number_samples}; n; --n )
                dataset.push_back ( std::vector<double> {random_double(), random_double() } );

        std::vector<double> alphas;
        double delta_alpha = 0.05;
        for ( double alpha = 0.05; alpha < 1; alpha += delta_alpha ) {
                alphas.push_back ( alpha );
        }
        alphas.push_back ( 1 );
        for ( double alpha = 1+ delta_alpha; alpha < 5; alpha += delta_alpha ) {
                alphas.push_back ( alpha );
        }        

        renyi_entropy::renyi_entropy<double> calculator;
        std::vector<unsigned int> indices{1, 2, 3, 4, 5, 6 ,7, 8, 9, 10};

        calculator.SetAlpha ( alphas );
        calculator.SetIndices ( indices );
        calculator.SetExp ( exp );
        calculator.SetLog ( log );
        calculator.SetPower ( pow );

        auto result = calculator.renyi_entropy_LeonenkoProzanto ( dataset, 2 );

        std::stringstream ss;
        boost::filesystem::path myFile = boost::filesystem::current_path() / "myfile.dat";
        boost::filesystem::ofstream ofs ( myFile );
        boost::archive::binary_oarchive oarch ( ofs );
        oarch << result;
        std::cout << ss.str();

        Eigen::MatrixXd sigma{{sigma_gaussion_distribution*sigma_gaussion_distribution, 0}, {0, sigma_gaussion_distribution*sigma_gaussion_distribution}};
        for( const auto & item_alpha: calculator.GetAlphas())
        {
                std::cout << item_alpha << " ";
                for(const auto & item_index: calculator.GetIndices())
                {
                        auto result_entropy = result[std::make_tuple(item_index, item_alpha)];
                        std::cout << result_entropy << " ";
                }
                std::cout << renyi_entropy::Renyi_entropy_normal_distribution(1, static_cast<double>( item_alpha ), sigma)  << std::endl;
        }
}
