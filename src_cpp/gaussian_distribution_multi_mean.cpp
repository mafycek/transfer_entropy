//
// Created by hynek on 11.4.25.
//

#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <string>

#include <gtest/gtest.h>

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/float128.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/number.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/tommath.hpp>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/math/special_functions/expm1.hpp>

#include "renyi_entropy.h"

using namespace std;

namespace boost { namespace serialization { // insert this code to the appropriate namespaces


/**
 Save a mpfr_float type to a boost archive.
 */
template <typename Archive>
void save(Archive& ar, boost::multiprecision::cpp_bin_float<100 > const& r, unsigned /*version*/)
{
    std::string tmp = r.str(0, std::ios::fixed);// 0 indicates use full precision
    ar & tmp;
}

/**
 Load a mpfr_float type from a boost archive.
 */
template <typename Archive>
void load(Archive& ar, boost::multiprecision::cpp_bin_float<100 >& r, unsigned /*version*/)
{
    std::string tmp;
    ar & tmp;
    r = tmp.c_str();
}

} } // re: namespaces


int main(int argc, char *argv[])
{
        double mean_gaussion_distribution = 0;
        double sigma_gaussion_distribution = 10;
        typedef boost::multiprecision::number <
        boost::multiprecision::cpp_bin_float<20 >>
                        multiprecision_type;

        auto exp_m = [&] ( multiprecision_type argument ) {
                return boost::math::expm1 ( argument ) + 1;
        };
        auto log_m = [&] ( multiprecision_type argument ) {
                return boost::math::log1p ( argument ) - 1;
        };
        auto power_m = [&] ( multiprecision_type argument,
        multiprecision_type exponent ) {
                return boost::multiprecision::pow ( argument, exponent );
        };
        auto number_samples = 100;
        std::random_device random{};
        std::mt19937 generator{random() };
        std::normal_distribution distribution{mean_gaussion_distribution, sigma_gaussion_distribution};

        auto random_double = [&distribution, &generator] {
                return distribution ( generator );
        };

        std::vector<std::vector<multiprecision_type>> dataset;
        for ( auto n{number_samples}; n; --n )
                dataset.push_back ( std::vector<multiprecision_type> {random_double() } );

        std::vector<multiprecision_type> alphas;
        for ( multiprecision_type alpha = 0.05; alpha <= 5; alpha += multiprecision_type(0.01) ) {
                alphas.push_back ( alpha );
        }

        renyi_entropy::renyi_entropy<multiprecision_type> calculator;
        std::vector<unsigned int> indices{1, 2, 3, 4, 5};

        calculator.SetAlpha ( alphas );
        calculator.SetIndices ( indices );
        calculator.SetExp ( exp_m );
        calculator.SetLog ( log_m );
        calculator.SetPower ( power_m );

        auto result = calculator.renyi_entropy_LeonenkoProzanto ( dataset, 2 );

        std::stringstream ss;
        boost::filesystem::path myFile =
                boost::filesystem::current_path() / "myfile.dat";
        boost::filesystem::ofstream ofs ( myFile );
        boost::archive::binary_oarchive oarch ( ofs );
        //oarch << result;
        std::cout << ss.str();

        std::cout << std::setprecision ( std::numeric_limits<multiprecision_type>::max_digits10 ) << std::endl;
        Eigen::MatrixXd sigma{{sigma_gaussion_distribution}};
        for( const auto & item_alpha: calculator.GetAlphas())
        {
                std::cout << item_alpha << " ";
                for(const auto & item_index: calculator.GetIndices())
                {
                        auto result_entropy = result[std::make_tuple(item_index, item_alpha)];
                        std::cout << result_entropy << " ";
                }
                std::cout << renyi_entropy::Renyi_entropy_normal_distribution(1, static_cast<double>(item_alpha), sigma)  << std::endl;
        }
}
