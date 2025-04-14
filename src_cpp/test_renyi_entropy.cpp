//
// Created by hynek on 24.3.25.
//

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <gtest/gtest.h>
#include <iostream>

#include "renyi_entropy.h"

double distance(const std::vector<double> &a, const std::vector<double> &b)
{
    if (a.size() == b.size()) {
        double distance = 0;
        for (int i = 0; i < a.size(); ++i) {
            distance += pow(a[i] - b[i], 2);
        }
        return sqrt(distance);
    }
    return NAN;
}

TEST(RenyiEntropy, SimpleTest)
{
    std::vector<std::vector<double>> dataset{{1, 1, 1}, {2, 1, 1}, {1, 1, 2},
        {2, 1, 2}, {3, 1, 1}, {1, 1, 3},
        {3, 3, 1}, {3, 3, 3}};
    renyi_entropy::renyi_entropy<double> calculator;
    auto tree = calculator.create_KDtree(dataset);

    int nearest = 5;
    std::vector<std::vector<double>> distances;
    for (const auto &point : dataset) {
        auto points = tree.nearest_points(point, nearest);
        std::vector<double> distance_vector;
        for (int i = 0; i < nearest; ++i) {
            distance_vector.push_back(distance(point, points[i]));
        }
        distances.push_back(distance_vector);
    }
    const int use_index = 2;
    auto maximum_distance_element =
        std::max_element(distances.begin(), distances.end(),
    [use_index](const auto & first, const auto & second) {
        if (first[use_index] > second[use_index]) {
            return true;
        } else {
            return false;
        }
    });
    auto maximum_distance = (*maximum_distance_element)[use_index];

    std::cout << maximum_distance << std::endl;

    std::vector<double> alpha{0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2};
    std::vector<unsigned int> indices{1, 2, 3, 4, 5};
    calculator.SetAlpha(alpha);
    calculator.SetIndices(indices);
    calculator.SetExp(exp);
    calculator.SetLog(log);
    calculator.SetPower(pow);

    auto result = calculator.renyi_entropy_LeonenkoProzanto(dataset, 2);
    std::stringstream ss;
    boost::filesystem::path myFile =
        boost::filesystem::current_path() / "myfile.dat";
    boost::filesystem::ofstream ofs(myFile);
    boost::archive::binary_oarchive oarch(ofs);
    oarch << result;
    std::cout << ss.str();

    EXPECT_EQ(0, 0);
}
