
#include <gtest/gtest.h>
#include <iostream>
#include <string>

#include <eigen3/Eigen/Core>

#include "utils.h"

TEST ( Utils, ArrayExtraction )
{
    std::string test_string{"1.2 2 3 4 5 1e220"};
    std::vector<double> values;

    convert_array ( test_string, values, std::stod );
    for ( auto item : values )
    {
        std::cout << item << std::endl;
    }
}

TEST ( Utils, EigenMatrixToVector )
{
    unsigned int columns{5}, rows{5};
    Eigen::MatrixXd dataset ( columns, rows );
    for ( unsigned int i = 0; i < columns; ++i )
    {
        for ( unsigned int j = 0; j < rows; ++j )
        {
            dataset ( i, j ) = i + j * 2;
        }
    }

    auto item = dataset.col ( 0 );
    std::cout << item << std::endl;
}

TEST ( Utils, ChopArrayExtraction )
{
    std::string test_string{"1.2 2, 3 4, 5 1e220"};
    std::vector<std::vector<double>> values;

    chop_string_arrays ( test_string, ',', values, std::stod );
    for ( auto &item : values )
    {
        for ( auto &item2 : item )
        {
            std::cout << item2 << " ";
        }
        std::cout << std::endl;
    }
}
