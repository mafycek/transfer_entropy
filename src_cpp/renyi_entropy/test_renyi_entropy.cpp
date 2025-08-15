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

double distance ( const std::vector<double> &a, const std::vector<double> &b )
{
    if ( a.size() == b.size() )
    {
        double distance = 0;
        for ( int i = 0; i < a.size(); ++i )
        {
            distance += pow ( a[i] - b[i], 2 );
        }
        return sqrt ( distance );
    }
    return NAN;
}

TEST ( RenyiEntropy, SimpleTest )
{
    std::vector<std::vector<double>> dataset{{1, 1, 1}, {2, 1, 1}, {1, 1, 2},
        {2, 1, 2}, {3, 1, 1}, {1, 1, 3},
        {3, 3, 1}, {3, 3, 3}};
    renyi_entropy::renyi_entropy<double> calculator;
    auto tree = calculator.create_KDtree ( dataset );

    int nearest = 5;
    std::vector<std::vector<double>> distances;
    for ( const auto &point : dataset )
    {
        auto points = tree.nearest_points ( point, nearest );
        std::vector<double> distance_vector;
        for ( int i = 0; i < nearest; ++i )
        {
            distance_vector.push_back ( distance ( point, points[i] ) );
        }
        distances.push_back ( distance_vector );
    }
    const int use_index = 2;
    auto maximum_distance_element =
        std::max_element ( distances.begin(), distances.end(),
                           [use_index] ( const auto &first, const auto &second )
    {
        if ( first[use_index] > second[use_index] )
        {
            return true;
        }
        else
        {
            return false;
        }
    } );
    auto maximum_distance = ( *maximum_distance_element ) [use_index];

    std::cout << maximum_distance << std::endl;

    std::vector<double> alpha{0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2};
    std::vector<unsigned int> indices{1, 2, 3, 4, 5};
    calculator.SetAlpha ( alpha );
    calculator.SetIndices ( indices );
    calculator.SetExp ( [&] (double x) { return exp(x);} );
    calculator.SetLog ( [&] (double x) { return log(x);} );
    calculator.SetPower ( [&] (double x, double y) { return pow(x, y);} );

    auto result = calculator.renyi_entropy_LeonenkoProzanto ( dataset, 2 );
    std::stringstream ss;
    boost::filesystem::path myFile =
        boost::filesystem::current_path() / "myfile.dat";
    boost::filesystem::ofstream ofs ( myFile );
    boost::archive::binary_oarchive oarch ( ofs );
    //oarch << result;
    std::cout << ss.str();

    EXPECT_EQ ( 0, 0 );
}

TEST ( RenyiEntropy, EigenSimpleTest )
{
    unsigned int dimension{10};
    Eigen::MatrixXd dataset = Eigen::MatrixXd::Identity ( dimension, dimension );

    renyi_entropy::renyi_entropy<double> calculator;
    auto tree = calculator.create_KDtree ( dataset );
}

TEST ( RenyiEntropy, Volumes )
{
    unsigned int max_dimension = 10;
    for ( unsigned dimensions = 2; dimensions <= max_dimension; ++dimensions )
    {
        auto volume_unit_sphere =
            renyi_entropy::volume_of_hypersphere<double> ( 1, 2, dimensions );
        auto volume_unit_sphere2 =
            pow ( std::numbers::pi_v<double>, dimensions / 2.0 ) /
            boost::math::tgamma ( dimensions / 2.0 + 1 );

        std::cout << volume_unit_sphere << " " << volume_unit_sphere2 << std::endl;
        // EXPECT_EQ(volume_unit_sphere, volume_unit_sphere2);
    }
}

TEST ( RenyiEntropy, SamplesFromArrays )
{
    unsigned int columns{100}, rows{2};
    Eigen::MatrixXd dataset ( columns, rows );
    for ( unsigned int i = 0; i < columns; ++i )
    {
        for ( unsigned int j = 0; j < rows; ++j )
        {
            dataset ( i, j ) = i + j * 0.5;
        }
    }

    std::map<std::string, std::any> parameters;
    parameters["select_indices"] = std::vector<unsigned int> ( {0, 1, 2, 3, 4, 5} );

    auto prepared_renyi_dataset =
        renyi_entropy::renyi_entropy<double>::samples_from_arrays ( dataset,
                parameters );
    std::cout << prepared_renyi_dataset->transpose() << std::endl;
}

TEST ( RenyiEntropy, CoroutineSamplesFromArrays )
{
    unsigned int columns{100}, rows{2};
    Eigen::MatrixXd dataset ( columns, rows );
    for ( unsigned int i = 0; i < columns; ++i )
    {
        for ( unsigned int j = 0; j < rows; ++j )
        {
            dataset ( i, j ) = i + j * 0.5;
        }
    }

    std::map<std::string, std::any> parameters;
    parameters["select_indices"] = std::vector<unsigned int> ( {0, 1, 2, 3, 4, 5} );

    for ( auto item :
            renyi_entropy::renyi_entropy<double>::generator_samples_from_arrays (
                dataset, parameters ) )
    {
        std::cout << item->transpose() << std::endl;
    }
}

TEST ( RenyiEntropy, PreparesDataset2D )
{
    unsigned int columns{100}, rows{2};
    Eigen::MatrixXd dataset ( columns, rows );
    {
        int count = 0;
        for ( unsigned int i = 0; i < columns; ++i )
        {
            for ( unsigned int j = 0; j < rows; ++j )
            {
                dataset ( i, j ) = count;
                ++count;
            }
        }
    }
    // std::cout << dataset << std::endl;

    {
        auto [dataset1, dataset2] =
            renyi_entropy::renyi_entropy<double>::prepare_dataset ( dataset, false,
                    false, false, 1, 1 );
        for ( unsigned int i = 0; i < columns; ++i )
        {
            EXPECT_EQ ( dataset1 ( i, 0 ), dataset ( i, 0 ) );
            EXPECT_EQ ( dataset2 ( i, 0 ), dataset ( i, 1 ) );
        }
    }

    {
        auto [dataset_swap1, dataset_swap2] =
            renyi_entropy::renyi_entropy<double>::prepare_dataset ( dataset, true,
                    false, false, 1, 1 );
        for ( unsigned int i = 0; i < columns; ++i )
        {
            EXPECT_EQ ( dataset_swap2 ( i, 0 ), dataset ( i, 0 ) );
            EXPECT_EQ ( dataset_swap1 ( i, 0 ), dataset ( i, 1 ) );
        }
    }

    {
        auto [dataset_shuffled1, dataset_shuffled2] =
            renyi_entropy::renyi_entropy<double>::prepare_dataset ( dataset, false,
                    true, false, 1, 1 );
        for ( unsigned int i = 0; i < columns; ++i )
        {
            auto data1 = static_cast<unsigned int> ( dataset_shuffled1 ( i, 0 ) );
            auto data2 = static_cast<unsigned int> ( dataset_shuffled2 ( i, 0 ) );
            EXPECT_EQ ( data1 % 2, 0 );
            EXPECT_EQ ( data2 % 2, 1 );
        }
    }

    {
        auto [dataset_shuffled1, dataset_shuffled2] =
            renyi_entropy::renyi_entropy<double>::prepare_dataset ( dataset, true,
                    true, false, 1, 1 );
        for ( unsigned int i = 0; i < columns; ++i )
        {
            auto data1 = static_cast<unsigned int> ( dataset_shuffled1 ( i, 0 ) );
            auto data2 = static_cast<unsigned int> ( dataset_shuffled2 ( i, 0 ) );
            EXPECT_EQ ( data2 % 2, 0 );
            EXPECT_EQ ( data1 % 2, 1 );
        }
    }
}

TEST ( RenyiEntropy, PreparesDatasetND )
{
    unsigned int maximum_number_rows{10};
    for ( unsigned int number_rows = 3; number_rows < maximum_number_rows;
            ++number_rows )
    {
        unsigned int columns{100}, rows{number_rows};
        Eigen::MatrixXd dataset ( columns, rows );
        {
            int count = 0;
            for ( unsigned int i = 0; i < columns; ++i )
            {
                for ( unsigned int j = 0; j < rows; ++j )
                {
                    dataset ( i, j ) = count;
                    ++count;
                }
            }
        }
        std::cout << dataset << std::endl;

        for ( unsigned int separation_row = 1; separation_row < number_rows - 1;
                ++separation_row )
        {
            {
                auto [dataset1, dataset2] =
                    renyi_entropy::renyi_entropy<double>::prepare_dataset (
                        dataset, false, false, false, separation_row,
                        number_rows - separation_row );
                for ( unsigned int i = 0; i < columns; ++i )
                {
                    for ( unsigned int row_number = 0; row_number < separation_row;
                            ++row_number )
                    {
                        EXPECT_EQ ( dataset1 ( i, row_number ), dataset ( i, row_number ) );
                    }
                    for ( unsigned int row_number = separation_row;
                            row_number < number_rows; ++row_number )
                    {
                        EXPECT_EQ ( dataset2 ( i, row_number - separation_row ),
                                    dataset ( i, row_number ) );
                    }
                }
            }

            {
                auto [dataset_swap1, dataset_swap2] =
                    renyi_entropy::renyi_entropy<double>::prepare_dataset (
                        dataset, true, false, false, separation_row,
                        number_rows - separation_row );
                for ( unsigned int i = 0; i < columns; ++i )
                {
                    for ( unsigned int row_number = 0; row_number < separation_row;
                            ++row_number )
                    {
                        EXPECT_EQ ( dataset_swap2 ( i, row_number ), dataset ( i, row_number ) );
                    }
                    for ( unsigned int row_number = separation_row;
                            row_number < number_rows; ++row_number )
                    {
                        EXPECT_EQ ( dataset_swap1 ( i, row_number - separation_row ),
                                    dataset ( i, row_number ) );
                    }
                }
            }

            {
                auto [dataset_shuffled1, dataset_shuffled2] =
                    renyi_entropy::renyi_entropy<double>::prepare_dataset (
                        dataset, false, true, false, separation_row,
                        number_rows - separation_row );
                for ( unsigned int i = 0; i < columns; ++i )
                {
                    for ( unsigned int row_number = 0; row_number < separation_row;
                            ++row_number )
                    {
                        auto data1 =
                            static_cast<unsigned int> ( dataset_shuffled1 ( i, row_number ) );
                        EXPECT_EQ ( data1 % number_rows, row_number );
                    }
                    for ( unsigned int row_number = separation_row;
                            row_number < number_rows; ++row_number )
                    {
                        auto data2 = static_cast<unsigned int> (
                                         dataset_shuffled2 ( i, row_number - separation_row ) );
                        EXPECT_EQ ( data2 % number_rows, row_number );
                    }
                }
            }

            {
                auto [dataset_shuffled1, dataset_shuffled2] =
                    renyi_entropy::renyi_entropy<double>::prepare_dataset (
                        dataset, true, true, false, separation_row,
                        number_rows - separation_row );
                for ( unsigned int i = 0; i < columns; ++i )
                {
                    for ( unsigned int row_number = 0; row_number < separation_row;
                            ++row_number )
                    {
                        auto data2 =
                            static_cast<unsigned int> ( dataset_shuffled2 ( i, row_number ) );
                        EXPECT_EQ ( data2 % number_rows, row_number );
                    }
                    for ( unsigned int row_number = separation_row;
                            row_number < number_rows; ++row_number )
                    {
                        auto data1 = static_cast<unsigned int> (
                                         dataset_shuffled1 ( i, row_number - separation_row ) );
                        EXPECT_EQ ( data1 % number_rows, row_number );
                    }
                }
            }
        }
    }
}
