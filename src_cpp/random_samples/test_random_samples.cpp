
#include <algorithm>
#include <iostream>

#include <gtest/gtest.h>

#include "random_samples.h"

TEST ( RandomSamples, SimpleNormalDistribution )
{
    Eigen::MatrixXd dataset;

    random_samples::samples_normal_distribution_uncorrelated ( dataset, 0, 1, 20,
            5 );
    std::cout << dataset.transpose() << std::endl;
}

TEST ( RandomSamples, CorrelatedNormalDistribution )
{
    Eigen::MatrixXd dataset2;
    Eigen::MatrixXd Sigma ( 5, 5 );
    Eigen::VectorXd mean_gaussion_distribution =
        Eigen::VectorXd::Zero ( Sigma.cols() );
    Sigma << 1, 0.5, 0, 0, 0.5, 0.5, 1, 0.5, 0, 0, 0, 0.5, 1, 0.5, 0, 0, 0, 0.5,
          1, 0.5, 0.5, 0, 0, 0.5, 1;

    random_samples::samples_normal_distribution_correlated (
        dataset2, mean_gaussion_distribution, Sigma, 1000 );
    std::cout << dataset2.transpose() << std::endl;
}

TEST ( RandomSamples, SimpleStudentTlDistribution )
{
    double degrees_of_freedom = 2;
    double sigma_gaussion_distribution = 1;
    unsigned int dimension = 5;
    Eigen::MatrixXd dataset;
    Eigen::MatrixXd sigma = Eigen::MatrixXd::Identity ( dimension, dimension ) *
                            sigma_gaussion_distribution *
                            sigma_gaussion_distribution;
    Eigen::VectorXd means = Eigen::VectorXd::Zero ( dimension );

    random_samples::sample_student_t_distribution ( dataset, sigma, means,
            degrees_of_freedom, 1000 );
    // std::cout << dataset.transpose()  << std::endl;

    auto random_dataset = dataset.row ( 0 );
    std::sort ( std::begin ( random_dataset ), std::end ( random_dataset ),
                std::less<double>() );
    std::cout << random_dataset.transpose() << std::endl;
}

TEST ( RandomSamples, AlphaStableDistributionSamples )
{
    unsigned int dimension = 5;
    Eigen::MatrixXd dataset;
    Eigen::MatrixXd Sigma = Eigen::MatrixXd::Identity ( dimension, dimension );
    random_samples::sample_alpha_stable_distribution ( dataset, Sigma, 1.9, 0, 0, 1,
            1000 );
    std::cout << dataset.transpose() << std::endl;
}


TEST ( RandomSamples, AB_process )
{
    std::random_device random{};
    std::mt19937 generator{random() };
    std::normal_distribution distributionA{0.,1.};
    std::normal_distribution distributionB{0.,1.};

    Eigen::MatrixXd dataset;
    random_samples::AB_process(dataset, 0.1, 0.1, 0.1, 100,
        std::make_tuple(
            [&]()
            {
                return distributionA ( generator );
            },
            [&]()
            {
                return distributionB ( generator );
            }
        ),
        std::make_tuple(0.,0.)
    );
    std::cout << dataset << std::endl;
}

TEST ( RandomSamples, CB_process )
{
    std::random_device random{};
    std::mt19937 generator{random() };
    std::normal_distribution distributionC{0.,1.};
    std::normal_distribution distributionB{0.,1.};

    Eigen::MatrixXd dataset;
    random_samples::CB_process(dataset, 0.1, 0.1, 0.1, 100,
        std::make_tuple(
            [&]()
            {
                return distributionC ( generator );
            },
            [&]()
            {
                return distributionB ( generator );
            }
        ),
        std::make_tuple(0.,0.)
    );
    std::cout << dataset << std::endl;
}

TEST ( RandomSamples, ADB_process )
{
    std::random_device random{};
    std::mt19937 generator{random() };
    std::normal_distribution distributionA{0.,1.};
    std::normal_distribution distributionD{0.,1.};
    std::normal_distribution distributionB{0.,1.};

    Eigen::MatrixXd dataset;
    random_samples::ADB_process(dataset, 0.1, 0.1, 0.1, 0.1, 0.1, 100,
        std::make_tuple(
            [&]()
            {
                return distributionA ( generator );
            },
            [&]()
            {
                return distributionD ( generator );
            },
            [&]()
            {
                return distributionB ( generator );
            }
        ),
        std::make_tuple(0., 0., 0.)
    );
    std::cout << dataset << std::endl;
}

TEST ( RandomSamples, ACB_process )
{
    std::random_device random{};
    std::mt19937 generator{random() };
    std::normal_distribution distributionA{0.,1.};
    std::normal_distribution distributionC{0.,1.};
    std::normal_distribution distributionB{0.,1.};

    Eigen::MatrixXd dataset;
    random_samples::ACB_process(dataset, 0.1, 0.1, 0.1, 0.1, 0.1, 100,
        std::make_tuple(
            [&]()
            {
                return distributionA ( generator );
            },
            [&]()
            {
                return distributionC ( generator );
            },
            [&]()
            {
                return distributionB ( generator );
            }
        ),
        std::make_tuple(0., 0., 0.)
    );
    std::cout << dataset << std::endl;
}
