
#include <iostream>
#include <algorithm>

#include <gtest/gtest.h>

#include "random_samples.h"

TEST(RandomSamples, SimpleNormalDistribution) {
	Eigen::MatrixXd dataset;

	random_samples::samples_normal_distribution_uncorrelated(dataset, 0, 1, 20, 5);
	std::cout << dataset.transpose()  << std::endl;
}

TEST(RandomSamples, CorrelatedNormalDistribution) {
	Eigen::MatrixXd dataset2;
	Eigen::MatrixXd Sigma(5,5);
	Eigen::VectorXd mean_gaussion_distribution = Eigen::VectorXd::Zero(Sigma.cols());
	Sigma << 
	1, 0.5, 0, 0, 0.5,
	0.5, 1 ,0.5, 0, 0,
	0, 0.5, 1, 0.5, 0,
	0, 0, 0.5, 1, 0.5,
	0.5, 0, 0, 0.5, 1;

	random_samples::samples_normal_distribution_correlated(dataset2, mean_gaussion_distribution, Sigma, 1000); 
	std::cout << dataset2.transpose()  << std::endl;	
}

TEST(RandomSamples, SimpleStudentTlDistribution) {
	double degrees_of_freedom = 2;
	double sigma_gaussion_distribution = 1;
	unsigned int dimension = 5;
	Eigen::MatrixXd dataset;
	Eigen::MatrixXd sigma = Eigen::MatrixXd::Identity ( dimension, dimension ) *
					sigma_gaussion_distribution * sigma_gaussion_distribution;
	Eigen::VectorXd means = Eigen::VectorXd::Zero( dimension );

	random_samples::sample_student_t_distribution(dataset, sigma, means, degrees_of_freedom, 1000);
	//std::cout << dataset.transpose()  << std::endl;
	
	auto random_dataset = dataset.row(0);
	std::sort(std::begin(random_dataset), std::end(random_dataset), std::less<double>());
	std::cout << random_dataset.transpose() << std::endl;
}
