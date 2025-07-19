
#include "random_samples.h"

#include <cmath>
#include <exception>
#include <list>
#include <random>

#include <boost/math/special_functions/gamma.hpp>

#include <eigen3/Eigen/Eigenvalues>

extern "C" {
#include "stable.h"
}

namespace random_samples
{

double full_correlation_matrix ( int dimension, double q )
{
    return pow ( 1 - q, dimension - 1 ) * ( 1 + ( dimension - 1 ) * q );
}

double tridiagonal_matrix_determinant ( int dimension, double q )
{
    if ( dimension == 1 )
    {
        return 1;
    }
    else if ( dimension == 2 )
    {
        return 1 - pow ( q, 2 );
    }
    else
    {
        std::vector<double> sample{1, 1 - pow ( q, 2 ) };
        for ( int n{2}; n < dimension; ++n )
        {
            auto result = sample[n - 1] - pow ( q, 2 ) * sample[n - 2];
            sample.push_back ( result );
        }
        return sample[dimension - 1];
    }
}

double Renyi_entropy_normal_distribution ( const double q,
        const Eigen::MatrixXd &Sigma )
{
    if ( Sigma.rows() == Sigma.cols() )
    {
        auto m = Sigma.cols();
        if ( q != 1 )
        {
            return log ( 2 * std::numbers::pi_v<double> ) * m / 2 +
                   log ( Sigma.eval().determinant() ) / 2 - m * log ( q ) / ( 1 - q ) / 2;
        }
        else
        {
            return log ( 2 * std::numbers::pi_v<double> * std::exp ( 1 ) ) * m / 2 +
                   log ( Sigma.eval().determinant() ) / 2;
        }
    }
    else
    {
        std::cerr << "Nonsymmetric matrix" << std::endl;
        return 0;
    }
}

void samples_normal_distribution_uncorrelated (
    std::vector<std::vector<double>> &dataset,
    double mean_gaussion_distribution, double sigma_gaussion_distribution,
    unsigned int number_samples, unsigned int dimension )
{
    std::random_device random{};
    std::mt19937 generator{random() };
    std::normal_distribution distribution{mean_gaussion_distribution,
                                          sigma_gaussion_distribution};

    auto random_double = [&distribution, &generator]
    {
        return distribution ( generator );
    };

    for ( auto n{number_samples}; n; --n )
    {
        dataset.push_back ( std::vector<double>() );
        for ( int i :
                std::ranges::iota_view<unsigned int, unsigned int> ( 0, dimension ) )
        {
            dataset.back().push_back ( random_double() );
        }
    }
}

void samples_normal_distribution_uncorrelated (
    Eigen::MatrixXd &dataset, double mean_gaussion_distribution,
    double sigma_gaussion_distribution, unsigned int number_samples,
    unsigned int dimension )
{
    std::random_device random{};
    std::mt19937 generator{random() };
    std::normal_distribution distribution{mean_gaussion_distribution,
                                          sigma_gaussion_distribution};

    auto random_double = [&distribution, &generator]
    {
        return distribution ( generator );
    };

    dataset.conservativeResize ( dimension, number_samples );
    for ( auto n{0}; n < number_samples; ++n )
    {
        Eigen::VectorXd random_vector ( dimension );
        for ( int i :
                std::ranges::iota_view<unsigned int, unsigned int> ( 0, dimension ) )
        {
            random_vector[i] = random_double();
        }
        dataset.col ( n ) = random_vector;
    }
}

void samples_normal_distribution_correlated (
    Eigen::MatrixXd &dataset, const Eigen::VectorXd &mean_gaussion_distribution,
    const Eigen::MatrixXd &Sigma, unsigned int number_samples )
{
    if ( Sigma.cols() == Sigma.rows() )
    {
        auto dimension = Sigma.cols();

        Eigen::MatrixXd normal_dataset;

        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> eigenValueSolver (
            Sigma, Eigen::MatrixXd::Identity ( Sigma.cols(), Sigma.cols() ) );
        const auto eigenvalueMatrix = eigenValueSolver.eigenvalues();
        const auto eigenvectorMatrix = eigenValueSolver.eigenvectors();
        const auto diagonal_matrix =
            eigenvalueMatrix.array().sqrt().matrix().asDiagonal();
        const auto transform = eigenvectorMatrix * diagonal_matrix;
        // std::cout << eigenvalueMatrix << std::endl;

        samples_normal_distribution_uncorrelated ( normal_dataset, 0., 1.,
                number_samples, dimension );
        // std::cout << transform << " " << eigenvalueMatrix << " " <<
        // eigenvectorMatrix << std::endl;
        dataset =
            ( transform * normal_dataset ).colwise() + mean_gaussion_distribution;
        // dataset = dataset.colwise() + mean_gaussion_distribution;
    }
    else
    {
        // throw Error("sigma parameter has wrong type");
    }
}

void sample_student_t_distribution ( Eigen::MatrixXd &dataset,
                                     const Eigen::MatrixXd &Sigma,
                                     const Eigen::VectorXd &mean,
                                     const double degrees_of_freedom,
                                     unsigned int number_samples )
{
    std::random_device random{};
    std::mt19937 generator{random() };
    std::chi_squared_distribution<double> chi2_generator{degrees_of_freedom};

    auto random_double = [&chi2_generator, &generator]
    {
        return chi2_generator ( generator );
    };

    Eigen::MatrixXd dataset_normal;
    samples_normal_distribution_correlated ( dataset_normal, mean, Sigma,
            number_samples );

    auto dimension = Sigma.cols();
    dataset.conservativeResize ( dimension, number_samples );
    for ( unsigned int i{0}; i < number_samples; ++i )
    {
        auto random_chi2_number = random_double();
        auto random_vector =
            mean +
            dataset_normal.col ( i ) * ( sqrt ( degrees_of_freedom / random_chi2_number ) );
        dataset.col ( i ) = random_vector;
    }
    // std::cout << dataset.transpose() << std::endl;
}

double Renyi_student_t_distribution ( double degrees_of_freedom, const double q,
                                      const Eigen::MatrixXd &sigma )
{
    double dimension = sigma.cols();
    if ( q * ( degrees_of_freedom + dimension ) / 2 - dimension / 2 <= 0 )
    {
        return std::nan ( "" );
    }

    if ( q == 1 )
    {
        return 0;
    }
    else
    {
        auto arg = pow ( sigma.eval().determinant(), ( 1.0 - q ) / 2 );

        auto arg2 = pow ( degrees_of_freedom * std::numbers::pi_v<double>,
                          dimension / 2 * ( 1 - q ) );
        auto arg31 =
            pow ( boost::math::tgamma ( ( degrees_of_freedom + dimension ) / 2 ), q );
        auto arg32 = boost::math::tgamma (
                         ( q * ( degrees_of_freedom + dimension ) - dimension ) / 2 );
        auto arg33 = boost::math::tgamma ( q * ( degrees_of_freedom + dimension ) / 2 );
        auto arg34 = pow ( boost::math::tgamma ( degrees_of_freedom / 2 ), q );
        auto arg3 = ( arg31 * arg32 ) / ( arg33 * arg34 );

        return 1 / ( 1 - q ) * log2 ( arg3 * arg * arg2 );
    }
}

void sample_alpha_stable_distribution ( Eigen::MatrixXd &dataset,
                                        const Eigen::MatrixXd &Sigma,
                                        const double alpha, const double beta,
                                        const double mu, const double param,
                                        unsigned int number_samples )
{
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MyMatrix;
    double sigma = 1;
    const unsigned int dimension = Sigma.cols();
    std::vector<StableDist *> dist ( dimension );
    std::fill ( dist.begin(), dist.end(), nullptr );

    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> eigenValueSolver (
        Sigma, Eigen::MatrixXd::Identity ( Sigma.cols(), Sigma.cols() ) );
    const auto eigenvalueMatrix = eigenValueSolver.eigenvalues();
    const auto eigenvectorMatrix = eigenValueSolver.eigenvectors();

    std::vector<double> raw_dataset ( dimension * number_samples );
    dataset.conservativeResize ( dimension, number_samples );

    int i = 0;
    for ( auto &item : dist )
    {
        if ( ( item = stable_create ( alpha, beta, eigenvalueMatrix[i], mu, param ) ) ==
                nullptr )
        {
            printf ( "Error while creatring ditribution. Please check introduced data" );
            exit ( 1 );
        }
        ++i;
    }

    for ( unsigned int i{0}; i < dimension * number_samples; ++i )
    {
        unsigned int column = i % dimension;
        stable_rnd ( dist[column], raw_dataset.data(), dimension * number_samples );
    }
    dataset = Eigen::Map<MyMatrix> ( raw_dataset.data(), dimension, number_samples );

    for ( auto &item : dist )
    {
        stable_free ( item );
    }
}

// A_{t+1} \ = \ a A_{t} \ + \ \varepsilon^{A}_t
// B^I_{t+1} \ = \ b B^I_{t} \ + \  \eta A_{t} \ + \ \varepsilon^{B^I}_t
void AB_process(Eigen::MatrixXd &dataset, double alpha, double beta, double eta, unsigned int number_samples, std::tuple<std::function<double ()>, std::function<double ()>> random_generator, std::tuple<double, double> initial_confition)
{
    constexpr unsigned int A = 0;
    constexpr unsigned int B = 1;
    dataset.conservativeResize ( 2, number_samples );
    dataset(A, 0) = std::get<A>(initial_confition);
    dataset(B, 0) = std::get<B>(initial_confition);
    for (unsigned int sample = 1; sample < number_samples; ++ sample)
    {
        auto epsilon_A = std::get<A>(random_generator)();
        auto epsilon_B = std::get<B>(random_generator)();
        dataset(A, sample) = alpha * dataset(A, sample - 1) + epsilon_A;
        dataset(B, sample) = beta * dataset(B, sample - 1) + eta * dataset(A, sample - 1) + epsilon_B;
    }
}

// C_{t+1} \ = \  a C_t \ + \ \varepsilon^C_t
// B^{II}_{t+1} \ = \  bB^{II}_t \ + \  eta C^3_t \ + \  \varepsilon^{B^{II}}_t
void CB_process(Eigen::MatrixXd &dataset, double alpha, double beta, double eta, unsigned int number_samples, std::tuple<std::function<double ()>, std::function<double ()>> random_generator, std::tuple<double, double> initial_confition, double lambda)
{
    constexpr unsigned int C = 0;
    constexpr unsigned int B = 1;
    dataset.conservativeResize ( 2, number_samples );
    dataset(C, 0) = std::get<C>(initial_confition);
    dataset(B, 0) = std::get<B>(initial_confition);
    for (unsigned int sample = 1; sample < number_samples; ++ sample)
    {
        auto epsilon_C = std::get<C>(random_generator)();
        auto epsilon_B = std::get<B>(random_generator)();
        dataset(C, sample) = alpha * dataset(C, sample - 1) + epsilon_C;
        dataset(B, sample) = beta * dataset(B, sample - 1) + eta * signbit(dataset(C, sample - 1)) * pow(abs(dataset(C, sample - 1)), lambda) + epsilon_B;
    }
}

// A_{t+1}\ = \ aA_{t} \ + \ \varepsilon^A_t
// D_{t+1} \ = \ aD_{t} \ + \  \varepsilon^D_t
// B^{III}_{t+1} \ = \ bB^{III}_{t} \ + \ \eta_1 A_{t} \ + \  \eta_2 D_{t} \ + \  \varepsilon^{B^{III}}_t
void ADB_process(Eigen::MatrixXd &dataset, double alpha, double beta, double gamma, double eta1, double eta2, unsigned int number_samples, std::tuple<std::function<double ()>, std::function<double ()>, std::function<double ()>> random_generator, std::tuple<double, double, double> initial_confition)
{
    constexpr unsigned int A = 0;
    constexpr unsigned int D = 1;
    constexpr unsigned int B = 2;
    dataset.conservativeResize ( 3, number_samples );
    dataset(A, 0) = std::get<A>(initial_confition);
    dataset(D, 0) = std::get<D>(initial_confition);
    dataset(B, 0) = std::get<B>(initial_confition);
    for (unsigned int sample = 1; sample < number_samples; ++ sample)
    {
        auto epsilon_A = std::get<A>(random_generator)();
        auto epsilon_D = std::get<D>(random_generator)();
        auto epsilon_B = std::get<B>(random_generator)();
        dataset(A, sample) = alpha * dataset(A, sample - 1) + epsilon_A;
        dataset(D, sample) = beta * dataset(D, sample - 1) + epsilon_D;
        dataset(B, sample) = gamma * dataset(B, sample - 1) + eta1 * dataset(A, sample - 1) + eta2 * dataset(D, sample - 1) + epsilon_B;
    }
}

// A_{t+1} \ = \  aA_{t} \ + \  \varepsilon^A_t
// C_{t+1} \ = \  cC_{t} \ + \  \varepsilon^C_t
// B^{IV}_{t+1}  \ = \ bB^{IV}_{t} \ + \  \eta A_{t} \ + \  \kappa C^3_{t} \ + \ \varepsilon^{B^{IV}}_t
void ACB_process(Eigen::MatrixXd &dataset, double alpha, double beta, double gamma, double eta1, double eta2, unsigned int number_samples, std::tuple<std::function<double ()>, std::function<double ()>, std::function<double ()>> random_generator, std::tuple<double, double, double> initial_confition, double lambda)
{
    constexpr unsigned int A = 0;
    constexpr unsigned int C = 1;
    constexpr unsigned int B = 2;
    dataset.conservativeResize ( 3, number_samples );
    dataset(A, 0) = std::get<A>(initial_confition);
    dataset(C, 0) = std::get<C>(initial_confition);
    dataset(B, 0) = std::get<B>(initial_confition);
    for (unsigned int sample = 1; sample < number_samples; ++ sample)
    {
        auto epsilon_A = std::get<A>(random_generator)();
        auto epsilon_C = std::get<C>(random_generator)();
        auto epsilon_B = std::get<B>(random_generator)();
        dataset(A, sample) = alpha * dataset(A, sample - 1) + epsilon_A;
        dataset(C, sample) = beta * dataset(C, sample - 1) + epsilon_C;
        dataset(B, sample) = gamma * dataset(B, sample - 1) + eta1 * dataset(A, sample - 1) + eta2 * signbit(dataset(C, sample - 1)) * pow(abs(dataset(C, sample - 1)), lambda) + epsilon_B;
    }
}


} // namespace random_samples
