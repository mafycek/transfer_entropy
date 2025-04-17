//
// Created by hynek on 5.3.25.
//

#include <eigen3/Eigen/Dense>

#include "renyi_entropy.h"

namespace renyi_entropy {

double Renyi_entropy_normal_distribution(double q,
                                         const Eigen::MatrixXd &Sigma) {
  if (Sigma.rows() == Sigma.cols()) {
    auto m = Sigma.cols();
    if (q != 1) {
      return log(2 * std::numbers::pi_v<double>) * m / 2 +
             log(Sigma.eval().determinant()) / 2 - m * log(q) / (1 - q) / 2;
    } else {
      return log(2 * std::numbers::pi_v<double> * std::exp(1)) * m / 2 +
             log(Sigma.eval().determinant()) / 2;
    }
  } else {
    std::cerr << "Nonsymmetric matrix" << std::endl;
    return 0;
  }
}

} // namespace renyi_entropy
