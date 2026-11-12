
#include "FunctionTableCache.h"

#include <boost/math/special_functions/gamma.hpp>
#include <gtest/gtest.h>

TEST(SampleStore, gamma) {
  FunctionTableCache<double, double> gamma(
      [&](double x) -> double { return boost::math::tgamma<double>(x); });
  auto g1 = gamma(1);
  auto g2 = gamma(2);
  auto g_pi = gamma(0.5);
  auto g3 = gamma(1);
}
