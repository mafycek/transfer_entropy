#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstdint>

#include "random_samples.h"

namespace py = pybind11;

int add(int a, int b) { return a + b; }

PYBIND11_MAKE_OPAQUE(std::vector<double>)

py::array_t<double> pyA_process(py::numpy_scalar<std::uint32_t> length,
                                py::numpy_scalar<double> rho1,
                                py::numpy_scalar<double> rho2,
                                py::numpy_scalar<double> variance,
                                py::numpy_scalar<double> alpha,
                                py::numpy_scalar<double> skewness_parameter,
                                py::numpy_scalar<std::uint32_t> cut_off_length,
                                py::numpy_scalar<std::uint32_t> skip,
                                py::numpy_scalar<std::uint32_t> seed) {
  auto result = py::array_t<double>(length.value);
  py::buffer_info buffer = result.request();
  double *result_ptr = static_cast<double *>(buffer.ptr);

  int index = 0;
  for (auto &&item : random_samples::A_process<double>(
           length.value, rho1.value, rho2.value, variance.value, alpha.value,
           skewness_parameter.value, cut_off_length.value, skip.value,
           seed.value)) {
    result_ptr[index] = item;

    ++index;
  }
  return result;
}

py::array_t<double>
pyB_process(py::numpy_scalar<std::uint32_t> length,
            py::numpy_scalar<double> rho1, py::numpy_scalar<double> lambda,
            py::numpy_scalar<double> rho2, py::numpy_scalar<double> rho3,
            py::numpy_scalar<double> variance, py::numpy_scalar<double> alpha,
            py::numpy_scalar<double> skewness_parameter,
            py::numpy_scalar<std::uint32_t> cut_off_length,
            py::numpy_scalar<std::uint32_t> skip,
            py::numpy_scalar<std::uint32_t> seed) {
  auto result = py::array_t<double>(length.value);
  py::buffer_info buffer = result.request();
  double *result_ptr = static_cast<double *>(buffer.ptr);

  int index = 0;
  for (auto &&item : random_samples::generalized_B_process<double>(
           length.value, rho1.value, lambda.value, rho2.value, rho3.value,
           variance.value, alpha.value, skewness_parameter.value,
           cut_off_length.value, skip.value, seed.value)) {
    result_ptr[index] = item;
    ++index;
  }
  return result;
}

PYBIND11_MODULE(pyrandom_samples, py_module) {
  py_module.doc() = "Random sample generator";
  py_module.def("add", &add, "Add two integers", py::arg("i"), py::arg("j"));

  py_module.def("A_process", &pyA_process,
                py::arg("length") = py::numpy_scalar<std::uint32_t>(1000),
                py::arg("rho1") = py::numpy_scalar<double>(0.0),
                py::arg("rho2") = py::numpy_scalar<double>(0.5),
                py::arg("variance") = py::numpy_scalar<double>(1.0),
                py::arg("alpha") = py::numpy_scalar<double>(2.0),
                py::arg("skewness_parameter") = py::numpy_scalar<double>(0.0),
                py::arg("cut_off_length") =
                    py::numpy_scalar<std::uint32_t>(1000U),
                py::arg("skip") = py::numpy_scalar<std::uint32_t>(0U),
                py::arg("seed") = py::numpy_scalar<std::uint32_t>(0U),
                "A process sampler");
  py_module.def("B_process", &pyB_process,
                py::arg("length") = py::numpy_scalar<std::uint32_t>(1000),
                py::arg("rho1") = py::numpy_scalar<double>(0.0),
                py::arg("lambda") = py::numpy_scalar<double>(0.0),
                py::arg("rho2") = py::numpy_scalar<double>(0.5),
                py::arg("rho3") = py::numpy_scalar<double>(0.0),
                py::arg("variance") = py::numpy_scalar<double>(1.0),
                py::arg("alpha") = py::numpy_scalar<double>(2.0),
                py::arg("skewness_parameter") = py::numpy_scalar<double>(0.0),
                py::arg("cut_off_length") =
                    py::numpy_scalar<std::uint32_t>(1000U),
                py::arg("skip") = py::numpy_scalar<std::uint32_t>(0U),
                py::arg("seed") = py::numpy_scalar<std::uint32_t>(0U),
                "generalized B process generator");
}
