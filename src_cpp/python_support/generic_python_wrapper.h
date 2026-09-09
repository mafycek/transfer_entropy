#pragma once

#include <iostream>

#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;
using namespace pybind11::literals;

namespace python_wrappers
{

/**
 * @todo write docs
 */
class generic_python_wrapper
{
protected:
  static bool _initialized;

public:
  generic_python_wrapper();

  ~generic_python_wrapper();
};

}
