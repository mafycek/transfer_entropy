#pragma once

#include "generic_python_wrapper.h"

namespace python_wrappers
{

/**
 * @todo Wrapper of pandas
 */
class pandas_wrapper : public generic_python_wrapper
{
public:
  py::object pandas;
  py::object pandas_dataframe;
  py::object pandas_multiindex;
  py::object pandas_multiindex_from_tuples;
  py::object pandas_series;
  py::object pandas_concat;

  pandas_wrapper();

  ~pandas_wrapper();

};

}
