
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
#include "pandas_wrapper.h"
#pragma GCC diagnostic pop

namespace python_wrappers
{

pandas_wrapper::pandas_wrapper()
: generic_python_wrapper()
{
  pandas = py::module_::import("pandas");

  pandas_dataframe = pandas.attr("DataFrame");
  pandas_multiindex = pandas.attr("MultiIndex");
  pandas_multiindex_from_tuples = pandas_multiindex.attr("from_tuples");
  pandas_series = pandas.attr("Series");
  pandas_concat = pandas.attr("concat");
}

pandas_wrapper::~pandas_wrapper()
{
}

}
