
#include <algorithm>
#include <iostream>

#include <gtest/gtest.h>

#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MAKE_OPAQUE(std::map<std::tuple<std::string, std::string>, double>)

TEST ( PyBind11, matplotlib )
{
  Py_Initialize();

  if ( !Py_IsInitialized() )
  {
    throw std::runtime_error("Unable to initialize Python interpreter!");
  }
  else
  {
    std::cout << Py_GetVersion() << std::endl;
    std::cout << Py_GetPlatform() << std::endl;
    std::cout << Py_GetCopyright() << std::endl;
    std::cout << Py_GetCompiler() << std::endl;
    std::cout << Py_GetBuildInfo() << std::endl;
  }

  py::object matplotlib = py::module_::import("matplotlib");
  py::object matplotlib_pyplot = py::module_::import("matplotlib.pyplot");
  py::object matplotlib_colormap = py::module_::import("matplotlib.cm");
  py::object matplotlib_color = py::module_::import("matplotlib.colors");
  py::object pylab = py::module_::import("pylab");

  py::object matplotlib_pyplot_show = matplotlib_pyplot.attr("show");
  py::object matplotlib_pyplot_close = matplotlib_pyplot.attr("close");
  py::object matplotlib_pyplot_draw = matplotlib_pyplot.attr("draw");
  py::object matplotlib_pyplot_pause = matplotlib_pyplot.attr("pause");
  py::object matplotlib_pyplot_figure = matplotlib_pyplot.attr("figure");
  py::object matplotlib_pyplot_fignum_exists = matplotlib_pyplot.attr("fignum_exists");
  py::object matplotlib_pyplot_plot = matplotlib_pyplot.attr("plot");
  py::object matplotlib_pyplot_quiver = matplotlib_pyplot.attr("quiver");
  py::object matplotlib_pyplot_contour = matplotlib_pyplot.attr("contour");
  py::object matplotlib_pyplot_axhline = matplotlib_pyplot.attr("axhline");
  py::object matplotlib_pyplot_axvline = matplotlib_pyplot.attr("axvline");
  py::object matplotlib_pyplot_semilogx = matplotlib_pyplot.attr("semilogx");
  py::object matplotlib_pyplot_semilogy = matplotlib_pyplot.attr("semilogy");
  py::object matplotlib_pyplot_loglog = matplotlib_pyplot.attr("loglog");
  py::object matplotlib_pyplot_fill = matplotlib_pyplot.attr("fill");
  py::object matplotlib_pyplot_fill_between = matplotlib_pyplot.attr("fill_between");
  py::object matplotlib_pyplot_hist = matplotlib_pyplot.attr("hist");
  py::object matplotlib_pyplot_scatter = matplotlib_pyplot.attr("scatter");
  py::object matplotlib_pyplot_spy = matplotlib_pyplot.attr("spy");
  py::object matplotlib_pyplot_subplot = matplotlib_pyplot.attr("subplot");
  py::object matplotlib_pyplot_legend = matplotlib_pyplot.attr("legend");
  py::object matplotlib_pyplot_xlim = matplotlib_pyplot.attr("xlim");
  py::object matplotlib_pyplot_ylim = matplotlib_pyplot.attr("ylim");
  py::object matplotlib_pyplot_title = matplotlib_pyplot.attr("title");
  py::object matplotlib_pyplot_axis = matplotlib_pyplot.attr("axis");
  py::object matplotlib_pyplot_xlabel = matplotlib_pyplot.attr("xlabel");
  py::object matplotlib_pyplot_ylabel = matplotlib_pyplot.attr("ylabel");
  py::object matplotlib_pyplot_xticks = matplotlib_pyplot.attr("xticks");
  py::object matplotlib_pyplot_yticks = matplotlib_pyplot.attr("yticks");
  py::object matplotlib_pyplot_xscale = matplotlib_pyplot.attr("xscale");
  py::object matplotlib_pyplot_yscale = matplotlib_pyplot.attr("yscale");
  py::object matplotlib_pyplot_grid = matplotlib_pyplot.attr("grid");

  py::object matplotlib_pyplot_ion = matplotlib_pyplot.attr("ion");
  py::object matplotlib_pyplot_ginput = matplotlib_pyplot.attr("ginput");
  py::object matplotlib_pyplot_savefig = matplotlib_pyplot.attr("savefig");
  py::object matplotlib_pyplot_annotate = matplotlib_pyplot.attr("annotate");
  py::object matplotlib_pyplot_clf = matplotlib_pyplot.attr("clf");
  py::object matplotlib_pyplot_errorbar = matplotlib_pyplot.attr("errorbar");
  py::object matplotlib_pyplot_tight_layout = matplotlib_pyplot.attr("tight_layout");
  py::object matplotlib_pyplot_stem = matplotlib_pyplot.attr("stem");
  py::object matplotlib_pyplot_xkcd = matplotlib_pyplot.attr("xkcd");
  py::object matplotlib_pyplot_suptitle = matplotlib_pyplot.attr("suptitle");
  py::object matplotlib_pyplot_bar = matplotlib_pyplot.attr("bar");
  py::object matplotlib_pyplot_subplots_adjust = matplotlib_pyplot.attr("subplots_adjust");
  py::object matplotlib_pyplot_imshow = matplotlib_pyplot.attr("imshow");
  py::object matplotlib_pyplot_colorbar = matplotlib_pyplot.attr("colorbar");
  py::object matplotlib_pyplot_rcParams = matplotlib_pyplot.attr("rcParams");
  py::object matplotlib_use = matplotlib.attr("use");
  py::object matplotlib_get_backend = matplotlib.attr("get_backend");
  py::object matplotlib_colormaps = matplotlib.attr("colormaps");

  std::vector<double> x{0, 1, 2, 3, 4};
  std::vector<double> y{0, 1, 5, 6, 1};
  std::vector<double> y1{0, 1, 5, 6, 1};
  std::vector<double> y2{0, 0.5, 3, 2, 1};
  std::vector<double> z1{0, 0.7, 4.5, 5.5, 1};
  std::vector<double> z2{0, 0.6, 3.5, 2.5, 1};

  auto text = matplotlib_pyplot_rcParams.attr("__setitem__");
  text("text.usetex", py::bool_(true));
  py::print(matplotlib_pyplot_rcParams["text.usetex"]);
  matplotlib_pyplot_rcParams["text.usetex"] = py::bool_(false);
  py::print(matplotlib_pyplot_rcParams["text.usetex"]);
  std::string title{"{\\Large \\textrm{RTE\\ dependence\\ on\\ }$\\alpha$ }"};
  matplotlib_pyplot_title(py::str(title));
  matplotlib_pyplot_plot(x, y);
  matplotlib_pyplot_show();

  auto fill_colormap = matplotlib_colormaps [py::str("jet")];
  auto color = fill_colormap( py::float_(0.1 ));
  py::dict kwargs ("color"_a = color, "linewidth"_a = py::float_(3.), "alpha"_a = py::float_(0.7));
  matplotlib_pyplot_fill_between(x, y1, y2, **kwargs);

  auto color2 = fill_colormap( py::float_(0.9 ));
  kwargs[py::str("color")] = color2;
  matplotlib_pyplot_fill_between(x, z1, z2, **kwargs);

  matplotlib_pyplot_show();
  matplotlib_pyplot_close();
}

TEST ( PyBind11, pandas )
{
  Py_Initialize();

  if ( !Py_IsInitialized() )
  {
    throw std::runtime_error("Unable to initialize Python interpreter!");
  }
  else
  {
    std::cout << Py_GetVersion() << std::endl;
    std::cout << Py_GetPlatform() << std::endl;
    std::cout << Py_GetCopyright() << std::endl;
    std::cout << Py_GetCompiler() << std::endl;
    std::cout << Py_GetBuildInfo() << std::endl;
  }

  py::object pandas = py::module_::import("pandas");
  py::object pandas_dataframe = pandas.attr("DataFrame");
  py::object pandas_multiindex = pandas.attr("MultiIndex");
  py::object pandas_multiindex_from_tuples = pandas_multiindex.attr("from_tuples");

  auto column_multiindex = pandas_multiindex_from_tuples("tuples"_a=std::list<py::tuple>{py::make_tuple(py::str("A"), py::make_tuple(0, 1, 2)), py::make_tuple(py::str("A"), py::make_tuple(0, 1, 2, 3)), py::make_tuple(py::str("B"), py::make_tuple(0, 1, 2, 3, 4))}, "names"_a=std::list<py::str>{py::str("AA"), py::str("BB")});
  auto indexn_multiindex = pandas_multiindex_from_tuples("tuples"_a=std::list<py::tuple>{py::make_tuple(py::str("A"), py::make_tuple(0, 1)), py::make_tuple(py::str("A"), py::make_tuple(0)), py::make_tuple(py::str("B"), py::make_tuple(0, 1, 2))}, "names"_a=std::list<py::str>{py::str("aa"), py::str("bb")});
  auto empty_dataframe = pandas_dataframe("index"_a=indexn_multiindex, "columns"_a=column_multiindex);
  py::print(empty_dataframe);

  std::map<py::tuple, py::float_ > column_dataset{{py::make_tuple("a", "b"), py::float_(1.)}};
  py::object dataset = py::cast( column_dataset );
  py::print(dataset);
  py::dict column;
  column [ py::make_tuple( 'a', 'v' ) ] = py::cast( column_dataset );
  py::object dataframe = pandas_dataframe();
  py::object dataframe_join = dataframe.attr("join");
  py::print(dataframe_join);
  dataframe_join(column);
  py::print(dataframe);

}

TEST ( PyBind11, vectorTest )
{
  Py_Initialize();

  if ( !Py_IsInitialized() )
  {
    throw std::runtime_error("Unable to initialize Python interpreter!");
  }
  else
  {
    std::cout << Py_GetVersion() << std::endl;
    std::cout << Py_GetPlatform() << std::endl;
    std::cout << Py_GetCopyright() << std::endl;
    std::cout << Py_GetCompiler() << std::endl;
    std::cout << Py_GetBuildInfo() << std::endl;
  }
  std::list<unsigned int> dataset;
  dataset.push_back(1);
  dataset.push_back(3);
  dataset.push_back(6);
  auto tuple = py::tuple(py::cast(dataset));
  py::print(tuple);

}
