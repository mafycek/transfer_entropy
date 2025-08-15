
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
#include "matplotlib_wrapper.h"
#pragma GCC diagnostic pop

namespace python_wrappers
{

matplotlib_wrapper::matplotlib_wrapper()
: generic_python_wrapper()
{
  matplotlib = py::module_::import("matplotlib");
  matplotlib_pyplot = py::module_::import("matplotlib.pyplot");
  matplotlib_colormap = py::module_::import("matplotlib.cm");
  matplotlib_color = py::module_::import("matplotlib.colors");
  pylab = py::module_::import("pylab");

  matplotlib_pyplot_show = matplotlib_pyplot.attr("show");
  matplotlib_pyplot_close = matplotlib_pyplot.attr("close");
  matplotlib_pyplot_draw = matplotlib_pyplot.attr("draw");
  matplotlib_pyplot_pause = matplotlib_pyplot.attr("pause");
  matplotlib_pyplot_figure = matplotlib_pyplot.attr("figure");
  matplotlib_pyplot_fignum_exists = matplotlib_pyplot.attr("fignum_exists");
  matplotlib_pyplot_plot = matplotlib_pyplot.attr("plot");
  matplotlib_pyplot_quiver = matplotlib_pyplot.attr("quiver");
  matplotlib_pyplot_contour = matplotlib_pyplot.attr("contour");
  matplotlib_pyplot_axhline = matplotlib_pyplot.attr("axhline");
  matplotlib_pyplot_axvline = matplotlib_pyplot.attr("axvline");
  matplotlib_pyplot_semilogx = matplotlib_pyplot.attr("semilogx");
  matplotlib_pyplot_semilogy = matplotlib_pyplot.attr("semilogy");
  matplotlib_pyplot_loglog = matplotlib_pyplot.attr("loglog");
  matplotlib_pyplot_fill = matplotlib_pyplot.attr("fill");
  matplotlib_pyplot_fill_between = matplotlib_pyplot.attr("fill_between");
  matplotlib_pyplot_hist = matplotlib_pyplot.attr("hist");
  matplotlib_pyplot_scatter = matplotlib_pyplot.attr("scatter");
  matplotlib_pyplot_spy = matplotlib_pyplot.attr("spy");
  matplotlib_pyplot_subplot = matplotlib_pyplot.attr("subplot");
  matplotlib_pyplot_legend = matplotlib_pyplot.attr("legend");
  matplotlib_pyplot_xlim = matplotlib_pyplot.attr("xlim");
  matplotlib_pyplot_ylim = matplotlib_pyplot.attr("ylim");
  matplotlib_pyplot_title = matplotlib_pyplot.attr("title");
  matplotlib_pyplot_axis = matplotlib_pyplot.attr("axis");
  matplotlib_pyplot_xlabel = matplotlib_pyplot.attr("xlabel");
  matplotlib_pyplot_ylabel = matplotlib_pyplot.attr("ylabel");
  matplotlib_pyplot_xticks = matplotlib_pyplot.attr("xticks");
  matplotlib_pyplot_yticks = matplotlib_pyplot.attr("yticks");
  matplotlib_pyplot_xscale = matplotlib_pyplot.attr("xscale");
  matplotlib_pyplot_yscale = matplotlib_pyplot.attr("yscale");
  matplotlib_pyplot_grid = matplotlib_pyplot.attr("grid");

  matplotlib_pyplot_ion = matplotlib_pyplot.attr("ion");
  matplotlib_pyplot_ginput = matplotlib_pyplot.attr("ginput");
  matplotlib_pyplot_savefig = matplotlib_pyplot.attr("savefig");
  matplotlib_pyplot_annotate = matplotlib_pyplot.attr("annotate");
  matplotlib_pyplot_clf = matplotlib_pyplot.attr("clf");
  matplotlib_pyplot_errorbar = matplotlib_pyplot.attr("errorbar");
  matplotlib_pyplot_tight_layout = matplotlib_pyplot.attr("tight_layout");
  matplotlib_pyplot_stem = matplotlib_pyplot.attr("stem");
  matplotlib_pyplot_xkcd = matplotlib_pyplot.attr("xkcd");
  matplotlib_pyplot_suptitle = matplotlib_pyplot.attr("suptitle");
  matplotlib_pyplot_bar = matplotlib_pyplot.attr("bar");
  matplotlib_pyplot_subplots_adjust = matplotlib_pyplot.attr("subplots_adjust");
  matplotlib_pyplot_imshow = matplotlib_pyplot.attr("imshow");
  matplotlib_pyplot_colorbar = matplotlib_pyplot.attr("colorbar");
  matplotlib_rcParams = matplotlib.attr("rcParams");
  matplotlib_use = matplotlib.attr("use");
  matplotlib_get_backend = matplotlib.attr("get_backend");
  matplotlib_colormaps = matplotlib.attr("colormaps");
}

matplotlib_wrapper::~matplotlib_wrapper()
{
}

}
