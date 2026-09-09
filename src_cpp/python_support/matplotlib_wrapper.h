#pragma once

#include "generic_python_wrapper.h"


/**
 * @todo Wrapper for matplotlib functions
 */
namespace python_wrappers
{

class matplotlib_wrapper : public generic_python_wrapper
{
public:
  // matplotlib imports
  py::object matplotlib;
  py::object matplotlib_pyplot;
  py::object matplotlib_colormap;
  py::object matplotlib_color;
  py::object pylab;

  // matplotlib_objects
  py::object matplotlib_pyplot_show;
  py::object matplotlib_pyplot_close;
  py::object matplotlib_pyplot_draw;
  py::object matplotlib_pyplot_pause;
  py::object matplotlib_pyplot_figure;
  py::object matplotlib_pyplot_fignum_exists;
  py::object matplotlib_pyplot_plot;
  py::object matplotlib_pyplot_quiver;
  py::object matplotlib_pyplot_contour;
  py::object matplotlib_pyplot_axhline;
  py::object matplotlib_pyplot_axvline;
  py::object matplotlib_pyplot_semilogx;
  py::object matplotlib_pyplot_semilogy;
  py::object matplotlib_pyplot_loglog;
  py::object matplotlib_pyplot_fill;
  py::object matplotlib_pyplot_fill_between;
  py::object matplotlib_pyplot_hist;
  py::object matplotlib_pyplot_scatter;
  py::object matplotlib_pyplot_spy;
  py::object matplotlib_pyplot_subplot;
  py::object matplotlib_pyplot_legend;
  py::object matplotlib_pyplot_xlim;
  py::object matplotlib_pyplot_ylim;
  py::object matplotlib_pyplot_title;
  py::object matplotlib_pyplot_axis;
  py::object matplotlib_pyplot_xlabel;
  py::object matplotlib_pyplot_ylabel;
  py::object matplotlib_pyplot_xticks;
  py::object matplotlib_pyplot_yticks;
  py::object matplotlib_pyplot_xscale;
  py::object matplotlib_pyplot_yscale;
  py::object matplotlib_pyplot_grid;
  py::object matplotlib_pyplot_ion;
  py::object matplotlib_pyplot_ginput;
  py::object matplotlib_pyplot_savefig;
  py::object matplotlib_pyplot_annotate;
  py::object matplotlib_pyplot_clf;
  py::object matplotlib_pyplot_errorbar;
  py::object matplotlib_pyplot_tight_layout;
  py::object matplotlib_pyplot_stem;
  py::object matplotlib_pyplot_xkcd;
  py::object matplotlib_pyplot_suptitle;
  py::object matplotlib_pyplot_bar;
  py::object matplotlib_pyplot_subplots_adjust;
  py::object matplotlib_pyplot_imshow;
  py::object matplotlib_pyplot_colorbar;
  py::object matplotlib_rcParams;
  py::object matplotlib_use;
  py::object matplotlib_get_backend;
  py::object matplotlib_colormaps;

public:
  matplotlib_wrapper();

  ~matplotlib_wrapper();

};

}
