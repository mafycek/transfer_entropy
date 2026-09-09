
#include <iostream>
#include <string>
#include <algorithm>
#include <array>
#include <cstdint> // <cstdint> requires c++11 support
#include <functional>
#include <map>
#include <numeric>
#include <stdexcept>
#include <vector>
#include <sstream>

#include <eigen3/Eigen/Core>

#include <Python.h>

#include <gtest/gtest.h>

#include "utils.h"


#ifndef WITHOUT_NUMPY
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#endif // WITHOUT_NUMPY

#if PY_MAJOR_VERSION >= 3
#define PyString_FromString PyUnicode_FromString
#define PyInt_FromLong PyLong_FromLong
#define PyString_FromString PyUnicode_FromString
#endif

TEST ( Utils, ArrayExtraction )
{
    std::string test_string{"1.2 2 3 4 5 1e220"};
    std::vector<double> values;

    convert_array ( test_string, values, std::stod );
    for ( auto item : values )
    {
        std::cout << item << std::endl;
    }
}

TEST ( Utils, EigenMatrixToVector )
{
    unsigned int columns{5}, rows{5};
    Eigen::MatrixXd dataset ( columns, rows );
    for ( unsigned int i = 0; i < columns; ++i )
    {
        for ( unsigned int j = 0; j < rows; ++j )
        {
            dataset ( i, j ) = i + j * 2;
        }
    }

    auto item = dataset.col ( 0 );
    std::cout << item << std::endl;
}

TEST ( Utils, ChopArrayExtraction )
{
    std::string test_string{"1.2 2, 3 4, 5 1e220"};
    std::vector<std::vector<double>> values;

    chop_string_arrays ( test_string, ',', values, std::stod );
    for ( auto &item : values )
    {
        for ( auto &item2 : item )
        {
            std::cout << item2 << " ";
        }
        std::cout << std::endl;
    }
}


class colormap
{
protected:
  PyObject *m_python_colormap{nullptr};

public:
  colormap (PyObject *s_python_class_colormaps, std::string colormap_scheme = "plasma")
  {
    PyObject * m_python_colormap_name = PyString_FromString(colormap_scheme.c_str());
    m_python_colormap = PyObject_GetItem(s_python_class_colormaps, m_python_colormap_name);

    if (!m_python_colormap)
    {
      throw std::runtime_error("couldnt instantiate colormap");
    }
    Py_DECREF(m_python_colormap_name);
  }

  PyObject * operator()(double param)
  {
    PyObject *param_object = PyFloat_FromDouble(param);
    PyObject *args = PyTuple_Pack(1, param_object);
    PyObject *color = PyObject_CallObject(m_python_colormap, args);
    Py_DECREF(param_object);
    Py_DECREF(args);
    Py_DECREF(param_object);

    if (color)
    {
      return color;
    }
    else
    {
      throw std::runtime_error("couldnt create color");
    }
  }

  ~colormap ()
  {
    if (m_python_colormap)
    {
      Py_DECREF(m_python_colormap);
    }
  }
};

TEST ( Utils, PythonMatplotlib)
{
    PyObject *s_python_colormap;
    PyObject *s_python_mpl_color;
    PyObject *s_python_class_colormaps;

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
      //std::cout << Py_GetPythonHome() << std::endl;
    }

    PyObject *matplotlibname = PyString_FromString("matplotlib");
    PyObject *pyplotname = PyString_FromString("matplotlib.pyplot");
    PyObject *cmname = PyString_FromString("matplotlib.cm");
    PyObject *colors_name = PyString_FromString("matplotlib.colors");
    PyObject *pylabname = PyString_FromString("pylab");
    if (!pyplotname || !pylabname || !matplotlibname || !cmname || !colors_name) {
      throw std::runtime_error("couldnt create string");
    }

    PyObject *matplotlib = PyImport_Import(matplotlibname);
    Py_DECREF(matplotlibname);
    if (!matplotlib) {
      PyErr_Print();
      throw std::runtime_error("Error loading module matplotlib!");
    }

    std::string s_backend;

    // matplotlib.use() must be called *before* pylab, matplotlib.pyplot,
    // or matplotlib.backends is imported for the first time
    if (!s_backend.empty()) {
      PyObject_CallMethod(matplotlib, const_cast<char *>("use"),
                          const_cast<char *>("s"), s_backend.c_str());
    }

    PyObject *pymod = PyImport_Import(pyplotname);
    Py_DECREF(pyplotname);
    if (!pymod) {
      throw std::runtime_error("Error loading module matplotlib.pyplot!");
    }

    s_python_colormap = PyImport_Import(cmname);
    Py_DECREF(cmname);
    if (!s_python_colormap) {
      throw std::runtime_error("Error loading module matplotlib.cm!");
    }

    s_python_mpl_color = PyImport_Import(colors_name);
    Py_DECREF(colors_name);
    if (!s_python_mpl_color) {
      throw std::runtime_error("Error loading module matplotlib.color!");
    }

    PyObject *pylabmod = PyImport_Import(pylabname);
    Py_DECREF(pylabname);
    if (!pylabmod) {
      throw std::runtime_error("Error loading module pylab!");
    }

    s_python_class_colormaps = PyObject_GetAttrString(matplotlib, "colormaps");
    //if (!PyDict_Check(s_python_class_colormaps))
    //{
    //  throw std::runtime_error("Python object is unexpectedly not a PyFunction.");
    //}

    colormap map_colors(s_python_class_colormaps);
    auto colore = map_colors(0.1);
}
