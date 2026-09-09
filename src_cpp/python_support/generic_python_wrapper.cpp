
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
#include "generic_python_wrapper.h"
#pragma GCC diagnostic pop

namespace python_wrappers
{

bool generic_python_wrapper::_initialized = false;

generic_python_wrapper::generic_python_wrapper()
{
  if (!_initialized)
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
  }
}

generic_python_wrapper::~generic_python_wrapper()
{
}

}
