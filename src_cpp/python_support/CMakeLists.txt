
set(SOURCES generic_python_wrapper.cpp matplotlib_wrapper.cpp pandas_wrapper.cpp)

add_library(python_wrapper_STATIC STATIC ${SOURCES})
add_library(python_wrapper_SHARED SHARED ${SOURCES})
set_target_properties(python_wrapper_SHARED PROPERTIES VERSION ${random_samples_VERSION} SOVERSION ${random_samples_VERSION_MAJOR} OUTPUT_NAME matplotlib_wrapper)
set_target_properties(python_wrapper_STATIC PROPERTIES OUTPUT_NAME matplotlib_wrapper)
target_include_directories(python_wrapper_SHARED PUBLIC ${CMAKE_SOURCE_DIR}/third-party/pybind11/include)
target_include_directories(python_wrapper_STATIC PUBLIC ${CMAKE_SOURCE_DIR}/third-party/pybind11/include)

target_link_libraries(python_wrapper_STATIC PUBLIC
        ${Boost_LIBRARIES}
)
target_link_libraries(python_wrapper_SHARED PUBLIC
        ${Boost_LIBRARIES}
)
