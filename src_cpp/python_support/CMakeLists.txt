
add_library(matplotlib_wrapper_STATIC STATIC matplotlib_wrapper.cpp)
add_library(matplotlib_wrapper_SHARED SHARED matplotlib_wrapper.cpp)
set_target_properties(matplotlib_wrapper_SHARED PROPERTIES VERSION ${random_samples_VERSION} SOVERSION ${random_samples_VERSION_MAJOR} OUTPUT_NAME matplotlib_wrapper)
set_target_properties(matplotlib_wrapper_STATIC PROPERTIES OUTPUT_NAME matplotlib_wrapper)
target_include_directories(matplotlib_wrapper_SHARED PUBLIC ${CMAKE_SOURCE_DIR}/third-party/pybind11/include)
target_include_directories(matplotlib_wrapper_STATIC PUBLIC ${CMAKE_SOURCE_DIR}/third-party/pybind11/include)

target_link_libraries(matplotlib_wrapper_STATIC PUBLIC
        ${Boost_LIBRARIES}
)
target_link_libraries(matplotlib_wrapper_SHARED PUBLIC
        ${Boost_LIBRARIES}
)
