
add_executable(signal_tracer signal_tracer.cpp)
target_include_directories(signal_tracer PUBLIC ${CMAKE_SOURCE_DIR}/third-party/cpptrace/include)
target_link_libraries(signal_tracer PUBLIC cpptrace::cpptrace)

add_executable(test_fftw3 test_fftw3.cpp)
target_include_directories(test_fftw3 PUBLIC ${CMAKE_SOURCE_DIR}/third-party/cpptrace/include)
target_link_libraries(test_fftw3 PUBLIC cpptrace::cpptrace fft 2Dfft GTest::gtest GTest::gtest_main)

add_executable(test_dataframe test_dataframe.cpp)
target_include_directories(test_dataframe PUBLIC ${CMAKE_SOURCE_DIR}/third-party/cpptrace/include ${CMAKE_SOURCE_DIR}/third-party/msgpack-c/include ${CMAKE_SOURCE_DIR}/third-party/DataFrame/include ${CMAKE_SOURCE_DIR}/third-party/fmt/include ${CMAKE_SOURCE_DIR}/third-party/matplotlib-cpp ${Python_INCLUDE_DIRS})
target_link_libraries(test_dataframe PUBLIC ${Boost_LIBRARIES} GTest::gtest GTest::gtest_main Eigen3::Eigen cpptrace::cpptrace )

add_executable(test_python test_python.cpp)
target_include_directories(test_python PUBLIC ${CMAKE_SOURCE_DIR}/third-party/pybind11/include)
target_link_libraries(test_python PUBLIC
        ${Boost_LIBRARIES} GTest::gtest GTest::gtest_main
)
