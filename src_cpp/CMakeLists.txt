
file(GLOB SOURCES *.cpp)

project(random_samples VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 23)

add_subdirectory(prototype)
add_subdirectory(python_support)
add_subdirectory(random_samples)
add_subdirectory(renyi_entropy)

add_executable(test_renyi_entropy_gaussian_distribution gaussian_distribution_mean.cpp)
target_include_directories(test_renyi_entropy_gaussian_distribution PUBLIC ${PROJECT_SOURCE_DIR}/third-party/ndarray/include ${PROJECT_SOURCE_DIR}/third-party/KDTree ${PROJECT_SOURCE_DIR}/third-party/cpp-btree renyi_entropy random_samples)
target_link_libraries(test_renyi_entropy_gaussian_distribution PUBLIC
        ${Boost_LIBRARIES} KDTree GTest::gtest GTest::gtest_main Eigen3::Eigen random_samples_STATIC renyi_entropy_STATIC
)

add_executable(test_renyi_entropy_multi_gaussian_distribution gaussian_distribution_multi_mean.cpp)
target_include_directories(test_renyi_entropy_multi_gaussian_distribution PUBLIC ${PROJECT_SOURCE_DIR}/third-party/ndarray/include ${PROJECT_SOURCE_DIR}/third-party/KDTree ${PROJECT_SOURCE_DIR}/third-party/cpp-btree renyi_entropy random_samples)
target_link_libraries(test_renyi_entropy_multi_gaussian_distribution PUBLIC
        ${Boost_LIBRARIES} KDTree GTest::gtest GTest::gtest_main ${GMP_LIBRARY} ${MPFR_LIBRARY} Eigen3::Eigen random_samples_STATIC renyi_entropy_STATIC
)

add_executable(binary_file_test read_binary_dataset_mean.cpp)
target_include_directories(binary_file_test PUBLIC ${CMAKE_SOURCE_DIR}/third-party/ndarray/include ${CMAKE_SOURCE_DIR}/third-party/KDTree ${CMAKE_SOURCE_DIR}/third-party/cpp-btree ${CMAKE_SOURCE_DIR}/third-party/msgpack-c/include renyi_entropy)
target_link_libraries(binary_file_test PUBLIC
        ${Boost_LIBRARIES} KDTree GTest::gtest GTest::gtest_main Eigen3::Eigen
)

add_executable(renyi_entropy_calculator renyi_entropy_calculator_main.cpp utils.cpp cpptrace_helper.cpp)
target_include_directories(renyi_entropy_calculator PUBLIC ${CMAKE_SOURCE_DIR}/third-party/ndarray/include ${CMAKE_SOURCE_DIR}/third-party/KDTree ${CMAKE_SOURCE_DIR}/third-party/cpp-btree ${CMAKE_SOURCE_DIR}/third-party/msgpack-c/include ${CMAKE_SOURCE_DIR}/third-party/cpptrace/include ${CMAKE_SOURCE_DIR}/third-party/pybind11/include renyi_entropy random_samples python_support)
target_link_libraries(renyi_entropy_calculator PUBLIC 
        ${Boost_LIBRARIES} KDTree stable_SHARED random_samples_STATIC renyi_entropy_STATIC python_wrapper_STATIC cpptrace::cpptrace
)

add_executable(test_utils test_utils.cpp utils.cpp)
target_include_directories(test_utils PUBLIC ${CMAKE_SOURCE_DIR}/third-party/ndarray/include ${CMAKE_SOURCE_DIR}/third-party/KDTree ${CMAKE_SOURCE_DIR}/third-party/cpp-btree)
target_link_libraries(test_utils PUBLIC
        ${Boost_LIBRARIES} KDTree GTest::gtest GTest::gtest_main Eigen3::Eigen random_samples_SHARED
)

add_executable(test_samples_store test_samples_store.cpp)
target_include_directories(test_samples_store PUBLIC ${CMAKE_SOURCE_DIR}/third-party/ndarray/include ${CMAKE_SOURCE_DIR}/third-party/KDTree ${CMAKE_SOURCE_DIR}/third-party/cpp-btree ${CMAKE_SOURCE_DIR}/third-party/msgpack-c/include renyi_entropy)
target_link_libraries(test_samples_store PUBLIC
        ${Boost_LIBRARIES} KDTree GTest::gtest GTest::gtest_main Eigen3::Eigen cpptrace::cpptrace random_samples_SHARED
)

add_executable(renyi_entropy_figures renyi_entropy_figures_main.cpp cpptrace_helper.cpp)
target_include_directories(renyi_entropy_figures PUBLIC ${CMAKE_SOURCE_DIR}/third-party/cpptrace/include ${CMAKE_SOURCE_DIR}/third-party/KDTree ${CMAKE_SOURCE_DIR}/third-party/msgpack-c/include ${CMAKE_SOURCE_DIR}/third-party/DataFrame/include ${CMAKE_SOURCE_DIR}/third-party/fmt/include ${CMAKE_SOURCE_DIR}/third-party/matplotlib-cpp ${Python_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/third-party/pybind11/include python_support renyi_entropy)
target_link_libraries(renyi_entropy_figures PUBLIC python_wrapper_STATIC ${Boost_LIBRARIES} GTest::gtest GTest::gtest_main Eigen3::Eigen cpptrace::cpptrace )
