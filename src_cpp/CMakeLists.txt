
file(GLOB SOURCES *.cpp)

project(random_samples VERSION 1.0.0)

add_library(random_samples_STATIC STATIC random_samples.cpp)
add_library(random_samples_SHARED SHARED random_samples.cpp)
set_target_properties(random_samples_SHARED PROPERTIES
        VERSION ${random_samples_VERSION}
        SOVERSION ${random_samples_VERSION_MAJOR} OUTPUT_NAME random_samples)
set_target_properties(random_samples_STATIC PROPERTIES
        OUTPUT_NAME random_samples)

target_include_directories(random_samples_STATIC PUBLIC ../third-party/libstable/stable/src)
target_include_directories(random_samples_SHARED PUBLIC ../third-party/libstable/stable/src)
target_link_libraries(random_samples_STATIC PUBLIC
        ${Boost_LIBRARIES} KDTree stable_STATIC
)
target_link_libraries(random_samples_SHARED PUBLIC
        ${Boost_LIBRARIES} KDTree stable_SHARED
)

add_executable(transfer_entropy main.cpp renyi_entropy.cpp)
target_include_directories(transfer_entropy PUBLIC ../third-party/ndarray/include ../third-party/KDTree ../third-party/cpp-btree)
target_link_libraries(transfer_entropy PUBLIC
        ${Boost_LIBRARIES} KDTree stable_SHARED
)

add_executable(test_renyi_entropy test_renyi_entropy.cpp renyi_entropy.cpp)
target_include_directories(test_renyi_entropy PUBLIC ../third-party/ndarray/include ../third-party/KDTree ../third-party/cpp-btree)
target_link_libraries(test_renyi_entropy PUBLIC
        ${Boost_LIBRARIES} KDTree GTest::gtest GTest::gtest_main Eigen3::Eigen
)

add_executable(test_random_samples test_random_samples.cpp)
target_include_directories(test_random_samples PUBLIC ../third-party/ndarray/include ../third-party/KDTree ../third-party/cpp-btree)
target_link_libraries(test_random_samples PUBLIC
        ${Boost_LIBRARIES} KDTree GTest::gtest GTest::gtest_main Eigen3::Eigen random_samples_SHARED
)

add_executable(test_renyi_entropy_gaussian_distribution gaussian_distribution_mean.cpp renyi_entropy.cpp)
target_include_directories(test_renyi_entropy_gaussian_distribution PUBLIC ../third-party/ndarray/include ../third-party/KDTree ../third-party/cpp-btree)
target_link_libraries(test_renyi_entropy_gaussian_distribution PUBLIC
        ${Boost_LIBRARIES} KDTree GTest::gtest GTest::gtest_main Eigen3::Eigen random_samples_STATIC
)

add_executable(test_renyi_entropy_multi_gaussian_distribution gaussian_distribution_multi_mean.cpp renyi_entropy.cpp)
target_include_directories(test_renyi_entropy_multi_gaussian_distribution PUBLIC ../third-party/ndarray/include ../third-party/KDTree ../third-party/cpp-btree)
target_link_libraries(test_renyi_entropy_multi_gaussian_distribution PUBLIC
        ${Boost_LIBRARIES} KDTree GTest::gtest GTest::gtest_main ${GMP_LIBRARY} ${MPFR_LIBRARY} Eigen3::Eigen random_samples_SHARED
)

add_executable(binary_file_test read_binary_dataset_mean.cpp)
target_include_directories(binary_file_test PUBLIC ../third-party/ndarray/include ../third-party/KDTree ../third-party/cpp-btree)
target_link_libraries(binary_file_test PUBLIC
        ${Boost_LIBRARIES} KDTree GTest::gtest GTest::gtest_main Eigen3::Eigen
)

add_executable(renyi_entropy_calculator renyi_entropy_calculator_main.cpp renyi_entropy.cpp utils.cpp cpptrace_helper.cpp)
target_include_directories(renyi_entropy_calculator PUBLIC ../third-party/ndarray/include ../third-party/KDTree ../third-party/cpp-btree ../third-party/msgpack-c/include ../third-party/cpptrace/include)
target_link_libraries(renyi_entropy_calculator PUBLIC 
        ${Boost_LIBRARIES} KDTree stable_SHARED random_samples_SHARED cpptrace::cpptrace
)

add_executable(test_utils test_utils.cpp utils.cpp)
target_include_directories(test_utils PUBLIC ../third-party/ndarray/include ../third-party/KDTree ../third-party/cpp-btree)
target_link_libraries(test_utils PUBLIC
        ${Boost_LIBRARIES} KDTree GTest::gtest GTest::gtest_main Eigen3::Eigen random_samples_SHARED
)

add_executable(test_samples_store test_samples_store.cpp)
target_include_directories(test_samples_store PUBLIC ../third-party/ndarray/include ../third-party/KDTree ../third-party/cpp-btree ../third-party/msgpack-c/include)
target_link_libraries(test_samples_store PUBLIC
        ${Boost_LIBRARIES} KDTree GTest::gtest GTest::gtest_main Eigen3::Eigen random_samples_SHARED
)

add_executable(signal_tracer signal_tracer.cpp)
target_include_directories(signal_tracer PUBLIC ../third-party/cpptrace/include)
target_link_libraries(signal_tracer PUBLIC cpptrace::cpptrace)

