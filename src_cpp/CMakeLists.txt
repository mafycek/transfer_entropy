
file(GLOB SOURCES *.cpp)

add_library(random_samples SHARED random_samples.cpp)
target_include_directories(random_samples PUBLIC ../third-party/libstable/stable/src)
target_link_libraries(random_samples PUBLIC
        ${Boost_LIBRARIES} KDTree stable
)

add_executable(transfer_entropy main.cpp renyi_entropy.cpp)
target_include_directories(transfer_entropy PUBLIC ../third-party/ndarray/include ../third-party/KDTree ../third-party/cpp-btree)
target_link_libraries(transfer_entropy PUBLIC
        ${Boost_LIBRARIES} KDTree stable
)

add_executable(test_renyi_entropy test_renyi_entropy.cpp renyi_entropy.cpp)
target_include_directories(test_renyi_entropy PUBLIC ../third-party/ndarray/include ../third-party/KDTree ../third-party/cpp-btree)
target_link_libraries(test_renyi_entropy PUBLIC
        ${Boost_LIBRARIES} KDTree GTest::gtest GTest::gtest_main Eigen3::Eigen
)

add_executable(test_random_samples test_random_samples.cpp)
target_include_directories(test_random_samples PUBLIC ../third-party/ndarray/include ../third-party/KDTree ../third-party/cpp-btree)
target_link_libraries(test_random_samples PUBLIC
        ${Boost_LIBRARIES} KDTree GTest::gtest GTest::gtest_main Eigen3::Eigen random_samples 
)

add_executable(test_renyi_entropy_gaussian_distribution gaussian_distribution_mean.cpp renyi_entropy.cpp)
target_include_directories(test_renyi_entropy_gaussian_distribution PUBLIC ../third-party/ndarray/include ../third-party/KDTree ../third-party/cpp-btree)
target_link_libraries(test_renyi_entropy_gaussian_distribution PUBLIC
        ${Boost_LIBRARIES} KDTree GTest::gtest GTest::gtest_main Eigen3::Eigen random_samples
)

add_executable(test_renyi_entropy_multi_gaussian_distribution gaussian_distribution_multi_mean.cpp renyi_entropy.cpp)
target_include_directories(test_renyi_entropy_multi_gaussian_distribution PUBLIC ../third-party/ndarray/include ../third-party/KDTree ../third-party/cpp-btree)
target_link_libraries(test_renyi_entropy_multi_gaussian_distribution PUBLIC
        ${Boost_LIBRARIES} KDTree GTest::gtest GTest::gtest_main ${GMP_LIBRARY} ${MPFR_LIBRARY} Eigen3::Eigen random_samples 
)

add_executable(binary_file_test read_binary_dataset_mean.cpp)
target_include_directories(binary_file_test PUBLIC ../third-party/ndarray/include ../third-party/KDTree ../third-party/cpp-btree)
target_link_libraries(binary_file_test PUBLIC
        ${Boost_LIBRARIES} KDTree GTest::gtest GTest::gtest_main Eigen3::Eigen
)

add_executable(renyi_entropy_calculator renyi_entropy_calculator_main.cpp renyi_entropy.cpp utils.cpp)
target_include_directories(renyi_entropy_calculator PUBLIC ../third-party/ndarray/include ../third-party/KDTree ../third-party/cpp-btree)
target_link_libraries(renyi_entropy_calculator PUBLIC 
        ${Boost_LIBRARIES} KDTree stable random_samples
)

add_executable(test_utils test_utils.cpp utils.cpp)
target_include_directories(test_utils PUBLIC ../third-party/ndarray/include ../third-party/KDTree ../third-party/cpp-btree)
target_link_libraries(test_utils PUBLIC
        ${Boost_LIBRARIES} KDTree GTest::gtest GTest::gtest_main Eigen3::Eigen random_samples 
)
