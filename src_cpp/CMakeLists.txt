
file(GLOB SOURCES *.cpp)

add_executable(transfer_entropy main.cpp renyi_entropy.cpp)

target_link_libraries(transfer_entropy PUBLIC
        ${Boost_LIBRARIES} KDTree
)

target_include_directories(transfer_entropy PUBLIC ../third-party/ndarray/include ../third-party/KDTree ../third-party/cpp-btree)

add_executable(test_renyi_entropy test_renyi_entropy.cpp renyi_entropy.cpp)
target_link_libraries(test_renyi_entropy PUBLIC
        ${Boost_LIBRARIES} KDTree GTest::gtest GTest::gtest_main Eigen3::Eigen
)
target_include_directories(test_renyi_entropy PUBLIC ../third-party/ndarray/include ../third-party/KDTree ../third-party/cpp-btree)

add_executable(test_random_samples test_random_samples.cpp random_samples.cpp)
target_link_libraries(test_random_samples PUBLIC
        ${Boost_LIBRARIES} KDTree GTest::gtest GTest::gtest_main Eigen3::Eigen
)
target_include_directories(test_random_samples PUBLIC ../third-party/ndarray/include ../third-party/KDTree ../third-party/cpp-btree)

add_executable(test_renyi_entropy_gaussian_distribution gaussian_distribution_mean.cpp renyi_entropy.cpp random_samples.cpp)
target_link_libraries(test_renyi_entropy_gaussian_distribution PUBLIC
        ${Boost_LIBRARIES} KDTree GTest::gtest GTest::gtest_main Eigen3::Eigen
)
target_include_directories(test_renyi_entropy_gaussian_distribution PUBLIC ../third-party/ndarray/include ../third-party/KDTree ../third-party/cpp-btree)

add_executable(test_renyi_entropy_multi_gaussian_distribution gaussian_distribution_multi_mean.cpp renyi_entropy.cpp random_samples.cpp)
target_link_libraries(test_renyi_entropy_multi_gaussian_distribution PUBLIC
        ${Boost_LIBRARIES} KDTree GTest::gtest GTest::gtest_main ${GMP_LIBRARY} ${MPFR_LIBRARY} Eigen3::Eigen
)
target_include_directories(test_renyi_entropy_multi_gaussian_distribution PUBLIC ../third-party/ndarray/include ../third-party/KDTree ../third-party/cpp-btree)

add_executable(binary_file_test read_binary_dataset_mean.cpp)
target_link_libraries(binary_file_test PUBLIC
        ${Boost_LIBRARIES} KDTree GTest::gtest GTest::gtest_main Eigen3::Eigen
)
target_include_directories(binary_file_test PUBLIC ../third-party/ndarray/include ../third-party/KDTree ../third-party/cpp-btree)
