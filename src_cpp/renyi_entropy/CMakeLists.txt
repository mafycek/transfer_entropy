
add_library(renyi_entropy_STATIC STATIC renyi_entropy.cpp)
add_library(renyi_entropy_SHARED SHARED renyi_entropy.cpp)
set_target_properties(renyi_entropy_SHARED PROPERTIES VERSION ${random_samples_VERSION} SOVERSION ${random_samples_VERSION_MAJOR} OUTPUT_NAME matplotlib_wrapper)
set_target_properties(renyi_entropy_STATIC PROPERTIES OUTPUT_NAME matplotlib_wrapper)
target_include_directories(renyi_entropy_SHARED PUBLIC ${CMAKE_SOURCE_DIR}/third-party/KDTree)
target_include_directories(renyi_entropy_STATIC PUBLIC ${CMAKE_SOURCE_DIR}/third-party/KDTree)

target_link_libraries(renyi_entropy_STATIC PUBLIC
        ${Boost_LIBRARIES}
)
target_link_libraries(renyi_entropy_SHARED PUBLIC
        ${Boost_LIBRARIES}
)

add_executable(transfer_entropy test_RTE.cpp)
target_include_directories(transfer_entropy PUBLIC ${CMAKE_SOURCE_DIR}/third-party/ndarray/include ${CMAKE_SOURCE_DIR}/third-party/KDTree ${CMAKE_SOURCE_DIR}/third-party/cpp-btree)
target_link_libraries(transfer_entropy PUBLIC
        ${Boost_LIBRARIES} KDTree stable_SHARED renyi_entropy_STATIC
)

add_executable(test_renyi_entropy test_renyi_entropy.cpp)
target_include_directories(test_renyi_entropy PUBLIC ../third-party/ndarray/include ../third-party/KDTree ../third-party/cpp-btree)
target_link_libraries(test_renyi_entropy PUBLIC
        ${Boost_LIBRARIES} KDTree GTest::gtest GTest::gtest_main Eigen3::Eigen renyi_entropy_STATIC
)
