
add_library(random_samples_STATIC STATIC random_samples.cpp)
add_library(random_samples_SHARED SHARED random_samples.cpp)
set_target_properties(random_samples_SHARED PROPERTIES
        VERSION ${random_samples_VERSION}
        SOVERSION ${random_samples_VERSION_MAJOR} OUTPUT_NAME random_samples)
set_target_properties(random_samples_STATIC PROPERTIES
        OUTPUT_NAME random_samples)

target_include_directories(random_samples_STATIC PUBLIC ${CMAKE_SOURCE_DIR}/third-party/libstable/stable/src)
target_include_directories(random_samples_SHARED PUBLIC ${CMAKE_SOURCE_DIR}/third-party/libstable/stable/src)
target_link_libraries(random_samples_STATIC PUBLIC
        ${Boost_LIBRARIES} KDTree stable_STATIC ${MPFR_LIBRARY}
)
target_link_libraries(random_samples_SHARED PUBLIC
        ${Boost_LIBRARIES} KDTree stable_SHARED ${MPFR_LIBRARY}
)

add_executable(test_random_samples test_random_samples.cpp)
target_include_directories(test_random_samples PUBLIC ${PROJECT_SOURCE_DIR}/third-party/ndarray/include ${PROJECT_SOURCE_DIR}/third-party/KDTree ${PROJECT_SOURCE_DIR}/third-party/cpp-btree)
target_link_libraries(test_random_samples PUBLIC
        ${Boost_LIBRARIES} KDTree GTest::gtest GTest::gtest_main Eigen3::Eigen ${MPFR_LIBRARY} random_samples_SHARED
)

add_subdirectory(pyrandom_samples)
