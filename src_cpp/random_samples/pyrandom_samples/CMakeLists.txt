
pybind11_add_module(pyrandom_samples MODULE pyrandom_samples.cpp)

target_include_directories(pyrandom_samples PUBLIC "..")

# Place the compiled extension (.so) in that folder
set_target_properties(pyrandom_samples
  PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "."
    SUFFIX ".so"
)

target_link_libraries(pyrandom_samples PUBLIC
        ${Boost_LIBRARIES} KDTree GTest::gtest GTest::gtest_main Eigen3::Eigen ${MPFR_LIBRARY} random_samples_SHARED
)



