
set(STABLE_SOURCES mcculloch.c  stable_cdf.c     stable_dist.c  stable_integration.c   stable_pdf.c  stable_rnd.c
methods.c stable_common.c  stable_fit.c   stable_koutrouvelis.c  stable_q.c)

add_library(stable_STATIC STATIC ${STABLE_SOURCES})
add_library(stable_SHARED SHARED ${STABLE_SOURCES})
set_target_properties(stable_SHARED PROPERTIES
                      VERSION ${libstable_VERSION}
                      SOVERSION ${libstable_VERSION_MAJOR} OUTPUT_NAME stable)
set_target_properties(stable_STATIC PROPERTIES OUTPUT_NAME stable)
target_link_libraries(stable_STATIC GSL::gsl GSL::gslcblas -lm -lpthread)
target_link_libraries(stable_SHARED GSL::gsl GSL::gslcblas -lm -lpthread)
set_target_properties(stable_STATIC PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(stable_SHARED PROPERTIES LINKER_LANGUAGE CXX)
