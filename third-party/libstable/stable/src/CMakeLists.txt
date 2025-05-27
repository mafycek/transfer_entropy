
set(STABLE_SOURCES mcculloch.c  stable_cdf.c     stable_dist.c  stable_integration.c   stable_pdf.c  stable_rnd.c
methods.c stable_common.c  stable_fit.c   stable_koutrouvelis.c  stable_q.c)

add_library(stable SHARED ${STABLE_SOURCES})
target_link_libraries(stable GSL::gsl GSL::gslcblas -lm -lpthread)
set_target_properties(stable PROPERTIES LINKER_LANGUAGE CXX)
