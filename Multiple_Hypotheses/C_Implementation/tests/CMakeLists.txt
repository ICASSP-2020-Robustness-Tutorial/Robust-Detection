find_library(MATH_LIB m)
find_library(GSL_LIB gsl)
find_library(BLAS_LIB openblas)

add_executable(lfds_tests lfds_tests.c
			  lfds_tests_objective_function.c
                          lfds_tests_helper_functions.c)

target_compile_options(lfds_tests PUBLIC -Wall -pedantic -O2 -march=native)

target_link_libraries(lfds_tests LINK_PUBLIC ${MATH_LIB} ${GSL_LIB} ${BLAS_LIB} lfds)
