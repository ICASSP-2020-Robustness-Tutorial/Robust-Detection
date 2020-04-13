find_library(MATH_LIB m)
find_library(GSL_LIB gsl)
find_library(BLAS_LIB openblas)

add_executable(cfpd_tests cfpd_tests.c
			  cfpd_tests_objective_function.c
                          cfpd_tests_helper_functions.c)

target_compile_options(cfpd_tests PUBLIC -Wall -pedantic -O2 -march=native)

target_link_libraries(cfpd_tests LINK_PUBLIC ${MATH_LIB} ${GSL_LIB} ${BLAS_LIB} cfpd)
