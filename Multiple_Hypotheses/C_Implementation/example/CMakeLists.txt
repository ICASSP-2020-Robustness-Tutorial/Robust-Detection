find_library(MATH_LIB m)
find_library(GSL_LIB gsl)
find_library(BLAS_LIB openblas)

add_executable(cfpd_example cfpd_example.c)

target_compile_options(cfpd_example PUBLIC -Wall -pedantic -O2 -march=native)

target_link_libraries(cfpd_example LINK_PUBLIC ${MATH_LIB} ${GSL_LIB} ${BLAS_LIB} cfpd)
