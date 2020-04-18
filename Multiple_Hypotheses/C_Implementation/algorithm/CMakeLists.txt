find_library(MATH_LIB m)
find_library(GSL_LIB gsl)
find_library(BLAS_LIB openblas)

add_library(lfds SHARED lfds_algorithm.c
                        lfds_opt_problem.c
                        lfds_opt_problem_checks.c
                        lfds_errors.c
                        lfds_helper_functions.c)

target_compile_options(lfds PUBLIC -Wall -pedantic -O2 -march=native -DHAVE_INLINE=1 -DGSL_RANGE_CHECK=0)

target_link_libraries(lfds ${MATH_LIB} ${GSL_LIB} ${BLAS_LIB})

set_target_properties(lfds PROPERTIES PUBLIC_HEADER lfds.h)

target_include_directories(lfds PUBLIC "/usr/local/include")

install(TARGETS lfds
    LIBRARY DESTINATION "/usr/local/lib"
    PUBLIC_HEADER DESTINATION "/usr/local/include")
