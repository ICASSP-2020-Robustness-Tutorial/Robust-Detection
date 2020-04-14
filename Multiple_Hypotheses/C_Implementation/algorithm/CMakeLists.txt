find_library(MATH_LIB m)
find_library(GSL_LIB gsl)
find_library(BLAS_LIB openblas)

add_library(cfpd SHARED cfpd_algorithm.c
                        cfpd_opt_problem.c
                        cfpd_opt_problem_checks.c
                        cfpd_errors.c
                        cfpd_helper_functions.c)

target_compile_options(cfpd PUBLIC -Wall -pedantic -O2 -march=native -DHAVE_INLINE=1 -DGSL_RANGE_CHECK=0)

target_link_libraries(cfpd ${MATH_LIB} ${GSL_LIB} ${BLAS_LIB})

set_target_properties(cfpd PROPERTIES PUBLIC_HEADER cfpd.h)

target_include_directories(cfpd PUBLIC "/home/mfauss/Projects/Cfpd/include")

install(TARGETS cfpd
    LIBRARY DESTINATION "/home/mfauss/Projects/Cfpd/lib"
    PUBLIC_HEADER DESTINATION "/home/mfauss/Projects/Cfpd/include")