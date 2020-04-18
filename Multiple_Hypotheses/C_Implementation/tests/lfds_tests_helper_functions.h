#ifndef CFPD_TESTS_HELPER_FUNCTIONS_H
#define CFPD_TESTS_HELPER_FUNCTIONS_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

double
lfds_fload(const char* filename);


int
lfds_vector_load(const char *filename,
                 gsl_vector *v);


int
lfds_vector_approx_equal(const gsl_vector *v1,
                         const gsl_vector *v2,
                         const double      eps);


int
lfds_matrix_load(const char *filename,
                 gsl_matrix *m);


int
lfds_matrix_approx_equal(const gsl_matrix *m1,
                         const gsl_matrix *m2,
                         const double      eps);

#endif  // CFPD_TESTS_HELPER_FUNCTIONS_H
