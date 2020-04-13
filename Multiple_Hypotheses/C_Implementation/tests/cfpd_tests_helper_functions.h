#ifndef CFPD_TESTS_HELPER_FUNCTIONS_H
#define CFPD_TESTS_HELPER_FUNCTIONS_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

double
cfpd_fload(const char* filename);


int
cfpd_vector_load(const char *filename,
                 gsl_vector *v);


int
cfpd_vector_approx_equal(const gsl_vector *v1,
                         const gsl_vector *v2,
                         const double      eps);


int
cfpd_matrix_load(const char *filename,
                 gsl_matrix *m);


int
cfpd_matrix_approx_equal(const gsl_matrix *m1,
                         const gsl_matrix *m2,
                         const double      eps);

#endif  // CFPD_TESTS_HELPER_FUNCTIONS_H
