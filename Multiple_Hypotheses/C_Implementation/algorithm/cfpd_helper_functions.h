#ifndef CFPD_HELPER_FUNCTIONS_H
#define CFPD_HELPER_FUNCTIONS_H


double
cfpd_vector_sum(const gsl_vector *v);


int
cfpd_vector_finite(const gsl_vector *v);


int
cfpd_matrix_geq(const gsl_matrix *m1,
                const gsl_matrix *m2);


int
cfpd_matrix_finite(const gsl_matrix *m);


#endif  //CFPD_HELPER_FUNCTIONS_H
