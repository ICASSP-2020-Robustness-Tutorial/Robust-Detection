#ifndef CFPD_TESTS_OBJECTIVE_FUNCTION_H
#define CFPD_TESTS_OBJECTIVE_FUNCTION_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

// Objective function
double
weighted_kl(const gsl_vector *x,
            const size_t      k,
            const void       *params);


// Objective function derivatives
double
weighted_kl_derivative(const size_t      n,
                       const gsl_vector *x,
                       const size_t      k,
                       const void       *params);


#endif  // CFPD_TESTS_OBJECTIVE_FUNCTION_H
