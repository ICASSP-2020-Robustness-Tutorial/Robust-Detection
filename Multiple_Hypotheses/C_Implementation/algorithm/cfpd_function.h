#ifndef CFPD_FUNCTION_H
#define CFPD_FUNCTION_H

#include <gsl/gsl_vector.h>


typedef double(*cfpd_function_t)(const gsl_vector *x,
                                 const size_t      k,
                                 const void       *params);


typedef double(*cfpd_derivative_t)(const size_t      n,
                                   const gsl_vector *x,
                                   const size_t      k,
                                   const void       *params);

#endif  //CFPD_FUNCTION_H
