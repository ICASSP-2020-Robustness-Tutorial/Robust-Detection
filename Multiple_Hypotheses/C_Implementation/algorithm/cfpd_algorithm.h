#ifndef CFPD_ALGORITHM_H
#define CFPD_ALGORITHM_H

#include "cfpd_function.h"
#include "cfpd_function_vector.h"
#include "cfpd_opt_problem.h"


double
cfpd_residual_df(double  x,
                 void   *params);


double
cfpd_update_residuals_objective(cfpd_opt_problem_t *opt_problem);


double
cfpd_update_residuals_densities(cfpd_opt_problem_t *opt_problem);


void
cfpd_optimize_density(cfpd_opt_problem_t *opt_problem);


double
cfpd_residual_density(double  c,
                      void   *params);


void
cfpd_update_density(cfpd_opt_problem_t *opt_problem,
                    const size_t        n);


int
cfpd_minimize(cfpd_opt_problem_t *opt_problem);


int
cfpd_minimize_proximal(cfpd_opt_problem_t *opt_problem);


#endif  // CFPD_ALGORITHM_H
