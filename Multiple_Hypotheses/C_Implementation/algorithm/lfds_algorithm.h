#ifndef CFPD_ALGORITHM_H
#define CFPD_ALGORITHM_H

#include "lfds_function.h"
#include "lfds_function_vector.h"
#include "lfds_opt_problem.h"


double
lfds_residual_df(double  x,
                 void   *params);


double
lfds_update_residuals_objective(lfds_opt_problem_t *opt_problem);


double
lfds_update_residuals_densities(lfds_opt_problem_t *opt_problem);


void
lfds_optimize_density(lfds_opt_problem_t *opt_problem);


double
lfds_residual_density(double  c,
                      void   *params);


void
lfds_update_density(lfds_opt_problem_t *opt_problem,
                    const size_t        n);


int
lfds_minimize(lfds_opt_problem_t *opt_problem);


int
lfds_minimize_proximal(lfds_opt_problem_t *opt_problem);


#endif  // CFPD_ALGORITHM_H
