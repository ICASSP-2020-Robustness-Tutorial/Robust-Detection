#ifndef CFPD_OPT_PROBLEM_H
#define CFPD_OPT_PROBLEM_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_roots.h>

#include "lfds_function.h"

typedef struct {
    lfds_function_t f;
    lfds_derivative_t df;
    void *f_params;
    gsl_matrix *P_min;
    gsl_matrix *P_max;
    gsl_matrix *P;
    gsl_matrix *P_proximal;
    gsl_vector *p_min_mass;
    gsl_vector *p_max_mass;
    gsl_vector *c_min;
    gsl_vector *c_max;
    gsl_vector *c;
    gsl_vector *res_objective;
    gsl_vector *res_densities;
    double mu;
    double eps_p;
    double eps_c;
    double eps_objective;
    double eps_densities;
    size_t n;
    size_t N;
    size_t k;
    size_t K;
    size_t itmax;
    size_t itmax_proximal;
    size_t iter;
    size_t iter_proximal;
    int user_defined_P;
    int user_defined_c;
    int proximal;
    int verbosity;
    int status;
    const gsl_root_fsolver_type *root_solver;
} lfds_opt_problem_t;


// Initialize struct


lfds_opt_problem_t*
lfds_opt_problem_new(const size_t N,
                     const size_t K,
                     const double mu);


void
lfds_opt_problem_free(lfds_opt_problem_t *opt_problem);


// Functions


void
lfds_opt_problem_reset(lfds_opt_problem_t *opt_problem);


void
lfds_opt_problem_update_P_proximal(lfds_opt_problem_t *opt_problem);


void
lfds_opt_problem_set_c_bounds_auto(lfds_opt_problem_t  *opt_problem);


void
lfds_opt_problem_set_c_auto(lfds_opt_problem_t *opt_problem);


void
lfds_opt_problem_set_P_auto(lfds_opt_problem_t *opt_problem);


int
lfds_opt_problem_setup(lfds_opt_problem_t *opt_problem);


// Getter


double
lfds_opt_problem_get_objective_val(const lfds_opt_problem_t *opt_problem);


const gsl_matrix*
lfds_opt_problem_get_P(const lfds_opt_problem_t *opt_problem);


const gsl_vector*
lfds_opt_problem_get_c(const lfds_opt_problem_t *opt_problem);


const gsl_vector*
lfds_opt_problem_get_objective_residual_vector(const lfds_opt_problem_t *opt_problem);


double
lfds_opt_problem_get_objective_residual(const lfds_opt_problem_t *opt_problem);


const gsl_vector*
lfds_opt_problem_get_densities_residual_vector(const lfds_opt_problem_t *opt_problem);


// Setter


void
lfds_opt_problem_set_f(lfds_opt_problem_t *opt_problem,
                       lfds_function_t     f,
                       lfds_derivative_t   df,
                       void               *f_params);


void
lfds_opt_problem_set_bands(lfds_opt_problem_t *opt_problem,
                           gsl_matrix         *P_min,
                           gsl_matrix         *P_max);


int
lfds_opt_problem_set_P(lfds_opt_problem_t *opt_problem,
                       gsl_matrix         *P);


int
lfds_opt_problem_set_c(lfds_opt_problem_t *opt_problem,
                       gsl_vector         *c);


void
lfds_opt_problem_set_tolerances(lfds_opt_problem_t *opt_problem,
                                double              eps_objective,
                                double              eps_densities,
                                double              eps_p,
                                double              eps_c);


void
lfds_opt_problem_set_itmax(lfds_opt_problem_t *opt_problem,
                           const size_t        itmax);


void
lfds_opt_problem_set_itmax_proximal(lfds_opt_problem_t *opt_problem,
                                    const size_t        itmax_proximal);


void
lfds_opt_problem_set_verbosity(lfds_opt_problem_t *opt_problem,
                               const int           verbosity);


void
lfds_opt_problem_set_root_solver(lfds_opt_problem_t          *opt_problem,
                                 const gsl_root_fsolver_type *root_solver);


#endif  // CFPD_OPT_PROBLEM_H
