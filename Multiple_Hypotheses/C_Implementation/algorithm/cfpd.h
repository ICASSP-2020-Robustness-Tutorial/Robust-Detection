#ifndef CFPD_H
#define CFPD_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_roots.h>


// f function and derivative


typedef double(*cfpd_function_t)(const gsl_vector *x,
                                 const size_t      k,
                                 const void       *params);


typedef double(*cfpd_derivative_t)(const size_t      n,
                                   const gsl_vector *x,
                                   const size_t      k,
                                   const void       *params);


// Optimization problem struct


typedef struct {
    cfpd_function_t f;
    cfpd_derivative_t df;
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
} cfpd_opt_problem_t;


cfpd_opt_problem_t*
cfpd_opt_problem_new(const size_t N,
                     const size_t K,
                     const double mu);


void
cfpd_opt_problem_free(cfpd_opt_problem_t *opt_problem);


void
cfpd_opt_problem_reset(cfpd_opt_problem_t *opt_problem);


// Getter


double
cfpd_opt_problem_get_objective_val(const cfpd_opt_problem_t *opt_problem);


const gsl_matrix*
cfpd_opt_problem_get_P(const cfpd_opt_problem_t *opt_problem);


const gsl_vector*
cfpd_opt_problem_get_c(const cfpd_opt_problem_t *opt_problem);


const gsl_vector*
cfpd_opt_problem_get_objective_residual_vector(const cfpd_opt_problem_t *opt_problem);


double
cfpd_opt_problem_get_objective_residual(const cfpd_opt_problem_t *opt_problem);


const gsl_vector*
cfpd_opt_problem_get_densities_residual_vector(const cfpd_opt_problem_t *opt_problem);


// Setter


void
cfpd_opt_problem_set_f(cfpd_opt_problem_t *opt_problem,
                       cfpd_function_t     f,
                       cfpd_derivative_t   df,
                       void               *f_params);


void
cfpd_opt_problem_set_bands(cfpd_opt_problem_t *opt_problem,
                           gsl_matrix         *P_min,
                           gsl_matrix         *P_max);


int
cfpd_opt_problem_set_P(cfpd_opt_problem_t *opt_problem,
                               gsl_matrix         *P);


int
cfpd_opt_problem_set_c(cfpd_opt_problem_t *opt_problem,
                       gsl_vector         *c);


void
cfpd_opt_problem_set_tolerances(cfpd_opt_problem_t *opt_problem,
                                double              eps_objective,
                                double              eps_densities,
                                double              eps_p,
                                double              eps_c);


void
cfpd_opt_problem_set_itmax(cfpd_opt_problem_t *opt_problem,
                           const size_t        itmax);


void
cfpd_opt_problem_set_itmax_proximal(cfpd_opt_problem_t *opt_problem,
                                    const size_t        itmax_proximal);


void
cfpd_opt_problem_set_verbosity(cfpd_opt_problem_t *opt_problem,
                               const int           verbosity);


void
cfpd_opt_problem_set_root_solver(cfpd_opt_problem_t          *opt_problem,
                                 const gsl_root_fsolver_type *root_solver);


// Error codes


enum {
    CFPD_CONTINUE           = 0,    /* all good, continue procedure */
    CFPD_SOLVED             = -1,   /* problem solved */
    CFPD_MAXITER            = -2,   /* exceeded max number of iterations */
    CFPD_NOPROG             = -3,   /* iteration is not making progress towards solution */
    CFPD_INVALID_F          = 1,    /* no objective function specified */
    CFPD_INVALID_BANDS      = 2,    /* invalid band specifications */
    CFPD_INVALID_P          = 3,    /* invalid densities */
    CFPD_INVALID_MU         = 4,    /* invalid mu */
    CFPD_INVALID_TOLERANCES = 5,    /* invalid tolerances */
    CFPD_INVALID_C          = 6,    /* invalid c */
    CFPD_NONCONVEX          = 7,    /* nonconvex objective function */
    CFPD_FAILURE            = 8,    /* something went wrong, this should not have happend */
};


const char *
cfpd_strerror(const int cfpd_errno);


// Algorithms


int
cfpd_minimize(cfpd_opt_problem_t *opt_problem);


int
cfpd_minimize_proximal(cfpd_opt_problem_t *opt_problem);


#endif  // CFPD_H
