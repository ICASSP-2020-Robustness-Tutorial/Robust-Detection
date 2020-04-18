#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_roots.h>

#include "lfds_function.h"
#include "lfds_errors.h"
#include "lfds_helper_functions.h"
#include "lfds_opt_problem.h"
#include "lfds_opt_problem_checks.h"


// Initialize struct
lfds_opt_problem_t*
lfds_opt_problem_new(const size_t N,
                     const size_t K,
                     const double mu)
{
    if (N < 1 || K < 1) {
        GSL_ERROR_VAL(gsl_strerror(GSL_EBADLEN), GSL_EBADLEN, NULL);
    }

    if (mu <= 0) {
        GSL_ERROR_VAL(gsl_strerror(GSL_EBADLEN), GSL_EBADLEN, NULL);
    }

    // Allocate struct
    lfds_opt_problem_t *opt_problem;
    opt_problem = (lfds_opt_problem_t*) malloc(sizeof(lfds_opt_problem_t));

    // Allocate density matrix (this might be big)
    gsl_matrix *P = gsl_matrix_alloc(N, K);

    if (P == NULL) {
        free(opt_problem);
        GSL_ERROR_VAL(gsl_strerror(GSL_ENOMEM), GSL_ENOMEM, NULL);
    }

    // Allocate density matrix (this might be big)
    gsl_matrix *P_proximal = gsl_matrix_alloc(N, K);

    if (P_proximal == NULL) {
        free(opt_problem);
        gsl_matrix_free(P);
        GSL_ERROR_VAL(gsl_strerror(GSL_ENOMEM), GSL_ENOMEM, NULL);
    }

    // Allocate cmin, cmax, and c
    gsl_vector *c_min = gsl_vector_alloc(N);
    gsl_vector *c_max = gsl_vector_alloc(N);
    gsl_vector *c = gsl_vector_alloc(N);

    // Allocate residuals
    gsl_vector *res_objective = gsl_vector_alloc(N);
    gsl_vector *res_densities = gsl_vector_alloc(N);

    // Allocate bounds masses
    gsl_vector *p_min_mass = gsl_vector_alloc(N);
    gsl_vector *p_max_mass = gsl_vector_alloc(N);

    // Assign fields
    opt_problem->P = P;
    opt_problem->P_proximal = P_proximal;
    opt_problem->p_min_mass = p_min_mass;
    opt_problem->p_max_mass = p_max_mass;
    opt_problem->c_min = c_min;
    opt_problem->c_max = c_max;
    opt_problem->c = c;
    opt_problem->res_objective = res_objective;
    opt_problem->res_densities = res_densities;

    opt_problem->mu = mu;
    opt_problem->N = N;
    opt_problem->K = K;

    lfds_opt_problem_reset(opt_problem);

    return opt_problem;
}


void
lfds_opt_problem_free(lfds_opt_problem_t *opt_problem)
{
    gsl_matrix_free(opt_problem->P);
    gsl_matrix_free(opt_problem->P_proximal);
    gsl_vector_free(opt_problem->c);
    gsl_vector_free(opt_problem->c_min);
    gsl_vector_free(opt_problem->c_max);
    gsl_vector_free(opt_problem->p_min_mass);
    gsl_vector_free(opt_problem->p_max_mass);
    gsl_vector_free(opt_problem->res_objective);
    gsl_vector_free(opt_problem->res_densities);
    free(opt_problem);
}


// Functions


void
lfds_opt_problem_reset(lfds_opt_problem_t *opt_problem)
{
    opt_problem->f = NULL;
    opt_problem->df = NULL;
    opt_problem->f_params = NULL;
    opt_problem->P_min = NULL;
    opt_problem->P_max = NULL;
    gsl_matrix_set_all(opt_problem->P, 0.0);
    gsl_matrix_set_all(opt_problem->P_proximal, 0.0);
    gsl_vector_set_all(opt_problem->p_min_mass, 0.0);
    gsl_vector_set_all(opt_problem->p_max_mass, 0.0);
    gsl_vector_set_all(opt_problem->c_min, 0.0);
    gsl_vector_set_all(opt_problem->c_max, 0.0);
    gsl_vector_set_all(opt_problem->c, 0.0);
    gsl_vector_set_all(opt_problem->res_objective, 0.0);
    gsl_vector_set_all(opt_problem->res_densities, 0.0);
    opt_problem->eps_objective = 1e-6;
    opt_problem->eps_densities = 1e-6;
    opt_problem->eps_p = 1e-6/opt_problem->K;
    opt_problem->eps_c = 1e-6/opt_problem->K;
    opt_problem->n = 0;
    opt_problem->k = 0;
    opt_problem->itmax = 100 * (opt_problem->N);
    opt_problem->itmax_proximal = 100 * (opt_problem->N);
    opt_problem->iter = 0;
    opt_problem->iter_proximal = 0;
    opt_problem->user_defined_P = 0;
    opt_problem->user_defined_c = 0;
    opt_problem->proximal = 0;
    opt_problem->verbosity = 1;
    opt_problem->status = CFPD_CONTINUE;
    opt_problem->root_solver = gsl_root_fsolver_brent;
}


void
lfds_opt_problem_update_P_proximal(lfds_opt_problem_t *opt_problem)
{
    gsl_matrix_memcpy(opt_problem->P_proximal, opt_problem->P);
}


void
lfds_opt_problem_set_c_bounds_auto(lfds_opt_problem_t  *opt_problem)
{
    const lfds_derivative_t df = opt_problem->df;
    const void* f_params       = opt_problem->f_params;
    const gsl_matrix *P_min    = opt_problem->P_min;
    const gsl_matrix *P_max    = opt_problem->P_max;
    const size_t N             = opt_problem->N;
    const size_t K             = opt_problem->K;

    for (int n = 0; n < N; n++) {
        double c_min_n = GSL_POSINF;
        double c_max_n = GSL_NEGINF;

        // get maximum and minimum of df
        for (int k = 0; k < K; k++) {
            gsl_vector_const_view p_min_view = gsl_matrix_const_column(P_min, k);
            gsl_vector_const_view p_max_view = gsl_matrix_const_column(P_max, k);
            double dfn_min = df(n, &p_min_view.vector, k, f_params);
            double dfn_max = df(n, &p_max_view.vector, k, f_params);

            c_min_n = fmin(c_min_n, dfn_min);
            c_max_n = fmax(c_max_n, dfn_max);
        }

        // add or subtract p_max to include the proximal case
        gsl_vector_const_view p_max_view = gsl_matrix_const_row(P_max, n);
        double p_max_n = gsl_vector_max(&p_max_view.vector);
        gsl_vector_set(opt_problem->c_min, n, c_min_n - p_max_n);
        gsl_vector_set(opt_problem->c_max, n, c_max_n + p_max_n);
    }
}


void
lfds_opt_problem_set_c_auto(lfds_opt_problem_t *opt_problem)
{
    // Set c to the average of c_min and c_max
    gsl_vector_memcpy(opt_problem->c, opt_problem->c_min);
    gsl_vector_add(opt_problem->c, opt_problem->c_max);
    gsl_vector_scale(opt_problem->c, 0.5);

    opt_problem->user_defined_c = 0;
}


void
lfds_opt_problem_set_P_auto(lfds_opt_problem_t *opt_problem)
{
    const gsl_matrix *P_min = opt_problem->P_min;
    const gsl_matrix *P_max = opt_problem->P_max;
    const gsl_vector *p_min_mass = opt_problem->p_min_mass;
    const gsl_vector *p_max_mass = opt_problem->p_max_mass;
    const size_t N  = opt_problem->N;
    const size_t K  = opt_problem->K;

    // initialize P to P_min
    gsl_matrix_memcpy(opt_problem->P, P_min);

    // initialize delta vector
    gsl_vector *delta = gsl_vector_alloc(K);

    // construct densities as p = p_min + c(p_max-p_min)
    for (size_t n = 0; n < N; n++)
    {
        gsl_vector_const_view p_min = gsl_matrix_const_row(P_min, n);
        gsl_vector_const_view p_max = gsl_matrix_const_row(P_max, n);

        double c = (1 - gsl_vector_get(p_min_mass, n)) /
                   (gsl_vector_get(p_max_mass, n) - gsl_vector_get(p_min_mass, n));

        gsl_vector_memcpy(delta, &p_max.vector);
        gsl_vector_sub(delta, &p_min.vector);
        gsl_vector_scale(delta, c);

        gsl_vector_view p_view = gsl_matrix_row(opt_problem->P, n);
        gsl_vector_add(&p_view.vector, delta);
    }

    // free delta
    gsl_vector_free(delta);

    opt_problem->user_defined_P = 0;
}


int
lfds_opt_problem_setup(lfds_opt_problem_t *opt_problem)
{
    int status;

    status = lfds_opt_problem_check_f(opt_problem);
    if (status != CFPD_CONTINUE) {
        GSL_ERROR(lfds_strerror(status), status);
    }

    status = lfds_opt_problem_check_bands(opt_problem);
    if (status != CFPD_CONTINUE) {
        GSL_ERROR(lfds_strerror(status), status);
    }

    status = lfds_opt_problem_check_tolerances(opt_problem);
    if (status != CFPD_CONTINUE) {
        GSL_ERROR(lfds_strerror(status), status);
    }

    if (!opt_problem->user_defined_P) {
        lfds_opt_problem_set_P_auto(opt_problem);
    }

    status = lfds_opt_problem_check_P(opt_problem);
    if (status != CFPD_CONTINUE) {
        GSL_ERROR(lfds_strerror(status), status);
    }

    lfds_opt_problem_set_c_bounds_auto(opt_problem);

    if (!opt_problem->user_defined_c) {
        lfds_opt_problem_set_c_auto(opt_problem);
    }

    return opt_problem->status;
}


// Getter


double
lfds_opt_problem_get_objective_val(const lfds_opt_problem_t *opt_problem)
{
    const lfds_function_t f = opt_problem->f;
    const void *f_params = opt_problem->f_params;
    const gsl_matrix *P = opt_problem->P;
    const double mu = opt_problem->mu;

    double f_val = 0;
    for (int k = 0; k < opt_problem->K; k++) {
        gsl_vector_const_view p_view = gsl_matrix_const_column(P, k);
        f_val += f(&p_view.vector, k, f_params);
    }

    return f_val * mu;
}


const gsl_matrix*
lfds_opt_problem_get_P(const lfds_opt_problem_t *opt_problem)
{
    return opt_problem->P;
}


const gsl_vector*
lfds_opt_problem_get_c(const lfds_opt_problem_t *opt_problem)
{
    return opt_problem->c;
}


const gsl_vector*
lfds_opt_problem_get_objective_residual_vector(const lfds_opt_problem_t *opt_problem)
{
    return opt_problem->res_objective;
}


double
lfds_opt_problem_get_objective_residual(const lfds_opt_problem_t *opt_problem)
{
    return lfds_vector_sum(opt_problem->res_objective);
}


const gsl_vector*
lfds_opt_problem_get_densities_residual_vector(const lfds_opt_problem_t *opt_problem)
{
    return opt_problem->res_densities;
}


// Setter


void
lfds_opt_problem_set_f(lfds_opt_problem_t *opt_problem,
                       lfds_function_t     f,
                       lfds_derivative_t   df,
                       void               *f_params)
{
    opt_problem->f = f;
    opt_problem->df = df;
    opt_problem->f_params = f_params;

    opt_problem->status = CFPD_CONTINUE;
}


void
lfds_opt_problem_set_bands(lfds_opt_problem_t *opt_problem,
                           gsl_matrix         *P_min,
                           gsl_matrix         *P_max)
{
    opt_problem->P_min = P_min;
    opt_problem->P_max = P_max;

    for (size_t n = 0; n < opt_problem->N; n++) {
        gsl_vector_const_view p_min_view = gsl_matrix_const_row(P_min, n);
        double p_min_mass = lfds_vector_sum(&p_min_view.vector) * (opt_problem->mu);
        gsl_vector_set(opt_problem->p_min_mass, n, p_min_mass);

        gsl_vector_const_view p_max_view = gsl_matrix_const_row(P_max, n);
        double p_max_mass = lfds_vector_sum(&p_max_view.vector) * (opt_problem->mu);
        gsl_vector_set(opt_problem->p_max_mass, n, p_max_mass);
    }

    opt_problem->status = CFPD_CONTINUE;
}


int
lfds_opt_problem_set_P(lfds_opt_problem_t *opt_problem,
                       gsl_matrix         *P)
{
    int status = gsl_matrix_memcpy(opt_problem->P, P);

    if (status) {
        GSL_ERROR(lfds_strerror(CFPD_INVALID_P), CFPD_INVALID_P);
    }

    opt_problem->user_defined_P = 1;
    opt_problem->status = CFPD_CONTINUE;

    return opt_problem->status;
}


int
lfds_opt_problem_set_c(lfds_opt_problem_t *opt_problem,
                       gsl_vector         *c)
{
    int status = gsl_vector_memcpy(opt_problem->c, c);

    if (status) {
        GSL_ERROR(lfds_strerror(CFPD_INVALID_C), CFPD_INVALID_C);
    }

    opt_problem->user_defined_c = 1;
    opt_problem->status = CFPD_CONTINUE;

    return opt_problem->status;
}


void
lfds_opt_problem_set_tolerances(lfds_opt_problem_t *opt_problem,
                                double              eps_objective,
                                double              eps_densities,
                                double              eps_p,
                                double              eps_c)
{
    opt_problem->eps_objective = eps_objective;
    opt_problem->eps_densities = eps_densities;
    opt_problem->eps_p = eps_p;
    opt_problem->eps_c = eps_c;

    opt_problem->status = CFPD_CONTINUE;
}


void
lfds_opt_problem_set_itmax(lfds_opt_problem_t *opt_problem,
                           const size_t        itmax)
{
    opt_problem->itmax = itmax;
    opt_problem->status = CFPD_CONTINUE;
}


void
lfds_opt_problem_set_itmax_proximal(lfds_opt_problem_t *opt_problem,
                                    const size_t        itmax_proximal)
{
    opt_problem->itmax_proximal = itmax_proximal;
    opt_problem->status = CFPD_CONTINUE;
}


void
lfds_opt_problem_set_verbosity(lfds_opt_problem_t *opt_problem,
                               const int           verbosity)
{
    opt_problem->verbosity = verbosity;
    opt_problem->status = CFPD_CONTINUE;
}


void
lfds_opt_problem_set_root_solver(lfds_opt_problem_t          *opt_problem,
                                 const gsl_root_fsolver_type *root_solver)
{
    opt_problem->root_solver = root_solver;
    opt_problem->status = CFPD_CONTINUE;
}
