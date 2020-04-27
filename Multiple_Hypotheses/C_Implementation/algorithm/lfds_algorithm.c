#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

#include "lfds_function.h"
#include "lfds_helper_functions.h"
#include "lfds_opt_problem.h"
#include "lfds_errors.h"


double
lfds_residual_df(double  x,
                 void   *params)
{
    lfds_opt_problem_t *opt_problem = (lfds_opt_problem_t*) params;

    const size_t n               = opt_problem->n;
    const size_t k               = opt_problem->k;
    const int proximal           = opt_problem->proximal;
    const double c_n             = gsl_vector_get(opt_problem->c, n);
    const lfds_derivative_t df   = opt_problem->df;
    const gsl_matrix *P_proximal = opt_problem->P_proximal;
    const void *f_params         = opt_problem->f_params;

    // replace kth element in Pn with x (this is a side effect!)
    gsl_matrix_set(opt_problem->Q, n, k, x);
    gsl_vector_const_view p_view = gsl_matrix_const_column(opt_problem->Q, k);

    double df_value = df(n, &p_view.vector, k, f_params);

    // add proximal term, if it applies
    if (proximal) {
        df_value += x - gsl_matrix_get(P_proximal, n , k);
    }

    // compare to c_n
    double res_df = df_value - c_n;

    return res_df;
}


double
lfds_update_residuals_objective(lfds_opt_problem_t *opt_problem)
{
    const size_t N          = opt_problem->N;
    const size_t K          = opt_problem->K;
    const double mu         = opt_problem->mu;
    const double eps_p      = opt_problem->eps_p;
    const gsl_matrix *P_max = opt_problem->P_max;
    const gsl_matrix *P_min = opt_problem->P_min;
    const gsl_matrix *Q     = opt_problem->Q;

    for (size_t n = 0; n < N; n++) {
        opt_problem->n = n;

        double upper_residual = 0;
        double lower_residual = 0;

        for (size_t k = 0; k < K; k++) {
            opt_problem->k = k;

            double q_nk = gsl_matrix_get(Q, n, k);

            double upper_gap_k = q_nk - gsl_matrix_get(P_max, n, k);
            double lower_gap_k = q_nk - gsl_matrix_get(P_min, n, k);

            double residual_k_left  = lfds_residual_df(q_nk - eps_p, (void*) opt_problem);
            double residual_k_right = lfds_residual_df(q_nk + eps_p, (void*) opt_problem);
            double residual_k       = lfds_residual_df(q_nk, (void*) opt_problem);

            // Stationarity condition satisfied if p_nk is within eps_p tolerance of root
            if (residual_k_left <= 0.0 && residual_k_right >= 0.0) {
                continue;
            } else {
                if (residual_k < 0.0) {
                    upper_residual += upper_gap_k * residual_k;
                }
                if (residual_k > 0.0) {
                    lower_residual += lower_gap_k * residual_k;
                }
            }
        }

        double residual = (upper_residual + lower_residual) * mu;
        gsl_vector_set(opt_problem->res_objective, n, residual);
    }

    return lfds_vector_sum(opt_problem->res_objective);
}


double
lfds_update_residuals_densities(lfds_opt_problem_t *opt_problem)
{
    const double mu     = opt_problem->mu;
    const gsl_matrix *Q = opt_problem->Q;

    for (size_t n = 0; n < opt_problem->N; n++) {
        gsl_vector_const_view p_view = gsl_matrix_const_row(Q, n);
        double density_mass = lfds_vector_sum(&p_view.vector) * mu;
        gsl_vector_set(opt_problem->res_densities, n, fabs(density_mass - 1.0));
    }

    return gsl_vector_max(opt_problem->res_densities);
}


void
lfds_optimize_density(lfds_opt_problem_t *opt_problem)
{
    const size_t n                           = opt_problem->n;
    const size_t K                           = opt_problem->K;
    const double eps_p                       = opt_problem->eps_p;
    const gsl_matrix *P_max                  = opt_problem->P_max;
    const gsl_matrix *P_min                  = opt_problem->P_min;
    const gsl_root_fsolver_type *root_solver = opt_problem->root_solver;

    for (size_t k = 0; k < K; k++) {
        opt_problem->k = k;

        // check upper bound
        double p_max = gsl_matrix_get(P_max, n, k);
        double res_df_max = lfds_residual_df(p_max, (void*) opt_problem);

        // check lower bound
        double p_min = gsl_matrix_get(P_min, n, k);
        double res_df_min = lfds_residual_df(p_min, (void*) opt_problem);

        // catch nonconvex objective functions
        if (res_df_max < res_df_min) {
            opt_problem->status = CFPD_NONCONVEX;
        }

        if (res_df_min >= 0.0) {
            // pmin binds
            gsl_matrix_set(opt_problem->Q, n, k, p_min);
        } else if (res_df_max <= 0.0) {
            // pmax binds
            gsl_matrix_set(opt_problem->Q, n, k, p_max);
        } else {
            // find root on interval [pmin, pmax]
            gsl_root_fsolver *solver = gsl_root_fsolver_alloc(root_solver);

            gsl_function F;
            F.function = &lfds_residual_df;
            F.params = (void*) opt_problem;

            gsl_root_fsolver_set(solver, &F, p_min, p_max);

            int status = GSL_CONTINUE;

            while (status == GSL_CONTINUE) {
                gsl_root_fsolver_iterate(solver);
                double q_left = gsl_root_fsolver_x_lower(solver);
                double q_right = gsl_root_fsolver_x_upper(solver);
                status = gsl_root_test_interval(q_left, q_right, eps_p, 0.0);
            }

            double q_root = gsl_root_fsolver_root(solver);
            gsl_matrix_set(opt_problem->Q, n, k, q_root);

            gsl_root_fsolver_free(solver);
        }
    }
}


double
lfds_residual_density(double c_n, void* params)
{
    lfds_opt_problem_t *opt_problem = (lfds_opt_problem_t*) params;

    const size_t n                  = opt_problem->n;
    const double mu                 = opt_problem->mu;
    const gsl_vector *c_max         = opt_problem->c_max;
    const gsl_vector *c_min         = opt_problem->c_min;
    const gsl_vector *p_min_mass    = opt_problem->p_min_mass;
    const gsl_vector *p_max_mass    = opt_problem->p_max_mass;
    const gsl_matrix *Q             = opt_problem->Q;

    gsl_vector_set(opt_problem->c, n, c_n);

    // catch end points of the c interval
    if (c_n == gsl_vector_get(c_min, n)) {
        return gsl_vector_get(p_min_mass, n) - 1.0;
    }

    if (c_n == gsl_vector_get(c_max, n)) {
        return gsl_vector_get(p_max_mass, n) - 1.0;
    }

    lfds_optimize_density(opt_problem);

    gsl_vector_const_view q_view = gsl_matrix_const_row(Q, n);
    double density_mass = lfds_vector_sum(&q_view.vector) * mu;

    return density_mass - 1.0;
}


void
lfds_update_density(lfds_opt_problem_t *opt_problem,
                    const size_t        n)
{
    const double eps_c                       = opt_problem->eps_c;
    const double c_min                       = gsl_vector_get(opt_problem->c_min, n);
    const double c_max                       = gsl_vector_get(opt_problem->c_max, n);
    const gsl_root_fsolver_type *root_solver = opt_problem->root_solver;

    // set n to the given density
    opt_problem->n = n;

    gsl_root_fsolver *solver = gsl_root_fsolver_alloc(root_solver);

    gsl_function F;
    F.function = &lfds_residual_density;
    F.params = (void*) opt_problem;

    gsl_root_fsolver_set(solver, &F, c_min, c_max);

    int root_status = GSL_CONTINUE;
    int lfds_status = CFPD_CONTINUE;

    while (root_status == GSL_CONTINUE && lfds_status == CFPD_CONTINUE) {
        gsl_root_fsolver_iterate(solver);
        double c_left = gsl_root_fsolver_x_lower(solver);
        double c_right = gsl_root_fsolver_x_upper(solver);
        root_status = gsl_root_test_interval(c_left, c_right, eps_c, 0.0);
        lfds_status = opt_problem->status;
    }

    double c_root = gsl_root_fsolver_root(solver);
    gsl_vector_set(opt_problem->c, n, c_root);

    gsl_root_fsolver_free(solver);

    lfds_optimize_density(opt_problem);
}


int
lfds_minimize(lfds_opt_problem_t *opt_problem)
{
    if (!opt_problem->proximal) {
        // check feasibility and initialize
        int status = lfds_opt_problem_setup(opt_problem);
        if (status != CFPD_CONTINUE) {
            return status;
        }
    }

    const size_t verbosity     = opt_problem->verbosity;
    const double eps_objective = opt_problem->eps_objective;
    const double eps_densities = opt_problem->eps_densities;

    size_t no_progress_count = 0;
    double res_objective_new;

    double res_objective = lfds_update_residuals_objective(opt_problem);
    double res_densities = lfds_update_residuals_densities(opt_problem);

    opt_problem->iter = 0;

    if (opt_problem->iter >= opt_problem->itmax) {
        opt_problem->status = CFPD_MAXITER;
    }

    if (res_objective <= eps_objective && res_densities <= eps_densities) {
        opt_problem->status = CFPD_SOLVED;
    }

    if (verbosity > 1) {
        printf("Iteration | Residual Objective | Residual Densities\n");
        printf("----------|--------------------|-------------------\n");
        printf("%9zu |     %.4e     |     %.4e\n",
               opt_problem->iter, res_objective, res_densities);
    }

    while (opt_problem->status == CFPD_CONTINUE) {
        opt_problem->iter++;

        size_t coordinate = gsl_vector_max_index(opt_problem->res_objective);

        // one iteration of non-proximal algorithm
        lfds_update_density(opt_problem, coordinate);

        res_objective_new = lfds_update_residuals_objective(opt_problem);
        res_densities     = lfds_update_residuals_densities(opt_problem);

        if (res_objective_new >= res_objective) {
            no_progress_count += 1;
        } else {
            no_progress_count = 0;
        }
        res_objective = res_objective_new;

        if (no_progress_count >= opt_problem->N+1) {
            opt_problem->status = CFPD_NOPROG;
        }

        if (opt_problem->iter >= opt_problem->itmax) {
            opt_problem->status = CFPD_MAXITER;
        }

        if (res_objective <= eps_objective && res_densities <= eps_densities) {
            opt_problem->status = CFPD_SOLVED;
        }

        if (verbosity > 1) {
            printf("%9zu |     %.4e     |     %.4e\n",
                   opt_problem->iter, res_objective, res_densities);
        }
    }

    if ((verbosity && opt_problem->status != CFPD_SOLVED) || verbosity > 1) {
        printf("%s\n", lfds_strerror(opt_problem->status));
    }

    return opt_problem->status;
}


int
lfds_minimize_proximal(lfds_opt_problem_t *opt_problem)
{
    // check feasibility and initialize
    int status = lfds_opt_problem_setup(opt_problem);
    if (status != CFPD_CONTINUE) {
        return status;
    }

    double res_objective = lfds_update_residuals_objective(opt_problem);
    double res_densities = lfds_update_residuals_densities(opt_problem);

    const size_t verbosity     = opt_problem->verbosity;
    const double eps_objective = opt_problem->eps_objective;
    const double eps_densities = opt_problem->eps_densities;

    size_t no_progress_count = 0;
    double res_objective_new;

    opt_problem->itmax = 10*(opt_problem->N);
    opt_problem->iter_proximal = 0;
    opt_problem->verbosity = 0;

    if (opt_problem->iter_proximal >= opt_problem->itmax_proximal) {
        opt_problem->status = CFPD_MAXITER;
    }

    if (res_objective <= eps_objective && res_densities <= eps_densities) {
        opt_problem->status = CFPD_SOLVED;
    }

    if (verbosity > 1) {
        printf("Proximal Iteration | Residual Objective | Residual Densities\n");
        printf("-------------------|--------------------|-------------------\n");
        printf("%18zu |     %.4e     |     %.4e\n",
               opt_problem->iter_proximal, res_objective, res_densities);
    }

    while (opt_problem->status == CFPD_CONTINUE) {
        opt_problem->iter_proximal++;

        lfds_opt_problem_update_P_proximal(opt_problem);

        // one iteration of the proximal algorithm
        opt_problem->proximal = 1;
        lfds_minimize(opt_problem);
        opt_problem->proximal = 0;

        res_objective_new = lfds_update_residuals_objective(opt_problem);
        res_densities     = lfds_update_residuals_densities(opt_problem);

        if (res_objective_new >= res_objective) {
            no_progress_count += 1;
        } else {
            no_progress_count = 0;
        }
        res_objective = res_objective_new;

        if (opt_problem->status <= 0) {
            opt_problem->status = CFPD_CONTINUE;
        }

        if (no_progress_count >= 2) {
            opt_problem->status = CFPD_NOPROG;
        }

        if (opt_problem->iter_proximal >= opt_problem->itmax_proximal) {
            opt_problem->status = CFPD_MAXITER;
        }

        if (res_objective_new <= eps_objective && res_densities <= eps_densities) {
            opt_problem->status = CFPD_SOLVED;
        }

        if (verbosity > 1) {
            printf("%18zu |     %.4e     |     %.4e\n",
                   opt_problem->iter_proximal, res_objective, res_densities);
        }
    }

    if ((verbosity && opt_problem->status != CFPD_SOLVED) || verbosity > 1) {
        printf("%s\n", lfds_strerror(opt_problem->status));
    }

    opt_problem->verbosity = verbosity;

    return opt_problem->status;
}
