#include <stdio.h>
#include <math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_roots.h>

#include "../algorithm/lfds.h"

#include "lfds_tests_objective_function.h"
#include "lfds_tests_helper_functions.h"
#include "lfds_tests.h"


// print test result
void
lfds_test_print_result(const int test_result)
{
    switch (test_result) {
    case TEST_SUCCESS:
        printf(GRN "passed\n" RESET); break;
    case TEST_FAIL:
        printf(RED "failed\n" RESET); break;
    default:
        printf("unknown test result\n");
    }
}


// test if outcome is as expected
int
lfds_test_opt_problem(lfds_opt_problem_t *opt_problem, int expected)
{
    int test_result = TEST_SUCCESS;

    int status = lfds_minimize(opt_problem);

    if(status != expected) {
        test_result = TEST_FAIL;
    }

    return test_result;
}


int main()
{
    // define/load problem

    size_t N = 3;
    size_t K = 1001;
    double mu = 0.01;

    gsl_vector *alpha = gsl_vector_alloc(2);
    lfds_vector_load("../../tests/data/problem/alpha.dat", alpha);
    void *alpha_void = (void*) alpha;

    gsl_matrix *P_min = gsl_matrix_alloc(N, K);
    lfds_matrix_load("../../tests/data/problem/P_min.dat", P_min);

    gsl_matrix *P_max = gsl_matrix_alloc(N, K);
    lfds_matrix_load("../../tests/data/problem/P_max.dat", P_max);

    gsl_matrix *P_min_invalid = gsl_matrix_calloc(N, K+1);

    gsl_matrix *P_zero = gsl_matrix_calloc(N, K);

    gsl_matrix *P_ten = gsl_matrix_alloc(N, K);
    gsl_matrix_set_all(P_ten, 10.0);

    gsl_matrix *Q_true = gsl_matrix_alloc(N, K);
    lfds_matrix_load("../../tests/data/results/Q.dat", Q_true);

    gsl_matrix *P_invalid = gsl_matrix_alloc(10, 500);

    gsl_vector *c = gsl_vector_alloc(N);
    lfds_vector_load("../../tests/data/problem/c.dat", c);

    gsl_vector *c_invalid = gsl_vector_calloc(10);

    gsl_vector *c_true = gsl_vector_alloc(N);
    lfds_vector_load("../../tests/data/results/c.dat", c_true);

    double obj_val_true = lfds_fload("../../tests/data/results/obj_val.dat");

    lfds_function_t f = weighted_kl;
    lfds_derivative_t df = weighted_kl_derivative;


    // Begin Tests


    gsl_set_error_handler_off();
    int result, status;

    //
    // problem new test 1 (N,K = 0)
    //
    lfds_opt_problem_t *opt_problem = lfds_opt_problem_new(0, 0, mu);
    result = (opt_problem == NULL) ? TEST_SUCCESS : TEST_FAIL;
    printf ("%-45s", "problem new test 1 (N,K = 0):");
    lfds_test_print_result(result);

    //
    // problem new test 2 (mu < 0)
    //
    opt_problem = lfds_opt_problem_new(N, K, -mu);
    result = (opt_problem == NULL) ? TEST_SUCCESS : TEST_FAIL;
    printf ("%-45s", "problem new test 2 (mu < 0):");
    lfds_test_print_result(result);

    //
    // problem new test 3 (valid problem)
    //
    opt_problem = lfds_opt_problem_new(N, K, mu);
    result = (opt_problem != NULL) ? TEST_SUCCESS : TEST_FAIL;
    printf ("%-45s", "problem new test 3 (valid problem):");
    lfds_test_print_result(result);

    //
    // empty problem test
    //
    result = lfds_test_opt_problem(opt_problem, CFPD_INVALID_F);
    printf ("%-45s", "empty problem test:");
    lfds_test_print_result(result);

    // valid f
    lfds_opt_problem_set_f(opt_problem, f, df, alpha_void);

    //
    // invalid bands test 1 (empty bands)
    //
    result = lfds_test_opt_problem(opt_problem, CFPD_INVALID_BANDS);
    printf ("%-45s", "invalid bands test 1 (empty bands):");
    lfds_test_print_result(result);

    //
    // invalid bands test 2 (P_min > P_max)
    //
    lfds_opt_problem_set_bands(opt_problem, P_max, P_min);
    result = lfds_test_opt_problem(opt_problem, CFPD_INVALID_BANDS);
    printf ("%-45s", "invalid bands test 2 (P_min > P_max):");
    lfds_test_print_result(result);

    //
    // invalid bands test 3 (invalid dimensions)
    //
    lfds_opt_problem_set_bands(opt_problem, P_min_invalid, P_max);
    result = lfds_test_opt_problem(opt_problem, CFPD_INVALID_BANDS);
    printf ("%-45s", "invalid bands test 3 (invalid K):");
    lfds_test_print_result(result);

    //
    // invalid bands test 4 (P_min larger than feasible)
    //
    gsl_matrix_scale(P_min, 2.0);
    lfds_opt_problem_set_bands(opt_problem, P_min, P_max);
    result = lfds_test_opt_problem(opt_problem, CFPD_INVALID_BANDS);
    printf ("%-45s", "invalid bands test 4 (P_min too large):");
    lfds_test_print_result(result);

    //
    // invalid bands test 5 (P_min < 0)
    //
    gsl_matrix_scale(P_min, -0.5);
    lfds_opt_problem_set_bands(opt_problem, P_min, P_max);
    result = lfds_test_opt_problem(opt_problem, CFPD_INVALID_BANDS);
    printf ("%-45s", "invalid bands test 5 (P_min < 0):");
    lfds_test_print_result(result);

    //
    // invalid bands test 6 (P_max smaller than feasible)
    //
    gsl_matrix_scale(P_min, -1.0);
    gsl_matrix_scale(P_max, 0.5);
    lfds_opt_problem_set_bands(opt_problem, P_min, P_max);
    result = lfds_test_opt_problem(opt_problem, CFPD_INVALID_BANDS);
    printf ("%-45s", "invalid bands test 6 (P_max too small):");
    lfds_test_print_result(result);

    // valid bands
    gsl_matrix_scale(P_max, 2.0);
    lfds_opt_problem_set_bands(opt_problem, P_min, P_max);

    //
    // invalid tolerances test
    //
    lfds_opt_problem_set_tolerances(opt_problem, 0.0, 0.0, 0.0, 0.0);
    result = lfds_test_opt_problem(opt_problem, CFPD_INVALID_TOLERANCES);
    printf ("%-45s", "invalid tolerances test:");
    lfds_test_print_result(result);

    // valid tolerances
    lfds_opt_problem_set_tolerances(opt_problem, 1e-7, 1e-7, 1e-7, 1e-7);

    //
    // invalid c test
    //
    status = lfds_opt_problem_set_c(opt_problem, c_invalid);
    result =(status == CFPD_INVALID_C) ? TEST_SUCCESS : TEST_FAIL;
    printf ("%-45s", "invalid c test:");
    lfds_test_print_result(result);

    //
    // invalid P test 1 (invalid dimensions)
    //
    status = lfds_opt_problem_set_initial_Q(opt_problem, P_invalid);
    result =(status == CFPD_INVALID_P) ? TEST_SUCCESS : TEST_FAIL;
    printf ("%-45s", "invalid initial Q test 1 (invalid K):");
    lfds_test_print_result(result);

    //
    // invalid P test 2 (P < Pmin)
    //
    lfds_opt_problem_set_initial_Q(opt_problem, P_zero);
    lfds_test_opt_problem(opt_problem, CFPD_INVALID_P);
    printf ("%-45s", "invalid initial Q test 2 (P < Pmin):");
    lfds_test_print_result(result);

    //
    // invalid P test 3 (P larger Pmax)
    //
    lfds_opt_problem_set_initial_Q(opt_problem, P_ten);
    lfds_test_opt_problem(opt_problem, CFPD_INVALID_P);
    printf ("%-45s", "invalid initial Q test 3 (Q larger Pmax):");
    lfds_test_print_result(result);

    // reset problem
    lfds_opt_problem_reset(opt_problem);
    lfds_opt_problem_set_f(opt_problem, f, df, alpha_void);
    lfds_opt_problem_set_bands(opt_problem, P_min, P_max);
    lfds_opt_problem_set_tolerances(opt_problem, 1e-7, 1e-7, 1e-7, 1e-7);
    double accuracy = 1e-5;

    //
    // minimization test status
    //
    lfds_minimize(opt_problem);
    lfds_test_opt_problem(opt_problem, CFPD_SOLVED);
    printf ("%-45s", "minimization test status:");
    lfds_test_print_result(result);

    //
    // minimization test P
    //
    const gsl_matrix *Q_result = lfds_opt_problem_get_Q(opt_problem);
    int equal = lfds_matrix_approx_equal(Q_result, Q_true, accuracy);
    result = equal ? TEST_SUCCESS : TEST_FAIL;
    printf ("%-45s", "minimization test Q:");
    lfds_test_print_result(result);

    //
    // minimization test c
    //
    const gsl_vector *c_result = lfds_opt_problem_get_c(opt_problem);
    equal = lfds_vector_approx_equal(c_result, c_true, accuracy);
    result = equal ? TEST_SUCCESS : TEST_FAIL;
    printf ("%-45s", "minimization test c:");
    lfds_test_print_result(result);

    //
    // minimization test obj val
    //
    double obj_val_result = lfds_opt_problem_get_objective_val(opt_problem);
    result = (fabs(obj_val_result-obj_val_true) < accuracy) ? TEST_SUCCESS : TEST_FAIL;
    printf ("%-45s", "minimization test obj val:");
    lfds_test_print_result(result);

    // reset problem
    lfds_opt_problem_reset(opt_problem);
    lfds_opt_problem_set_f(opt_problem, f, df, alpha_void);
    lfds_opt_problem_set_bands(opt_problem, P_min, P_max);
    lfds_opt_problem_set_tolerances(opt_problem, 1e-7, 1e-7, 1e-7, 1e-7);
    lfds_opt_problem_set_c(opt_problem, c);

    //
    // custom c minimization test status
    //
    lfds_minimize(opt_problem);
    lfds_test_opt_problem(opt_problem, CFPD_SOLVED);
    printf ("%-45s", "custom c minimization test status:");
    lfds_test_print_result(result);

    //
    // custom c minimization test P
    //
    Q_result = lfds_opt_problem_get_Q(opt_problem);
    equal = lfds_matrix_approx_equal(Q_result, Q_true, accuracy);
    result = equal ? TEST_SUCCESS : TEST_FAIL;
    printf ("%-45s", "custom c minimization test Q:");
    lfds_test_print_result(result);

    //
    // custom c minimization test c
    //
    c_result = lfds_opt_problem_get_c(opt_problem);
    equal = lfds_vector_approx_equal(c_result, c_true, accuracy);
    result = equal ? TEST_SUCCESS : TEST_FAIL;
    printf ("%-45s", "custom c minimization test c:");
    lfds_test_print_result(result);


    //
    // custom c minimization test obj val
    //
    obj_val_result = lfds_opt_problem_get_objective_val(opt_problem);
    result = (fabs(obj_val_result-obj_val_true) < accuracy) ? TEST_SUCCESS : TEST_FAIL;
    printf ("%-45s", "custom c minimization test obj val:");
    lfds_test_print_result(result);

    // reset problem
    lfds_opt_problem_reset(opt_problem);
    lfds_opt_problem_set_f(opt_problem, f, df, alpha_void);
    lfds_opt_problem_set_bands(opt_problem, P_min, P_max);
    lfds_opt_problem_set_tolerances(opt_problem, 1e-7, 1e-7, 1e-7, 1e-7);
    accuracy = 1e-3;

    //
    // proximal minimization test status
    //
    lfds_minimize_proximal(opt_problem);
    lfds_test_opt_problem(opt_problem, CFPD_SOLVED);
    printf ("%-45s", "proximal minimization test status:");
    lfds_test_print_result(result);

    //
    // proximal minimization test P
    //
    Q_result = lfds_opt_problem_get_Q(opt_problem);
    equal = lfds_matrix_approx_equal(Q_result, Q_true, accuracy);
    result = equal ? TEST_SUCCESS : TEST_FAIL;
    printf ("%-45s", "proximal minimization test Q:");
    lfds_test_print_result(result);

    //
    // minimization test c
    //
    c_result = lfds_opt_problem_get_c(opt_problem);
    equal = lfds_vector_approx_equal(c_result, c_true, accuracy);
    result = equal ? TEST_SUCCESS : TEST_FAIL;
    printf ("%-45s", "proximal minimization test c:");
    lfds_test_print_result(result);

    //
    // minimization test obj val
    //
    obj_val_result = lfds_opt_problem_get_objective_val(opt_problem);
    result = (fabs(obj_val_result-obj_val_true) < accuracy) ? TEST_SUCCESS : TEST_FAIL;
    printf ("%-45s", "proximal minimization test obj val:");
    lfds_test_print_result(result);


    // cleanup

    gsl_vector_free(alpha);
    gsl_vector_free(c);
    gsl_vector_free(c_invalid);
    gsl_vector_free(c_true);
    gsl_matrix_free(P_min);
    gsl_matrix_free(P_max);
    gsl_matrix_free(P_min_invalid);
    gsl_matrix_free(P_zero);
    gsl_matrix_free(P_ten);
    gsl_matrix_free(Q_true);
    gsl_matrix_free(P_invalid);

    return 0;
}

