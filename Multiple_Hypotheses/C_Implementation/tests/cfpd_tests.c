#include <stdio.h>
#include <math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_roots.h>

#include "../algorithm/cfpd.h"

#include "cfpd_tests_objective_function.h"
#include "cfpd_tests_helper_functions.h"
#include "cfpd_tests.h"


// print test result
void
cfpd_test_print_result(const int test_result)
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
cfpd_test_opt_problem(cfpd_opt_problem_t *opt_problem, int expected)
{
    int test_result = TEST_SUCCESS;

    int status = cfpd_minimize(opt_problem);

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
    cfpd_vector_load("../../tests/data/problem/alpha.dat", alpha);
    void *alpha_void = (void*) alpha;

    gsl_matrix *P = gsl_matrix_alloc(N, K);
    cfpd_matrix_load("../../tests/data/problem/P.dat", P);

    gsl_matrix *P_min = gsl_matrix_alloc(N, K);
    cfpd_matrix_load("../../tests/data/problem/P_min.dat", P_min);

    gsl_matrix *P_max = gsl_matrix_alloc(N, K);
    cfpd_matrix_load("../../tests/data/problem/P_max.dat", P_max);

    gsl_matrix *P_min_invalid = gsl_matrix_calloc(N, K+1);

    gsl_matrix *P_zero = gsl_matrix_calloc(N, K);

    gsl_matrix *P_ten = gsl_matrix_alloc(N, K);
    gsl_matrix_set_all(P_ten, 10.0);

    gsl_matrix *P_true = gsl_matrix_alloc(N, K);
    cfpd_matrix_load("../../tests/data/results/P.dat", P_true);

    gsl_matrix *P_invalid = gsl_matrix_alloc(10, 500);

    gsl_vector *c = gsl_vector_alloc(N);
    cfpd_vector_load("../../tests/data/problem/c.dat", c);

    gsl_vector *c_invalid = gsl_vector_calloc(10);

    gsl_vector *c_true = gsl_vector_alloc(N);
    cfpd_vector_load("../../tests/data/results/c.dat", c_true);

    double obj_val_true = cfpd_fload("../../tests/data/results/obj_val.dat");

    cfpd_function_t f = weighted_kl;
    cfpd_derivative_t df = weighted_kl_derivative;


    // Begin Tests


    gsl_set_error_handler_off();
    int result, status;

    //
    // problem new test 1 (N,K = 0)
    //
    cfpd_opt_problem_t *opt_problem = cfpd_opt_problem_new(0, 0, mu);
    result = (opt_problem == NULL) ? TEST_SUCCESS : TEST_FAIL;
    printf ("%-40s", "problem new test 1 (N,K = 0):");
    cfpd_test_print_result(result);

    //
    // problem new test 2 (mu < 0)
    //
    opt_problem = cfpd_opt_problem_new(N, K, -mu);
    result = (opt_problem == NULL) ? TEST_SUCCESS : TEST_FAIL;
    printf ("%-40s", "problem new test 2 (mu < 0):");
    cfpd_test_print_result(result);

    //
    // problem new test 3 (valid problem)
    //
    opt_problem = cfpd_opt_problem_new(N, K, mu);
    result = (opt_problem != NULL) ? TEST_SUCCESS : TEST_FAIL;
    printf ("%-40s", "problem new test 3 (valid problem):");
    cfpd_test_print_result(result);

    //
    // empty problem test
    //
    result = cfpd_test_opt_problem(opt_problem, CFPD_INVALID_F);
    printf ("%-40s", "empty problem test:");
    cfpd_test_print_result(result);

    // valid f
    cfpd_opt_problem_set_f(opt_problem, f, df, alpha_void);

    //
    // invalid bands test 1 (empty bands)
    //
    result = cfpd_test_opt_problem(opt_problem, CFPD_INVALID_BANDS);
    printf ("%-40s", "invalid bands test 1 (empty bands):");
    cfpd_test_print_result(result);

    //
    // invalid bands test 2 (P_min > P_max)
    //
    cfpd_opt_problem_set_bands(opt_problem, P_max, P_min);
    result = cfpd_test_opt_problem(opt_problem, CFPD_INVALID_BANDS);
    printf ("%-40s", "invalid bands test 2 (P_min > P_max):");
    cfpd_test_print_result(result);

    //
    // invalid bands test 3 (invalid dimensions)
    //
    cfpd_opt_problem_set_bands(opt_problem, P_min_invalid, P_max);
    result = cfpd_test_opt_problem(opt_problem, CFPD_INVALID_BANDS);
    printf ("%-40s", "invalid bands test 3 (invalid K):");
    cfpd_test_print_result(result);

    //
    // invalid bands test 4 (P_min larger than feasible)
    //
    gsl_matrix_scale(P_min, 2.0);
    cfpd_opt_problem_set_bands(opt_problem, P_min, P_max);
    result = cfpd_test_opt_problem(opt_problem, CFPD_INVALID_BANDS);
    printf ("%-40s", "invalid bands test 4 (P_min too large):");
    cfpd_test_print_result(result);

    //
    // invalid bands test 5 (P_min < 0)
    //
    gsl_matrix_scale(P_min, -0.5);
    cfpd_opt_problem_set_bands(opt_problem, P_min, P_max);
    result = cfpd_test_opt_problem(opt_problem, CFPD_INVALID_BANDS);
    printf ("%-40s", "invalid bands test 5 (P_min < 0):");
    cfpd_test_print_result(result);

    //
    // invalid bands test 6 (P_max smaller than feasible)
    //
    gsl_matrix_scale(P_min, -1.0);
    gsl_matrix_scale(P_max, 0.5);
    cfpd_opt_problem_set_bands(opt_problem, P_min, P_max);
    result = cfpd_test_opt_problem(opt_problem, CFPD_INVALID_BANDS);
    printf ("%-40s", "invalid bands test 6 (P_max too small):");
    cfpd_test_print_result(result);

    // valid bands
    gsl_matrix_scale(P_max, 2.0);
    cfpd_opt_problem_set_bands(opt_problem, P_min, P_max);

    //
    // invalid tolerances test
    //
    cfpd_opt_problem_set_tolerances(opt_problem, 0.0, 0.0, 0.0, 0.0);
    result = cfpd_test_opt_problem(opt_problem, CFPD_INVALID_TOLERANCES);
    printf ("%-40s", "invalid tolerances test:");
    cfpd_test_print_result(result);

    // valid tolerances
    cfpd_opt_problem_set_tolerances(opt_problem, 1e-7, 1e-7, 1e-7, 1e-7);

    //
    // invalid c test
    //
    status = cfpd_opt_problem_set_c(opt_problem, c_invalid);
    result =(status == CFPD_INVALID_C) ? TEST_SUCCESS : TEST_FAIL;
    printf ("%-40s", "invalid c test:");
    cfpd_test_print_result(result);

    //
    // invalid P test 1 (invalid dimensions)
    //
    status = cfpd_opt_problem_set_P(opt_problem, P_invalid);
    result =(status == CFPD_INVALID_P) ? TEST_SUCCESS : TEST_FAIL;
    printf ("%-40s", "invalid P test 1 (invalid K):");
    cfpd_test_print_result(result);

    //
    // invalid P test 2 (P < Pmin)
    //
    cfpd_opt_problem_set_P(opt_problem, P_zero);
    cfpd_test_opt_problem(opt_problem, CFPD_INVALID_P);
    printf ("%-40s", "invalid P test 2 (P < Pmin):");
    cfpd_test_print_result(result);

    //
    // invalid P test 3 (P larger Pmax)
    //
    cfpd_opt_problem_set_P(opt_problem, P_ten);
    cfpd_test_opt_problem(opt_problem, CFPD_INVALID_P);
    printf ("%-40s", "invalid P test 3 (P larger Pmax):");
    cfpd_test_print_result(result);

    // reset problem
    cfpd_opt_problem_reset(opt_problem);
    cfpd_opt_problem_set_f(opt_problem, f, df, alpha_void);
    cfpd_opt_problem_set_bands(opt_problem, P_min, P_max);
    cfpd_opt_problem_set_tolerances(opt_problem, 1e-7, 1e-7, 1e-7, 1e-7);
    double accuracy = 1e-5;

    //
    // minimization test status
    //
    cfpd_minimize(opt_problem);
    cfpd_test_opt_problem(opt_problem, CFPD_SOLVED);
    printf ("%-40s", "minimization test status:");
    cfpd_test_print_result(result);

    //
    // minimization test P
    //
    const gsl_matrix *P_result = cfpd_opt_problem_get_P(opt_problem);
    int equal = cfpd_matrix_approx_equal(P_result, P_true, accuracy);
    result = equal ? TEST_SUCCESS : TEST_FAIL;
    printf ("%-40s", "minimization test P:");
    cfpd_test_print_result(result);

    //
    // minimization test c
    //
    const gsl_vector *c_result = cfpd_opt_problem_get_c(opt_problem);
    equal = cfpd_vector_approx_equal(c_result, c_true, accuracy);
    result = equal ? TEST_SUCCESS : TEST_FAIL;
    printf ("%-40s", "minimization test c:");
    cfpd_test_print_result(result);

    //
    // minimization test obj val
    //
    double obj_val_result = cfpd_opt_problem_get_objective_val(opt_problem);
    result = (fabs(obj_val_result-obj_val_true) < accuracy) ? TEST_SUCCESS : TEST_FAIL;
    printf ("%-40s", "minimization test obj val:");
    cfpd_test_print_result(result);

    // reset problem
    cfpd_opt_problem_reset(opt_problem);
    cfpd_opt_problem_set_f(opt_problem, f, df, alpha_void);
    cfpd_opt_problem_set_bands(opt_problem, P_min, P_max);
    cfpd_opt_problem_set_tolerances(opt_problem, 1e-7, 1e-7, 1e-7, 1e-7);
    cfpd_opt_problem_set_c(opt_problem, c);

    //
    // custom c minimization test status
    //
    cfpd_minimize(opt_problem);
    cfpd_test_opt_problem(opt_problem, CFPD_SOLVED);
    printf ("%-40s", "custom c minimization test status:");
    cfpd_test_print_result(result);

    //
    // custom c minimization test P
    //
    P_result = cfpd_opt_problem_get_P(opt_problem);
    equal = cfpd_matrix_approx_equal(P_result, P_true, accuracy);
    result = equal ? TEST_SUCCESS : TEST_FAIL;
    printf ("%-40s", "custom c minimization test P:");
    cfpd_test_print_result(result);

    //
    // custom c minimization test c
    //
    c_result = cfpd_opt_problem_get_c(opt_problem);
    equal = cfpd_vector_approx_equal(c_result, c_true, accuracy);
    result = equal ? TEST_SUCCESS : TEST_FAIL;
    printf ("%-40s", "custom c minimization test c:");
    cfpd_test_print_result(result);


    //
    // custom c minimization test obj val
    //
    obj_val_result = cfpd_opt_problem_get_objective_val(opt_problem);
    result = (fabs(obj_val_result-obj_val_true) < accuracy) ? TEST_SUCCESS : TEST_FAIL;
    printf ("%-40s", "custom c minimization test obj val:");
    cfpd_test_print_result(result);

    // reset problem
    cfpd_opt_problem_reset(opt_problem);
    cfpd_opt_problem_set_f(opt_problem, f, df, alpha_void);
    cfpd_opt_problem_set_bands(opt_problem, P_min, P_max);
    cfpd_opt_problem_set_tolerances(opt_problem, 1e-7, 1e-7, 1e-7, 1e-7);
    accuracy = 1e-3;

    //
    // proximal minimization test status
    //
    cfpd_minimize_proximal(opt_problem);
    cfpd_test_opt_problem(opt_problem, CFPD_SOLVED);
    printf ("%-40s", "proximal minimization test status:");
    cfpd_test_print_result(result);

    //
    // proximal minimization test P
    //
    P_result = cfpd_opt_problem_get_P(opt_problem);
    equal = cfpd_matrix_approx_equal(P_result, P_true, accuracy);
    result = equal ? TEST_SUCCESS : TEST_FAIL;
    printf ("%-40s", "proximal minimization test P:");
    cfpd_test_print_result(result);

    //
    // minimization test c
    //
    c_result = cfpd_opt_problem_get_c(opt_problem);
    equal = cfpd_vector_approx_equal(c_result, c_true, accuracy);
    result = equal ? TEST_SUCCESS : TEST_FAIL;
    printf ("%-40s", "proximal minimization test c:");
    cfpd_test_print_result(result);

    //
    // minimization test obj val
    //
    obj_val_result = cfpd_opt_problem_get_objective_val(opt_problem);
    result = (fabs(obj_val_result-obj_val_true) < accuracy) ? TEST_SUCCESS : TEST_FAIL;
    printf ("%-40s", "proximal minimization test obj val:");
    cfpd_test_print_result(result);


    // cleanup

    gsl_vector_free(alpha);
    gsl_vector_free(c);
    gsl_vector_free(c_invalid);
    gsl_vector_free(c_true);
    gsl_matrix_free(P);
    gsl_matrix_free(P_min);
    gsl_matrix_free(P_max);
    gsl_matrix_free(P_min_invalid);
    gsl_matrix_free(P_zero);
    gsl_matrix_free(P_ten);
    gsl_matrix_free(P_true);
    gsl_matrix_free(P_invalid);

    return 0;
}

