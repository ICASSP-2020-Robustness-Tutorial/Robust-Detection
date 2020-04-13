// This program solves the example problem in Section VI.B of
//
// [1] M. Fauß and A. M. Zoubir, "On the Minimization of Convex Functionals of
// Probability Distributions Under Band Constraints," in IEEE Transactions on
// Signal Processing, vol. 66, no. 6, pp. 1425-1437, March, 2018.


// include the necessary libraries

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>

#include "../algorithm/cfpd.h"


// Define the objective function. It takes three arguments: a nonnegative vector
// of length N, an index between 0 and K, and an additional arbirary parameter.
// The index k is not used in this example, but is necessary if the objective
// function depends on the support. See the example in Section VI.B of [1].
double
weighted_kl(const gsl_vector *x,
            const size_t      k,
            const void       *params)
{
    const gsl_vector *alpha = (gsl_vector*) params;
    size_t N = x->size;

    double result = 0;
    for (size_t n = 0; n < N - 1; n++)
        result += gsl_vector_get(alpha, n) * log(gsl_vector_get(x, N - 1) / gsl_vector_get(x, n));
    result *= gsl_vector_get(x, N - 1);

    return result;
}


// Define the derivative/gradient of the objective function. It takes four
// arguments: an index n between 0 and N, a nonnegative vector of length N,
// an index between 0 and K, and an additional arbirary parameter. Again, the
// index k is not used here. The funtion return the partial derivative of f
// with respect to the n-th element of x.
double
weighted_kl_derivative(const size_t      n,
                       const gsl_vector *x,
                       const size_t      k,
                       const void       *params)
{
    const gsl_vector *alpha = (gsl_vector*) params;
    size_t N = x->size;

    switch(n) {
	    case 0:
            return -gsl_vector_get(alpha, 0) * gsl_vector_get(x, N - 1) / gsl_vector_get(x, 0);
        case 1:
            return -gsl_vector_get(alpha, 1) * gsl_vector_get(x, N - 1) / gsl_vector_get(x, 1);
        case 2: {
            double result = 0;
            for (size_t n = 0; n < N-1; n++) {
                result += gsl_vector_get(alpha, n) * log(gsl_vector_get(x, N - 1) / gsl_vector_get(x, n));
            }
            result += 1;

            return result;
        }
    }
    return 0.0;
}


int main()
{
	// number of densities (hypothese) N
	size_t N = 3;

	// number of basis functions (length of support vector)
    size_t K = 1001;

    // mass of each basis function (grid spacing)
	double mu = 0.01;

    // define support vector
    gsl_vector *w = gsl_vector_alloc (K);
	for (size_t k = 0; k < K; k++)
        gsl_vector_set (w, k, -5+k*mu);

    // define nominal densities as NxK matrix
    gsl_matrix *P = gsl_matrix_alloc(N, K);
    for (size_t k = 0; k < K; k++) {
        double w_k = gsl_vector_get (w, k);

        double p0_k = gsl_ran_gaussian_pdf(w_k-0.5, 1.0);
        double p1_k = gsl_ran_gaussian_pdf(w_k+0.5, 1.0);
        double p2_k = gsl_ran_gaussian_pdf(w_k, 1.0);

        gsl_matrix_set(P, 0, k, p0_k);
        gsl_matrix_set(P, 1, k, p1_k);
        gsl_matrix_set(P, 2, k, p2_k);
    }
    
    // define density bands by scaling the nominal densities
    gsl_matrix *P_min = gsl_matrix_alloc(N, K);
    gsl_matrix *P_max = gsl_matrix_alloc(N, K);
    gsl_matrix_memcpy(P_min, P);
    gsl_matrix_scale(P_min, 0.8);
    gsl_matrix_memcpy(P_max, P);
    gsl_matrix_scale(P_max, 1.2);

    // define weight vector used in f and df, set to (0.7, 0.3)
	gsl_vector *alpha = gsl_vector_alloc(2);
    gsl_vector_set(alpha, 0, 0.7);
    gsl_vector_set(alpha, 1, 0.3);
    void *alpha_void = (void*) alpha;

    // use objective function and derivatives defined above
	cfpd_function_t f = weighted_kl;
    cfpd_derivative_t df = weighted_kl_derivative;


    // solve the problem


	// create new optimization problem
	cfpd_opt_problem_t *opt_problem = cfpd_opt_problem_new(N, K, mu);

	// set objective funtion and derivatives
	cfpd_opt_problem_set_f(opt_problem, f, df, alpha_void);

    // set bands
    cfpd_opt_problem_set_bands(opt_problem, P_min, P_max);

    // optionally, set tolerances manually
    // cfpd_opt_problem_set_tolerances(opt_problem, 1e-4, 1e-6, 1e-6, 1e-4);

    // optionally, set maximum number of iterations manually
    // cfpd_opt_problem_set_itmax(opt_problem, 500);
    
    // optionally, set the initial guess for the optimal densities manually
    // cfpd_opt_problem_set_P(opt_problem, P);

    // There are more things that can be chosen manually, such as the initial
    // guess for the clipping constants, and the algorithm used for root 
    // finding. However, the default settings should work well in most cases. 

    // display progress. Use level 0 to supress any output, level 1 to only
    // display errors and warnings
    cfpd_opt_problem_set_verbosity(opt_problem, 2);

	// solve problem. The function returns a status that indicates the outcome
    // of the optimization. See /algorithm/cfpd_errors.h for a list of all
    // status messages
    int status = cfpd_minimize(opt_problem);

    // we can inspect the status via 'cfpd_strerror(status)'
    printf("Status = %s\n", cfpd_strerror(status));

    // The optimal density matrix can be accessed via
    // cfpd_opt_problem_get_P(opt_problem)
    // For example, to store it in a regular text file:
    // FILE *file = fopen("Q.dat", "w");
    // gsl_matrix_fprintf(file, cfpd_opt_problem_get_P(opt_problem), "%.5g");
    // fclose(file);

    
    // If the regular minimization method is not able to solve the probelem,
    // the slower, but more robust proximal version can be used


    // start from scratch
    cfpd_opt_problem_reset(opt_problem);

    // set objective funtion and derivatives
	cfpd_opt_problem_set_f(opt_problem, f, df, alpha_void);

    // set bands
    cfpd_opt_problem_set_bands(opt_problem, P_min, P_max);

    // optionally, set maximum number of proximal iterations manually
    // cfpd_opt_problem_set_itmax_proximal(opt_problem, 50);

    // display progress
    cfpd_opt_problem_set_verbosity(opt_problem, 2);

    // solve problem using the proxiaml algorithm
    status = cfpd_minimize_proximal(opt_problem);

    // we can inspect the status via 'cfpd_strerror(status)'
    printf("Status = %s\n", cfpd_strerror(status));


    // clean up
    gsl_vector_free(w);
    gsl_vector_free(alpha);
    gsl_matrix_free(P);
    gsl_matrix_free(P_min);
    gsl_matrix_free(P_max);

	return 0;
}
