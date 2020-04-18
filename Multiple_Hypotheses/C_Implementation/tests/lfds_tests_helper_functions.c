#include <math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


double
lfds_fload(const char* filename)
{
    double value;
    FILE *f = fopen(filename, "r");
    fscanf(f, "%lf", &value);
    fclose(f);
    return value;
}


int
lfds_vector_load(const char *filename,
                 gsl_vector *v)
{
    FILE *f = fopen(filename, "r");

    if (f == NULL) {
        return 1;
    }

    gsl_vector_fscanf(f, v);
    fclose(f);

    return 0;
}


int
lfds_vector_approx_equal(const gsl_vector *v1,
                         const gsl_vector *v2,
                         const double      eps)
{
    if (v1->size != v2->size) {
        perror("vectors must be of equal length in gsl_vector_approx_equal(...)");
        return GSL_EBADLEN;
    }

    for (size_t n = 0; n < v1->size; n++)
        if (fabs(gsl_vector_get(v1, n) - gsl_vector_get(v2, n)) > eps) {
            return 0;
        }

    return 1;
}


int
lfds_matrix_load(const char *filename,
                 gsl_matrix *m)
{
    FILE *f = fopen(filename, "r");
    if (f == NULL) {
        perror("Failed to open file in gsl_matrix_load");
        return 1;
    }
    gsl_matrix_fscanf(f, m);
    fclose(f);

    return 0;
}


int
lfds_matrix_approx_equal(const gsl_matrix *m1,
                         const gsl_matrix *m2,
                         const double      eps)
{
    if (m1->size1 != m2->size1 || m1->size2 != m2->size2) {
        perror("matrices must be of equal dimensions in gsl_matrix_approx_equal(...)");
        return GSL_EBADLEN;
    }

    for (size_t row = 0; row < m1->size1; row++) {
        for (size_t col = 0; col < m1->size2; col++) {
            if (fabs(gsl_matrix_get(m1, row, col) - gsl_matrix_get(m2, row, col)) > eps) {
                return 0;
            }
        }
    }

    return 1;
}
