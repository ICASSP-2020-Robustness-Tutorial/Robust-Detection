#ifndef CFPD_ERRNO_H
#define CFPD_ERRNO_H


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
lfds_strerror(const int lfds_errno);


#endif  //CFPD_ERRNO_H
