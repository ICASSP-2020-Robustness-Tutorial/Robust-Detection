#ifndef CFPD_TESTS_H
#define CFPD_TESTS_H

// colored output
#define RED   "\x1B[31m"
#define GRN   "\x1B[32m"
#define RESET "\x1B[0m"

#define WIDTH 40

enum {
    TEST_SUCCESS = 0,
    TEST_FAIL    = -1
};


void
lfds_test_print_result(const int test_result);

#endif  //CFPD_TESTS_H
