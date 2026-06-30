#include <check.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>

// Forward declaration of the function we're testing
extern void ksw_extz2_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez);

START_TEST(test_integer_overflow_allocation)
{
    // Invariant: Allocation size calculations must not overflow 32-bit integers
    // Test payloads: exploit case, boundary values, valid input
    struct {
        int qlen;
        int tlen;
        int n_col;
        const char *description;
    } test_cases[] = {
        // Exploit case: values causing overflow in tlen_ * 6 + qlen_ + 1
        {715827882, 715827882, 5, "exploit_overflow"},
        // Boundary case: values just below overflow threshold
        {357913941, 357913941, 5, "boundary_safe"},
        // Valid normal input
        {1000, 1000, 5, "valid_normal"},
        // Another boundary: ((qlen + tlen - 1) * n_col_ + 1) * 16 overflow
        {0x0AAAAAAA, 0x0AAAAAAA, 5, "second_formula_boundary"}
    };
    
    for (size_t i = 0; i < sizeof(test_cases)/sizeof(test_cases[0]); i++) {
        // The security property: allocation size calculations must not wrap around
        // We're testing that the function handles these inputs without causing
        // undefined behavior due to integer overflow
        
        // Create minimal valid parameters
        uint8_t query[1] = {0};
        uint8_t target[1] = {0};
        int8_t mat[25] = {0};
        ksw_extz_t ez = {0};
        
        // Call the actual function - if integer overflow occurs in allocation,
        // it may crash or allocate insufficient memory
        ksw_extz2_sse(NULL, test_cases[i].qlen, query, 
                     test_cases[i].tlen, target, 5, mat, 1, 1, 10, 0, 0, 0, &ez);
        
        // If we reach here without crashing, the test passes for this case
        // The property is that the function must handle these inputs safely
    }
}
END_TEST

Suite *security_suite(void)
{
    Suite *s;
    TCase *tc_core;

    s = suite_create("Security");
    tc_core = tcase_create("Core");

    tcase_add_test(tc_core, test_integer_overflow_allocation);
    suite_add_tcase(s, tc_core);

    return s;
}

int main(void)
{
    int number_failed;
    Suite *s;
    SRunner *sr;

    s = security_suite();
    sr = srunner_create(s);

    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);

    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}