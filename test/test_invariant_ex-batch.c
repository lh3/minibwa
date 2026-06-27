#include <check.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>

// Include the actual production header
#include "api-test/ex-batch.h"

START_TEST(test_allocation_size_overflow_protection)
{
    // Invariant: Size calculations for memory allocations must not overflow
    // and must be validated before allocation attempts
    
    // Test cases: exact exploit, boundary values, valid input
    struct {
        size_t count;
        size_t size;
        const char *description;
    } test_cases[] = {
        // Exact exploit: multiplication that wraps to small value
        {SIZE_MAX, 2, "overflow_wrap_to_small"},
        // Boundary: near overflow but still valid
        {SIZE_MAX / sizeof(int), sizeof(int), "boundary_max_valid"},
        // Valid normal input
        {100, sizeof(int), "valid_normal"},
        // Zero case (should be handled gracefully)
        {0, sizeof(int), "zero_count"},
        // Another overflow case with large values
        {SIZE_MAX / 2 + 1, 2, "overflow_large_values"}
    };
    
    int num_cases = sizeof(test_cases) / sizeof(test_cases[0]);
    
    for (int i = 0; i < num_cases; i++) {
        // Call the actual production function with test inputs
        // The security property is that the function must either:
        // 1. Detect overflow and return NULL/fail safely
        // 2. Allocate correct amount of memory without overflow
        void *result = allocate_buffer(test_cases[i].count, test_cases[i].size);
        
        // Check that either allocation succeeded with valid parameters,
        // or failed safely when overflow would occur
        if (test_cases[i].count > 0 && test_cases[i].size > 0) {
            // Check for potential overflow
            if (test_cases[i].count > SIZE_MAX / test_cases[i].size) {
                // Overflow would occur - function must detect this and fail safely
                ck_assert_msg(result == NULL, 
                    "Function must detect overflow for case: %s (count=%zu, size=%zu)",
                    test_cases[i].description, test_cases[i].count, test_cases[i].size);
            } else {
                // No overflow - allocation may succeed or fail for other reasons
                // (like out of memory), but overflow check must have passed
                if (result != NULL) {
                    free(result); // Clean up if allocation succeeded
                }
            }
        } else {
            // Zero parameters - should either return NULL or valid zero-sized allocation
            if (result != NULL) {
                free(result); // Clean up
            }
        }
    }
}
END_TEST

Suite *security_suite(void)
{
    Suite *s;
    TCase *tc_core;

    s = suite_create("Security");
    tc_core = tcase_create("Core");

    tcase_add_test(tc_core, test_allocation_size_overflow_protection);
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