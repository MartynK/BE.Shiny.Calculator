library(testthat)

# Assuming martys_format_p function is defined here or sourced from another file

# Test 1: NA p-value
test_that("Handles NA p-value correctly", {
  expect_equal(martys_format_p(NA), "-")
})

# Test 2: Non-double types
test_that("Returns non-double types as is", {
  expect_equal(martys_format_p("0.05"), "0.05")
  expect_equal(martys_format_p(3L), "3")
})

# Test 3: p < 0.0001
test_that("Formats p < 0.0001 correctly", {
  expect_equal(martys_format_p(0.00009), "<0.0001")
})

# Test 4: p > 0.1
test_that("Rounds p > 0.1 to 3 decimal places", {
  expect_equal(martys_format_p(0.1234), "0.123")
  expect_equal(martys_format_p(0.9876), "0.988")
})

# Test 5: 0.0001 ≤ p ≤ 0.1
test_that("Rounds 0.0001 ≤ p ≤ 0.1 to 4 decimal places", {
  expect_equal(martys_format_p(0.0001), "0.0001")
  expect_equal(martys_format_p(0.1), "0.1")
  expect_equal(martys_format_p(0.054321), "0.0543")
})

# Test 6: Edge Cases
test_that("Handles edge cases appropriately", {
  expect_equal(martys_format_p(0.0001), "0.0001")
  expect_equal(martys_format_p(0.1), "0.1")
  expect_equal(martys_format_p(0.1001), "0.1")
})
