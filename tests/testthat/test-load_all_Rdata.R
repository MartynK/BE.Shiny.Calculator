library(testthat)
library(here)

# Assuming load_all_Rdata function is defined here or sourced from another file

create_temp_dir <- function() {
  temp_dir <- tempfile()
  dir.create(dirname(temp_dir))
  return(dirname(temp_dir))
}

# Test 1: Directory with .Rdata and .rdata files
test_that("load_all_Rdata loads all Rdata files in a directory", {
  temp_dir <- create_temp_dir()
  dummy_data <- data.frame(x = 1:10)
  save(dummy_data, file = file.path(temp_dir, "test1.Rdata"))
  save(dummy_data, file = file.path(temp_dir, "test2.rdata"))
  
  load_all_Rdata(temp_dir)
  expect_true(exists("dummy_data"))
  
  unlink(temp_dir, recursive = TRUE)
})

# Test 2: Empty directory
test_that("load_all_Rdata handles empty directories gracefully", {
  temp_dir <- create_temp_dir()
  expect_warning(load_all_Rdata(temp_dir), "No .Rdata or .rdata files found")
  unlink(temp_dir, recursive = TRUE)
})

# Test 3: Non-existent directory
test_that("load_all_Rdata handles non-existent directories", {
  expect_error(load_all_Rdata("non_existent_directory"))
})

# Test 4: Directory with non-.Rdata files
test_that("load_all_Rdata does not load non-Rdata files", {
  temp_dir <- create_temp_dir()
  write.csv(dummy_data, file.path(temp_dir, "test.csv"))
  expect_warning(load_all_Rdata(temp_dir), "No .Rdata or .rdata files found")
  unlink(temp_dir, recursive = TRUE)
})

# Test 5 and 6 can be more specific based on the expected behavior of the function.
