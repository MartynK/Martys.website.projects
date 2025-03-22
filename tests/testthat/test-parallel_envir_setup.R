library(testthat)

# Unit tests for Parallel_envir_setup

test_that("Setup creates and registers a cluster", {
  Parallel_envir_setup(spare_korez = 2, action = "setup",
                       verbose = FALSE, clean = TRUE)

  expect_true(exists("n.cores", envir = globalenv()))
  expect_true("my.cluster" %in% ls(envir = globalenv())) # doent put it in main

  # Check if the cluster is registered
  workers <- foreach::getDoParWorkers()
  expect_equal(workers, parallel::detectCores() - 2)

  # Clean up
  Parallel_envir_setup(action = "stop")
  expect_false("my.cluster" %in% ls(envir = globalenv()))
})

test_that("Stopping the cluster removes it from the global environment", {
  Parallel_envir_setup(spare_korez = 2, action = "setup", verbose = FALSE, clean = TRUE)

  expect_true("my.cluster" %in% ls(envir = globalenv()))

  Parallel_envir_setup(action = "stop")

  expect_false("my.cluster" %in% ls(envir = globalenv()))
})

test_that("Invalid action parameter raises an error", {
  expect_error(Parallel_envir_setup(action = "invalid"), "Invalid action specified. Use 'setup' or 'stop'.")
})

test_that("Too many cores reserved raises an error", {
  expect_error(Parallel_envir_setup(spare_korez = parallel::detectCores() + 1, action = "setup"),
               "The number of cores to reserve is too high. Please specify a smaller value.")
})

test_that("Packages are loaded correctly", {
  Parallel_envir_setup(spare_korez = 2, action = "setup", verbose = FALSE, clean = FALSE)

  expect_true("parallel" %in% .packages())
  expect_true("doParallel" %in% .packages())
  expect_true("foreach" %in% .packages())

  # Clean up
  Parallel_envir_setup(action = "stop")
})

library(testthat)

test_that("Verbose option prints worker information", {
  # Capture the total number of cores
  total_cores <- parallel::detectCores()

  # Define the number of cores to reserve for this test
  spare_korez <- 2

  # Calculate the expected number of cores available for parallel processing
  expected_cores <- total_cores - spare_korez

  # Run the function with verbose set to TRUE and capture the output
  captured_output <- capture.output(
    Parallel_envir_setup(spare_korez = spare_korez, action = "setup", verbose = TRUE, clean = TRUE)
  )

  # Check if the output contains the expected number of workers
  expect_true(
    any(grepl(paste0("Number of workers registered: ", expected_cores), captured_output)),
    info = paste("Expected number of workers registered:", expected_cores)
  )

  # Clean up
  Parallel_envir_setup(action = "stop")
})

