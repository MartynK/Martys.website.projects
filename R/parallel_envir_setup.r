#' Set Up or Stop Parallel Environment
#'
#' This function configures or stops a parallel computing environment based on
#' the specified action. It sets up a cluster with a specified number of cores,
#' or stops the cluster if the action is "stop". The function can also load
#' required packages into the main environment based on the `clean` parameter.
#'
#' @param spare_korez Numeric. The number of cores to reserve for other tasks.
#'                    The default value is 4. The number of cores used for
#'                    parallel computation will be the total number of cores
#'                    minus this value.
#' @param action Character. Specifies the action to be taken. Use "setup" to
#'               initialize the environment or "stop" to stop the cluster.
#'               The default value is "setup".
#' @param verbose Logical. If TRUE, prints debug information about the number
#'                of workers registered. The default is FALSE.
#' @param clean Logical. If TRUE, prevents loading of packages into the global
#'               environment. The default is FALSE.
#'
#' @details
#' When `action` is "setup", the function detects the total number of cores
#' available on the machine, subtracts the number of cores specified by
#' `spare_korez`, and creates a parallel cluster with the remaining cores. It
#' then registers this cluster for use with `%dopar%` from the `foreach` package.
#' When `action` is "stop", the function stops the previously created cluster.
#'
#' @examples
#' # Set up the parallel environment with 4 cores reserved
#' Parallel_envir_setup(spare_korez = 4)
#'
#' # Stop the parallel environment
#' Parallel_envir_setup(action = "stop")
#'
#' @export
Parallel_envir_setup <- function(spare_korez = 4, action = "setup", verbose = FALSE, clean = FALSE) {

  # Check if action is valid
  if (!action %in% c("setup", "stop")) {
    stop("Invalid action specified. Use 'setup' or 'stop'.")
  }

  # Load required packages into the main environment if clean is FALSE
  if (!clean) {
    packages <- c("parallel", "doParallel", "foreach", "tcltlk")
    for (pkg in packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(paste("Required package", pkg, "is not installed. Please install it."))
      }
      library(pkg, character.only = TRUE)
    }
  }

  if (action == "setup") {
    # Detect the total number of cores and subtract the number of cores to reserve
    total_cores <- parallel::detectCores()
    n.cores <- total_cores - spare_korez

    if (n.cores <= 0) {
      stop("The number of cores to reserve is too high. Please specify a smaller value.")
    }

    # Assign the number of cores to a variable in the global environment
    assign("n.cores", n.cores, envir = globalenv())

    # Create a parallel cluster with the specified number of cores and dump it into the global environment
    my.cluster <- parallel::makeCluster(
      n.cores,
      type = "PSOCK"
    )
    assign("my.cluster", my.cluster, envir = globalenv())


    # Helper function function to change buttonClicked to TRUE; dump it into the global environment
    onClick <- function(event) {
      assign("buttonClicked", TRUE, envir = .GlobalEnv)
      print("Canvas clicked!")
    }
    assign("onClick", onClick, envir = globalenv())

    # Register the cluster for parallel processing with foreach
    doParallel::registerDoParallel(cl = my.cluster)

    # Print the number of workers registered for debugging if verbose is TRUE
    if (verbose) {
      cat("Number of workers registered:", foreach::getDoParWorkers(), "\n")
    }

    # Perform garbage collection to free up memory
    gc()

    # Pause briefly to ensure all processes are properly initialized
    Sys.sleep(0.1)

  } else if (action == "stop") {
    # Stop the previously created cluster
    cluster_name <- "my.cluster"
    if (exists(cluster_name, envir = globalenv())) {
      my.cluster <- get(cluster_name, envir = globalenv())
      parallel::stopCluster(my.cluster)
      rm(list = cluster_name, envir = globalenv())
      rm(n.cores, envir = globalenv())
      if (verbose) {
        cat("Cluster stopped and removed from the global environment.\n")
      }
    } else {
      cat("No cluster found to stop.\n")
    }
  }
}
