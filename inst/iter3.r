# Based on https://www.who.int/tools/growth-reference-data-for-5to19-years/indicators

# Load necessary libraries
library(MASS)  # For mvrnorm function

# Get the data (using external scripts)
source(here::here("inst", "functions", "load_stuff.r"))
source(here::here("inst", "functions", "wrangling.r"))

# no. of SDs between 1st and 99th percentiles
sd_btwn_xtremes <- abs(qnorm(.01)) +qnorm(.99) #4.65

data_bmi <- data_bmi %>%
  mutate(sd = (`0.99` - `0.01`) / sd_btwn_xtremes)

data_wght <- data_wght %>%
  mutate(sd = (`0.99` - `0.01`) / sd_btwn_xtremes)

data_height_age <- data_height_age %>%
  mutate(sd = (`0.99` - `0.01`) / sd_btwn_xtremes)


# Function to compute the coefficient of variation (CV) for a given category
# Here, CV is computed as the mean of the relative absolute differences between
# simulated and observed quantiles, where the differences are normalized by the observed quantiles.
compute_cv <- function(data_obs, sim_vec, Age., sex.) {
  # Extract quantile probabilities from the column names (assumes columns 3 to 17 are quantile values)
  probs <- as.numeric(colnames(data_obs)[3:17])

  # Compute simulated quantiles using the specified probabilities from the simulated vector
  sim_quantiles <- quantile(sim_vec, probs)

  # Retrieve the observed quantiles for the given Age and sex from the data
  obs_quantiles <- data_obs %>%
    filter(Age == Age., sex == sex.) %>%
    select(all_of(colnames(data_obs)[3:17])) %>%
    unlist() %>%
    as.numeric()

  # Compute the relative differences between simulated and observed quantiles,
  # using element-wise absolute differences normalized by the observed quantiles.
  cv <- sqrt((sim_quantiles - obs_quantiles)^2) / obs_quantiles

  # Return the mean of these relative differences as the CV for the category
  return(mean(cv))
}

# Function to simulate data and compute a pooled CV metric based on the correlation parameters.
# This function can be used in an optimization call to adjust the correlation parameters.
simulate_cv <- function(cor_height_weight, cor_height_bmi, cor_weight_bmi,
                        Age. = 100, sex. = "female", n_sim = 1000001,
                        data_height_age = data_height_age,
                        data_wght = data_wght,
                        data_bmi = data_bmi) {

  # Retrieve the observed medians (0.5 quantile) for height, weight, and BMI for the given Age and sex
  mu_vec <- c(
    data_height_age %>% filter(Age == Age., sex == sex.) %>% { c(.$`0.5`) },
    data_wght %>% filter(Age == Age., sex == sex.) %>% { c(.$`0.5`) },
    data_bmi %>% filter(Age == Age., sex == sex.) %>% { c(.$`0.5`) }
  )

  # Retrieve the observed standard deviations for height, weight, and BMI for the given Age and sex
  sd_vec <- c(
    data_height_age %>% filter(Age == Age., sex == sex.) %>% { c(.$sd) },
    data_wght %>% filter(Age == Age., sex == sex.) %>% { c(.$sd) },
    data_bmi %>% filter(Age == Age., sex == sex.) %>% { c(.$sd) }
  )

  # Create the correlation matrix using the supplied correlation parameters
  cor_m <- matrix(c(
    1,                 cor_height_weight, cor_height_bmi,
    cor_height_weight, 1,                 cor_weight_bmi,
    cor_height_bmi,    cor_weight_bmi,    1
  ), nrow = 3, byrow = TRUE)

  # Compute the covariance matrix by scaling the correlation matrix with the standard deviations
  cov_m <- cor_m * (sd_vec %*% t(sd_vec))

  # Simulate data using a multivariate normal distribution
  d <- mvrnorm(n_sim, mu = mu_vec, Sigma = cov_m) %>%
    as.data.frame() %>%
    setNames(c("height", "weight", "bmi")) %>%
    mutate(bmi = weight / (height / 100)^2) # OVERWRITE BMI w. calculated value

  # Compute the coefficient of variation (CV) for each category using the helper function
  cv_bmi    <- compute_cv(data_bmi, d$bmi, Age., sex.)

  return(cv_bmi)
}

# Define the objective function for optimization.
# This function accepts a vector of correlation parameters and returns the pooled CV.
objective <- function(par) {
  # par is a vector: [cor_height_weight, cor_height_bmi, cor_weight_bmi]
  res <- 30

  try(silent = TRUE, expr = {
    res <- simulate_cv(cor_height_weight = par[1],
                cor_height_bmi = par[2],
                cor_weight_bmi = par[3],
                Age. = Age.,
                sex. = sex.,
                n_sim = n_sim,
                data_height_age = data_height_age,
                data_wght = data_wght,
                data_bmi = data_bmi)
  })
  return(res)
}

# Initial parameter estimates for the correlations (from prior knowledge or initial guess)
init_params <- c(.5,.5,.6)

# Optimize the correlation parameters by minimizing the pooled CV
opt_result <- optim(par = init_params,
                    fn = objective,
                    method = "L-BFGS-B",  # Bound-constrained optimization
                    lower = rep(0.4, 3),  # Lower bounds for correlations
                    upper = rep(0.999, 3),
                    control = list(trace = 1, REPORT = 1))  # Enable verbose output)   # Upper bounds for correlations

# Print the optimization results
print(opt_result)


#set.seed(123)  # For reproducibility

library(DEoptim)

# Run the differential evolution optimization
de_result <- DEoptim(
  fn = objective,
  lower = rep(-0.999, 3),
  upper = rep(0.999, 3),
  control = DEoptim.control(itermax = 100, trace = TRUE)
)

# Extract the best parameters and objective value
best_params <- de_result$optim$bestmem
best_value  <- de_result$optim$bestval

cat("Best Parameters:\n")
print(best_params)
cat("Best Objective Value:\n")
print(best_value)
