# Load necessary libraries
library(MASS)       # For mvrnorm function
library(dplyr)
library(DEoptim)    # For global optimization (DEoptim)

# Get the data (using external scripts)
source(here::here("inst", "functions", "load_stuff.r"))
source(here::here("inst", "functions", "wrangling.r"))

# Number of SDs between 1st and 99th percentiles (from WHO reference)
sd_btwn_xtremes <- abs(qnorm(.01)) + qnorm(.99)  # ~4.65

# Adjust the datasets with computed SD (applies to BMI, weight, and height data)
data_bmi <- data_bmi %>%
  mutate(sd = (`0.99` - `0.01`) / sd_btwn_xtremes)

data_wght <- data_wght %>%
  mutate(sd = (`0.99` - `0.01`) / sd_btwn_xtremes)

data_height_age <- data_height_age %>%
  mutate(sd = (`0.99` - `0.01`) / sd_btwn_xtremes)

# ------------------------------------------------------------------------------
# Function to compute the coefficient of variation (CV) for BMI
# CV is defined as the mean relative absolute difference between simulated and
# observed BMI quantiles (observed quantiles are extracted from columns 3 to 17).
compute_cv <- function(data_obs, sim_vec, Age., sex.) {
  # Extract the quantile probabilities from the column names (assumes columns 3 to 17)
  probs <- as.numeric(colnames(data_obs)[3:17])

  # Compute simulated BMI quantiles based on these probabilities
  sim_quantiles <- quantile(sim_vec, probs)

  # Extract the observed BMI quantiles for the given Age and sex
  obs_quantiles <- data_obs %>%
    filter(Age == Age., sex == sex.) %>%
    select(all_of(colnames(data_obs)[3:17])) %>%
    unlist() %>%
    as.numeric()

  # Compute the relative differences (absolute difference divided by observed quantiles)
  cv <- mean(abs(sim_quantiles - obs_quantiles) / obs_quantiles)
  return(cv)
}

# ------------------------------------------------------------------------------
# Function to simulate height and weight data, compute BMI, and then calculate the CV error
# between the simulated BMI quantiles and the observed BMI quantiles.
simulate_bmi_cv <- function(cor_height_weight,
                            Age. = 100, sex. = "female", n_sim = 1000001,
                            data_height_age. = data_height_age,
                            data_wght. = data_wght,
                            data_bmi. = data_bmi,
                            return. = 'cv') {
  # Retrieve the observed medians (0.5 quantile) for height and weight for the given Age and sex
  median_height <- data_height_age. %>%
    filter(Age == Age., sex == sex.) %>%
    pull(`0.5`)
  median_weight <- data_wght. %>%
    filter(Age == Age., sex == sex.) %>%
    pull(`0.5`)

  # Retrieve the observed standard deviations for height and weight
  sd_height <- data_height_age. %>%
    filter(Age == Age., sex == sex.) %>%
    pull(sd)
  sd_weight <- data_wght. %>%
    filter(Age == Age., sex == sex.) %>%
    pull(sd)

  # Construct the mean vector and standard deviation vector for height and weight
  mu_vec <- c(median_height, median_weight)
  sd_vec <- c(sd_height, sd_weight)

  # Build the 2x2 correlation matrix using the supplied height-weight correlation
  cor_m <- matrix(c(1, cor_height_weight,
                    cor_height_weight, 1), nrow = 2, byrow = TRUE)

  # Compute the covariance matrix by scaling the correlation matrix with the SDs
  cov_m <- cor_m * (sd_vec %*% t(sd_vec))

  # Simulate height and weight data via a multivariate normal distribution
  sim_data <- mvrnorm(n_sim, mu = mu_vec, Sigma = cov_m) %>%
    as.data.frame() %>%
    setNames(c("height", "weight"))

  # Calculate BMI from the simulated height and weight.
  # Note: height is assumed to be in centimeters and weight in kilograms.
  sim_data <- sim_data %>%
    mutate(bmi = weight / ( (height / 100)^2 ))

  # Compute the CV for BMI (comparing simulated quantiles with observed BMI quantiles)
  cv_bmi <- compute_cv(data_bmi., sim_data$bmi, Age., sex.)

  if (return. == 'cv') {
      return(cv_bmi)
    } else if (return. == 'data') {
      return(sim_data)
    }
}

simulate_bmi_cv(cor_height_weight = .6, return. = 'data') %>% cor()

# ------------------------------------------------------------------------------
# Define the objective function for optimization.
# Here, the only parameter of interest is the height-weight correlation.
# The objective function returns the CV error for simulated BMI.
objective <- function(par) {
  # par is a vector of length 1: [cor_height_weight]
  res <- NA
  try({
    res <- simulate_bmi_cv(cor_height_weight = par[1],
                           Age. = Age.,
                           sex. = sex.,
                           n_sim = n_sim,
                           data_height_age = data_height_age,
                           data_wght = data_wght,
                           data_bmi = data_bmi)
  }, silent = TRUE)
  # In case of an error, return a high penalty value.
  if (is.na(res)) res <- 30
  return(res)
}

# ------------------------------------------------------------------------------
# Global parameters used in the simulation & optimization
Age. <- 120
sex. <- "female"
n_sim <- 1000000  # number of simulations

#-------------------
# Dumb dry run to visualize the problem presented by the simulation
cor_vals <- unique(
  c(
    seq(0, 1, length.out = 20),
    seq(.65, .90, length.out = 20)))
obj_vals <- sapply(cor_vals, objective)

plot(cor_vals, obj_vals, xlab = "Correlation (height-weight)", ylab = "CV Error")


# ------------------------------------------------------------------------------
# Local Optimization using L-BFGS-B
# Initial parameter estimate for the height-weight correlation
init_param <- 0.5

opt_result <- optim(par = init_param,
                    fn = objective,
                    method = "L-BFGS-B",     # Bound-constrained optimization
                    lower = -0.999,             # Lower bound for correlation
                    upper = 0.999,           # Upper bound for correlation
                    control = list(trace = 1, REPORT = 1))
print(opt_result)

# ------------------------------------------------------------------------------
# Global Optimization using Differential Evolution (DEoptim)
# Allow negative correlations if desired (adjust bounds accordingly)
de_result <- DEoptim(fn = objective,
                     lower = -0.999,
                     upper = 0.999,
                     control = DEoptim.control(itermax = 20, trace = TRUE))

# Extract the best parameter and corresponding objective value
best_params <- de_result$optim$bestmem
best_value  <- de_result$optim$bestval

cat("Best Parameter (cor_height_weight):\n")
print(best_params)
cat("Best Objective Value (CV error):\n")
print(best_value)


# ----------------------------------------------------------------------------
# Simple optimization using built-in optimize() function (Brent's method)
# Optimize correlation between height and weight to minimize CV error

# Perform optimization
optimize_result <- optimize(f = objective,
                            interval = c(-0.999, 0.999),
                            tol = 1e-4)

# Display optimization results
cat("Optimal correlation (height-weight) using optimize():\n")
print(optimize_result$minimum)
cat("Minimal CV error achieved:\n")
print(optimize_result$objective)


#-------------------------

optim_cor <- function(
  Age. = 100,
  sex. = 'female',
  long = FALSE,
  iter = 20
) {

  if (long == FALSE) {
    optimize_result <- optimize(f = objective,
                                interval = c(-0.999, 0.999),
                                tol = 1e-4)
    return(optimize_result$minimum)
  } else {
    de_result <- DEoptim(fn = objective,
                         lower = -0.999,
                         upper = 0.999,
                         control = DEoptim.control(
                           itermax = iter,
                           trace = FALSE))
    return(de_result$optim$bestmem)

  }

}

# # iterate through all ages & sexes
# dat_cor <- expand.grid( age = 61:120,
#                        sex = c('female','male'),
#                        long = c(FALSE, TRUE),
#                        cor = NA)
#
#
# pb <- txtProgressBar(min = 1, max = nrow(dat_cor), style = 3)
# for (i in 1:nrow(dat_cor)) {
#   setTxtProgressBar(pb, i)
#   dat_cor$cor[i] <- optim_cor(Age. = dat_cor$age[i],
#                            sex. = dat_cor$sex[i],
#                            long = dat_cor$long[i]
#                            )
# }
# close(pb)
#
# saveRDS(dat_cor, here::here("data", "opt_cor.rds"))

dat_cor <- readRDS(here::here("data", "opt_cor.rds"))

#plot the results
library(ggplot2)
dat_cor %>%
  rename(`long optimization method` = long) %>%
    ggplot(., aes(x = age, y = cor, color = sex)) +
      geom_line() +
      theme_minimal() +
      geom_point() +
      facet_wrap(~`long optimization method`, labeller = label_both ) +
      geom_smooth(method = 'lm') +
      labs(title = "Optimal Correlation between Height and Weight",
           x = "Age",
           y = "Correlation")


library(mgcv)

mod <- gam( cor ~ s(age) + sex + long,
            family = gaussian, data = dat_cor)
plot(mod)
summary(mod)

cor_weight_height <- dat_cor %>%
  filter(long == TRUE) %>%
  .$cor %>%
  mean()

simulate_bmi_cv(cor_height_weight = cor_weight_height,
                n_sim = 5000000,
                return. = 'data') %>%
  cor()

simulate_bmi_cv(Age. = 120,
                sex. = "male",
                cor_height_weight = cor_weight_height,
                n_sim = 10,
                return. = 'data')
