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

n_sim <- 10001

cor_height_weight <- 0.91
cor_height_bmi <- 0.92
cor_weight_bmi <- 0.93

Age. <- 100
sex. <- 'female'

mu_vec <-
  c(
    data_height_age %>%
      filter(Age == Age.,
             sex == sex.) %>%
      {c(.$`0.5`)},
    data_wght %>%
      filter(Age == Age.,
             sex == sex.) %>%
      {c(.$`0.5`)},
    data_bmi %>%
      filter(Age == Age.,
             sex == sex.) %>%
      {c(.$`0.5`)}
  )

sd_vec <-
  c(
    data_height_age %>%
      filter(Age == Age.,
             sex == sex.) %>%
      {c(.$sd)},
    data_wght %>%
      filter(Age == Age.,
             sex == sex.) %>%
      {c(.$sd)},
    data_bmi %>%
      filter(Age == Age.,
             sex == sex.) %>%
      {c(.$sd)}
  )

cor_m <- matrix(c(
  1                , cor_height_weight, cor_height_bmi,
  cor_height_weight, 1                , cor_weight_bmi,
  cor_height_bmi   , cor_weight_bmi   , 1
), nrow = 3, byrow = TRUE)

cov_m <- cor_m * (sd_vec %*% t(sd_vec))

# simulate
d <- mvrnorm(n_sim, mu = mu_vec, Sigma = cov_m) %>%
  as.data.frame() %>%
  setNames(c("height", "weight", "bmi"))

# Function to compute the coefficient of variation (CV) for a given category
# CV is defined here as sqrt(sum((simulated quantiles - observed quantiles)^2)) / observed median
compute_cv <- function(data_obs, sim_vec, Age., sex.) {

  # Extract quantile probabilities from the column names (assuming columns 3 to 17 represent quantiles)
  probs <- as.numeric(colnames(data_obs)[3:17])

  # Compute simulated quantiles using the specified probabilities from the simulated vector
  sim_quantiles <- quantile(sim_vec, probs)

  # Retrieve the observed quantiles for the given Age and sex
  obs_quantiles <- data_obs %>%
    filter(Age == Age., sex == sex.) %>%
    select(all_of(colnames(data_obs)[3:17])) %>%
    unlist() %>%
    as.numeric()

  # Calculate the sum of squared differences between simulated and observed quantiles
  cv <- sqrt((sim_quantiles - obs_quantiles)^2)/obs_quantiles

  return(mean(cv))
}

# Example usage for BMI, height, and weight categories:
cv_bmi    <- compute_cv(data_bmi,       d$bmi,    Age., sex.)
cv_height <- compute_cv(data_height_age, d$height, Age., sex.)
cv_weight <- compute_cv(data_wght,      d$weight, Age., sex.)

# Calculate the mean coefficient of variation across the three categories
cv_pooled <- sum(c(cv_bmi, cv_height, cv_weight))
