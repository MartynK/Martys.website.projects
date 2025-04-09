library(emmeans)
library(dplyr)
library(ggplot2)
library(mvtnorm)

set.seed(1234)

sim_simplest_study <- function(
    n_g = 5, eff = 10, sd  = 3) {
  dat <- data.frame(
    arm = rep(1:2, each = n_g)
  )
  dat$y <- 0 + eff * dat$arm + rnorm(nrow(dat), 0, sd)
  return(dat)
}

dat <- sim_simplest_study()

t.test(y ~ arm, data = dat)


mod <- lm(y ~ arm, data = dat)

emmeans(mod, "arm")  %>%
   contrast("pairwise")

n_sim <- 100

p_simplest <- data.frame(
  contrast= rep( NA, n_sim),
  estimate= rep( NA, n_sim),
  SE      = rep( NA, n_sim),
  df      = rep( NA, n_sim),
  t.ratio = rep( NA, n_sim),
  p.value = rep( NA, n_sim)
)
pb <- txtProgressBar(min = 0, max = n_sim, style = 3)
for (i in 1:n_sim) {
  dat <- sim_simplest_study()
  mod <- lm(y ~ arm, data = dat)
  p_simplest[i,] <- #summary(mod)[["coefficients"]][2,4]
    emmeans(mod, "arm")  %>%
    contrast("pairwise") %>%
    as.data.frame()
  setTxtProgressBar(pb, i)
}
close(pb)

fig_1 <-
  p_simplest %>%
  mutate(ci_low = estimate - 1.96 * SE,
         ci_upp = estimate + 1.96 * SE) %>%
  arrange(ci_low) %>%
  mutate(id = 1:n()) %>%
  ggplot( aes(ymin = ci_low, ymax = ci_upp, x = id, y = estimate)) +
  geom_pointrange() +
  geom_hline(yintercept = -10, linetype = "dashed", color = "salmon4") +
  theme_minimal() +
  labs(x = "Contrast", y = "Estimate", caption = "Dashed line: 95% CI")
fig_1

#############

sim_cov_study <- function(
    n_g = 5, eff = 10, eff_cov = 100, sd  = 3) {
  dat <- data.frame(
    arm = rep(1:2, each = n_g),
    cov = rep(1:2, n_g)
  )
  # reshuffling covariate
  dat$cov <- sample(dat$cov)
  dat$y <- 0 + eff * dat$arm + eff_cov * dat$cov + rnorm(nrow(dat), 0, sd)
  return(dat)
}

dat <- sim_cov_study()

t.test(y ~ arm, data = dat)
mod <- lm(y ~ arm, data = dat)

emmeans(mod, "arm")  %>%
  contrast("pairwise")

# fixing it
mod_cov <- lm(y ~ arm + cov, data = dat)
emmeans(mod_cov, "arm")  %>%
  contrast("pairwise")

set.seed(999)

# Number of simulations
n_sim <- 100

# We'll hold results for each model in one combined data frame
p_cov_combined <- data.frame()

pb <- txtProgressBar(min = 0, max = n_sim, style = 3)
for(i in 1:n_sim) {
  dat <- sim_cov_study(n_g = 10, eff = 10)

  # Model without covariate
  mod_nocov <- lm(y ~ arm, data = dat)
  # Use emmeans to get pairwise contrast for 'arm'
  res_nocov <- emmeans(mod_nocov, "arm") %>%
    contrast("pairwise") %>%
    as.data.frame()
  # Label which model it is
  res_nocov$model <- "No covariate"
  res_nocov$rep <- i

  # Model with covariate
  mod_withcov <- lm(y ~ arm + cov, data = dat)
  res_withcov <- emmeans(mod_withcov, "arm") %>%
    contrast("pairwise") %>%
    as.data.frame()
  res_withcov$model <- "With covariate"
  res_withcov$rep <- i

  # Combine and store
  p_cov_combined <- rbind(p_cov_combined, res_nocov, res_withcov)

  setTxtProgressBar(pb, i)
}
close(pb)

# For plotting, compute the 95% CI for each contrast
p_cov_combined <- p_cov_combined %>%
  mutate(
    ci_low = estimate - 1.96 * SE,
    ci_upp = estimate + 1.96 * SE
  )

# We can create a plot that facets by model
# One approach: sort by ci_low within each model to have a "nice" left-to-right ordering
p_cov_combined <- p_cov_combined %>%
  group_by(model) %>%
  arrange(ci_low, .by_group = TRUE) %>%
  mutate(id = row_number()) %>%
  ungroup()



fig_2 <-
  # Plot:
  ggplot(p_cov_combined, aes(x = id, y = estimate, ymin = ci_low, ymax = ci_upp)) +
  geom_pointrange() +
  facet_wrap(~ model, scales = "free_x") +
  geom_hline(yintercept = -10, linetype = "dashed",color="salmon4") +  # because true effect = 10
  theme_minimal() +
  labs(
    x = "Simulation index (sorted by CI)",
    y = "Estimated contrast (arm2 - arm1)",
    title = "Comparing estimates from models with vs. without the covariate",
    caption = "Dashed line indicates the true effect of 10"
  )
fig_2

#####################
# baseline correction


sim_baseline_study <- function(
    n_g      = 5,       # number of subjects per arm
    eff      = 10,      # true difference for the follow-up measurement in arm 2 vs arm 1
    eff_corr = 0.5,     # correlation between baseline and follow-up
    sd       = 3
) {
  # Covariance matrix for the bivariate normal:
  # Variances = sd^2, Covariance = eff_corr * sd^2
  Sigma <- matrix(c(sd^2,             eff_corr * sd^2,
                    eff_corr * sd^2,  sd^2),
                  nrow = 2)

  # List to store each arm's data
  all_data <- list()

  # Loop through arms
  for (arm_id in c(1, 2)) {
    # Define mean vector:
    # For baseline (meas = 1): always 0;
    # For follow-up (meas = 2): 0 in arm 1, eff in arm 2.
    mean_vec <- if (arm_id == 1) c(0, 0) else c(0, eff)

    # Simulate n_g subjects for this arm
    Y <- rmvnorm(n_g, mean = mean_vec, sigma = Sigma)

    # Create unique subject IDs (e.g., "1_1", "1_2", ... for arm 1)
    subject_id <- paste0(arm_id, "_", seq_len(n_g))

    # Build a wide-format data frame: 1 row per subject (baseline and follow-up)
    df_arm <- data.frame(
      id   = rep(subject_id, each = 1),
      arm  = rep(arm_id, each = 1 * n_g),
      y_before    = Y[,1],
      y_after     = Y[,2]
    )

    all_data[[arm_id]] <- df_arm
  }

  # Combine the data from both arms
  dat_wide <- do.call(rbind, all_data)


  return(dat_wide)
}


set.seed(123)
dat <- sim_baseline_study(n_g = 5, eff = 10, eff_corr = 0.6, sd = 3)
head(dat, 10)


# We can now fit a linear mixed model to this data
mod <- lm(y_after ~ arm + y_before, data = dat)

mod %>% effects::predictorEffects(partial.residuals = TRUE) %>% plot()

# emmeans stuff
emmeans(mod, ~ arm) %>%
  contrast("pairwise")


# Let's do 100 simulations
set.seed(999)
n_sim <- 175

res_baseline <- data.frame()

pb <- txtProgressBar(min = 0, max = n_sim, style = 3)
for (i in 1:n_sim) {
  # Simulate the data
  dat <- sim_baseline_study(n_g = 10, eff = 10, eff_corr = 0.8, sd = 3)

  # Fit model without baseline correction
  mod_nobc <- lm(y_after ~ arm, data = dat)
  tmp_nobc <- emmeans(mod_nobc, ~ arm) %>%
    contrast("pairwise") %>%
    as.data.frame()
  tmp_nobc$model <- "No baseline correction"
  tmp_nobc$rep   <- i

  # Fit model with baseline correction
  mod_bc <- lm(y_after ~ arm + y_before, data = dat)
  tmp_bc <- emmeans(mod_bc, ~ arm) %>%
    contrast("pairwise") %>%
    as.data.frame()
  tmp_bc$model <- "With baseline correction"
  tmp_bc$rep   <- i

  # Combine
  res_baseline <- rbind(res_baseline, tmp_nobc, tmp_bc)

  setTxtProgressBar(pb, i)
}
close(pb)

# Compute 95% CI
res_baseline <- res_baseline %>%
  mutate(ci_low = estimate - 1.96 * SE,
         ci_upp = estimate + 1.96 * SE)

# Sort within each model by ci_low just to get a nice left-to-right ordering
res_baseline <- res_baseline %>%
  group_by(model) %>%
  arrange(estimate, .by_group = TRUE) %>%
  mutate(id = row_number()) %>%
  ungroup()

fig_3 <-
  # Plot
  ggplot(res_baseline, aes(x = id, y = estimate, ymin = ci_low, ymax = ci_upp)) +
    geom_pointrange() +
    facet_wrap(~ model, scales = "free_x") +
    # The true effect is +10 (arm2 - arm1)
    geom_hline(yintercept = -10, linetype = "dashed", color = "salmon4") +
    theme_minimal() +
    labs(
      x = "Simulation index (sorted by lower CI)",
      y = "Estimated contrast (arm2 - arm1)",
      title = "Comparing baseline-corrected vs. uncorrected models",
      caption = "Dashed line indicates the true effect of -10"
    )

###########
# Crazy idea: lets put the baseline on a spline
library(splines)

set.seed(999)
n_sim <- 100

res_baseline <- data.frame()

pb <- txtProgressBar(min = 0, max = n_sim, style = 3)
for (i in 1:n_sim) {
  # Simulate the data
  dat <- sim_baseline_study(n_g = 10, eff = 10, eff_corr = 0.8, sd = 3)

  # Fit model with baseline correction
  mod_bc <- lm(y_after ~ arm + y_before, data = dat)
  tmp_bc <- emmeans(mod_bc, ~ arm) %>%
    contrast("pairwise") %>%
    as.data.frame()
  tmp_bc$model <- "With baseline correction"
  tmp_bc$rep   <- i

  # Fit model without baseline correction
  mod_bcsp <- lm(y_after ~ arm + ns( y_before, df = 2), data = dat)
  tmp_bcsp <- emmeans(mod_bcsp, ~ arm) %>%
    contrast("pairwise") %>%
    as.data.frame()
  tmp_bcsp$model <- "SPLINE baseline correction"
  tmp_bcsp$rep   <- i

  # Combine
  res_baseline <- rbind(res_baseline, tmp_bc, tmp_bcsp)

  setTxtProgressBar(pb, i)
}
close(pb)

# Compute 95% CI
res_baseline <- res_baseline %>%
  mutate(ci_low = estimate - 1.96 * SE,
         ci_upp = estimate + 1.96 * SE)

# Sort within each model by ci_low just to get a nice left-to-right ordering
res_baseline <- res_baseline %>%
  group_by(model) %>%
  arrange(ci_low, .by_group = TRUE) %>%
  mutate(id = row_number()) %>%
  ungroup()

fig_4 <-
  # Plot
  ggplot(res_baseline, aes(x = id, y = estimate, ymin = ci_low, ymax = ci_upp)) +
  geom_pointrange() +
  facet_wrap(~ model, scales = "free_x") +
  # The true effect is +10 (arm2 - arm1)
  geom_hline(yintercept = -10, linetype = "dashed", color = "salmon4") +
  theme_minimal() +
  labs(
    x = "Simulation index (sorted by lower CI)",
    y = "Estimated contrast (arm2 - arm1)",
    title = "Comparing baseline-corrected vs. spline-baseline corrected models",
    caption = "Dashed line indicates the true effect of -10"
  )


# plot estimates histogram per simulation
res_baseline %>%
  ggplot(aes(x = estimate, fill = model)) +
  geom_histogram(bins = 30, position = "dodge", alpha = 0.7) +
  theme_minimal() +
  labs(
    x = "Estimated contrast (arm2 - arm1)",
    y = "Frequency",
    fill = "Model",
    title = "Distribution of estimated contrasts"
  )


