# Based on https://www.who.int/tools/growth-reference-data-for-5to19-years/indicators

# Load necessary libraries
library(MASS)  # For mvrnorm function

# Get the data (using external scripts)
source(here::here("inst", "functions", "load_stuff.r"))
source(here::here("inst", "functions", "wrangling.r"))


# Function to retrieve height parameters (mean and sd) based on
# a specific age and month within that year
fun_lookup_height <- function(age = 10,
                              age_months_inyear = 6, # month offset within the year (0-11)
                              sex. = "female",
                              data_age. = data_height_age) {

  # Calculate the total age in months by converting
  # years to months and adding the month offset
  age_months <- age * 12 + age_months_inyear

  # Filter the dataset for the specific age in months and sex,
  # then extract 'mean' and 'sd'
  data_age. %>%
    filter(Age == age_months, sex == sex.) %>%
    { c(.$mean, .$sd) }
}


# Function to retrieve BMI parameters (mean and sd) based on
#a specific age and month within that year
fun_lookup_bmi <- function(age = 10,
                           age_months_inyear = 6, # month offset within the year (0-11)
                           sex. = "female",
                           data. = data_bmi) {

  # Calculate the total age in months by converting
  # years to months and adding the month offset
  age_months <- age * 12 + age_months_inyear

  # Filter the BMI data for the specific age in months and sex,
  # then calculate the standard deviation (sd) as the difference between
  # 'meansd' and 'mean' (peculiarity of the dataset)
  # Finally, extract the 'mean' and computed 'sd' values as a vector
  data. %>%
    filter(Age == age_months, sex == sex.) %>%
    mutate(sd = meansd - mean) %>%
    { c(.$mean, .$sd) }
}

cor_height_weight <- 0.91
cor_height_bmi <- 0.92
cor_weight_bmi <- 0.93

age <- 10
age_month <- 3
sex <- 'female'

mu_vec <- data_height_age %>%
  filter(Age == age * 12 + age_month,
         sex == 'female') %>%
  {c(.$mean, .$sd)}

cor_m <- matrix(c(
  1                , cor_height_weight, cor_height_bmi,
  cor_height_weight, 1                , cor_weight_bmi,
  cor_height_bmi   , cor_weight_bmi   , 1
                  ), nrow = 3, byrow = TRUE)

# simulate 100 10 yo girls data

d <- mvrnorm(100, mu = c(10, 150, 20), Sigma = cor_m) %>%
  as.data.frame() %>%
  setNames(c("age", "height", "weight")) %>%
  mutate(
    bmi = weight / (height / 100)^2
  )






# Simulate data based on the corr. matrix & linreg equation in the paper
fun_generate_child_data <- function(N = 100, age = 13, sex = "female") {

  theta_height <- fun_lookup_height(age, sex)
  theta_bmi    <- fun_lookup_bmi(   age, sex)

  df_sim <- data.frame(
    age = age,
    sex = sex,
    id = 1:N,
    height = rnorm(N, theta_height[1], theta_height[2]),
    bmi = rnorm(N, theta_bmi[1], theta_bmi[2])
  ) %>%
    mutate(
      weight = bmi * (height / 100)^2,
      theta_mm = ifelse(sex == "male",
                        0 - 28.669 + 0.887 * age + 0.298 * weight + 0.255 * height,
                        0 - 16.264 + 0.182 * age + 0.302 * weight + 0.198 * height)
      )

  df_sim$mm <- rnorm(N, df_sim$theta_mm, 5.071)

  return(as.data.frame(df_sim))

}

ages <- rnorm(1000, 13, 2) %>% round() %>% pmin(18) %>% pmax(10)

age_tab <- ages %>% table() %>% as.data.frame()

for (i in 1:nrow(age_tab)) {

  n_per_sex <- ceiling( age_tab$Freq[i] / 2)

  m <-  fun_generate_child_data(age = age_tab[i,1] %>%
                                      as.character() %>%
                                      as.numeric(),
                                  N = n_per_sex,
                                sex = "male")

  f <-  fun_generate_child_data(age = age_tab[i,1] %>%
                                        as.character() %>%
                                        as.numeric(),
                                  N = n_per_sex,
                                sex = "female")
  if (i == 1) {
    d <- bind_rows(m, f)
  } else {
    d <- bind_rows(m,f,d)
  }
}

d$mm_perc <- {d$mm / d$weight} %>% pmin(1)
d$mm_htsq <- {d$mm / (d$height / 100)^2} #%>% pmin(1)
d$mm_bmi  <- {d$mm / d$bmi} #%>% pmin(1)

fig_cormat <-
  d %>%
    select(!(c("id","bmi","theta_mm"))) %>%
      GGally::ggpairs(mapping=aes(alpha=.01))

fig_cormat + theme_minimal()


#d %>% select(!(c("id","sex","bmi","theta_mm"))) %>% pairs
d %>% select(!(c("id","sex","bmi","theta_mm"))) %>% cor

