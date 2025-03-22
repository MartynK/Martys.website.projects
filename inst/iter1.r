# Load necessary libraries
library(MASS)  # For mvrnorm function

source(here::here("inst", "functions", "load_stuff.r"))
source(here::here("inst", "functions", "wrangling.r"))


fun_lookup_height <- function(age = 10, sex. = "female",
                              data_age. = data_height_age) {
  age_months <- runif(1, age * 12, age * 12 + 11) %>% round()

  data_age. %>%
    filter(Age == age_months
           ,sex == sex.
           )  %>%
    {c(.$mean, .$sd)}
}

fun_lookup_bmi <- function( age = 10, sex. = "female",
                            data. = data_bmi) {
  age_months <- runif(1, age * 12, age * 12 + 11) %>% round()

  data. %>%
    filter(Age == age_months
           ,sex == sex.
    ) %>%
    mutate( sd = meansd - mean ) %>%
    {c(.$mean, .$sd)}

}

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

