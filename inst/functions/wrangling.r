# Data wrangling

fun_read_excel <- function(fil, descr., type = "xlsx") {

  descriptor <-
    descr. %>%
    file.path() %>%
    readxl::read_excel(skip = 0)

  # make a list to be used as labels for 'labelled::'
  labs <-
    # take a 1xN matrix
    matrix( descriptor$description,
            ncol = length(descriptor$description)
    ) %>%
    as.list() %>%
    # add column names as 'names'
    `names<-`(descriptor$name_new)

  sheets <- readxl::excel_sheets(fil)

  for (dataset in sheets) {

    if (type == "xlsx") {
      datachunk <- fil |>
        readxl::read_xlsx(sheet = dataset)
    } else {
      datachunk <- fil |>
        readxl::read_xls(sheet = dataset)
    }

    datachunk <- datachunk |> # or read_xls() as appropriate
      # handle duplicate colnames
      janitor::clean_names() |>
      mutate( across( .cols = which( descriptor$trf == "factor"),
                      .fns = as.factor
      ),
      across( .cols = which( descriptor$trf == "numeric"),
              .fns = as.numeric # removing potential '?', 'NA', '.' etc.
      ),
      across( .cols = which( descriptor$trf == "date"),
              .fns = lubridate::as_datetime
      )) %>%
      `colnames<-`( descriptor$name_new) %>%
      # select only columnswhose name isnt NA due to duplicates or whatnot
      select( descriptor$name_new[!is.na(descriptor$name_new)])%>%
      labelled::`var_label<-`(   labs  ) %>%
      mutate(dataset = dataset)

    # binding the actual dataset to a 'master list'
    if ( dataset == sheets[1]) {
      data <- datachunk
    } else {
      data <- bind_rows(data, datachunk)
    }

  }


  # deselect the "not_relevant" columns
  try( {
    data <- data %>%
      select( -not_relevant)
  }, silent = TRUE)

}

fil <- here::here("inst","extdata","who",
                  "hfa-boys-perc-who2007-exp (1).xlsx"
)

descr <- here::here("inst","extdata","description.xlsx")

data_boys <- fun_read_excel(fil, descr) %>% mutate(sex = "male")


fil2 <- here::here("inst","extdata","who",
                  "hfa-girls-perc-who2007-exp (1).xlsx"
)

data_girls <- fun_read_excel(fil2, descr) %>% mutate(sex = "female")

data_height_age <- bind_rows( data_boys, data_girls)


###########
# BMI


fil3 <- here::here("inst","extdata","who",
                  "bmi-boys-perc-who2007-exp.xlsx"
)

descr <- here::here("inst","extdata","description_bmi.xlsx")


data_boys_bmi <- fun_read_excel(fil3, descr) %>% mutate(sex = "male")

fil4 <- here::here("inst","extdata","who",
                  "bmi-girls-perc-who2007-exp.xlsx"
)

data_girls_bmi <- fun_read_excel(fil4, descr) %>% mutate(sex = "female")

data_bmi <- bind_rows( data_boys_bmi, data_girls_bmi)

# cleanup
rm( data_boys, data_boys_bmi, data_girls, data_girls_bmi)


###########
# Weights

fil5 <- here::here("inst","extdata","who",
                   "hfa-boys-perc-who2007-exp.xlsx"
)


descr <- here::here("inst","extdata","description_weights.xlsx")


data_boys_wght <- fun_read_excel(fil5, descr) %>% mutate(sex = "male")


fil6 <- here::here("inst","extdata","who",
                   "hfa-girls-perc-who2007-exp.xlsx"
)


descr <- here::here("inst","extdata","description_weights.xlsx")


data_girls_wght <- fun_read_excel(fil6, descr) %>% mutate(sex = "female")

data_wght <- bind_rows( data_boys_wght, data_girls_wght)

# cleanup
rm( data_boys_wght, data_girls_wght)



