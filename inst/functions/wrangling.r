# Data wrangling

fil <- here::here("inst","extdata","rds.xls")


descriptor <-
  here::here("inst","extdata","description.xlsx") %>%
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

  datachunk <- fil |>
    readxl::read_xls(sheet = dataset) |> # or read_xlsx() as appropriate
    mutate( across( .cols = which( descriptor$trf == "factor"),
                    .fns = as.factor
    ),
    across( .cols = which( descriptor$trf == "numeric"),
            .fns = as.numeric # removing potential '?', 'NA', '.' etc.
    ),
    across( .cols = which( descriptor$trf == "date"),
            .fns = lubridate::as_datetime
    )) %>%
    .[,1:5] %>% # Some datasets have a "logPK" column which is superfluous
    `colnames<-`( descriptor$name_new) %>%
    labelled::`var_label<-`(   labs  ) %>%
    mutate(dataset = dataset)

  # binding the actual dataset to a 'master list'
  if ( dataset == sheets[1]) {
    data <- datachunk
  } else {
    data <- bind_rows(data, datachunk)
  }

}





