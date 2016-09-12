arabidopsis <- read.table('data-raw/arabidopsis.tsv', header = TRUE)
devtools::use_data(arabidopsis, overwrite = TRUE)
