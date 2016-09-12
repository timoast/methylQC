methylome <- methylQC::loadData("data-raw/methylome.CGmap.gz")
devtools::use_data(methylome, overwrite = TRUE)
