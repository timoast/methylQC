methylome <- methylQC::loadData(system.file("extdata", "methylome.CGmap.gz", package = "methylQC"))
devtools::use_data(methylome, overwrite = TRUE)
