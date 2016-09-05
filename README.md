# methylQC
Quality control for methylome data

## Install

Install from github

```R
# install.packages("devtools")
install_github("timoast/methylqc")
```

## Package overview

### Functions

* `loadData`: Loads a BSseeker2 CGmap file into memory. This just calls `data.table::fread` and gives each column the correct names.
* `cytosines`: Count the number of cytosines in each sequence context, for each chromosome, using a given fasta file.  
* `plotBrowser`: Generate a genome browser view of sequencing coverage for a given genomic region.  
* `coverageSurvival`: Calculate the diminishing percentage of all cytosines with increasing levels of sequencing depth.  
* `plotCoverage`: Plot the diminishing percentage of all cytosines with increasing levels of sequencing depth.  
* `methylomeStats`: Generate some summary statistics for the distribution of sequencing depth. This returns a dataframe with the quantiles, mean, and kurtosis (how "spiky" the coverage is).  
* `nonConversion`: Calculate the bisulfite non-conversion rates for each sequence context.  
* `strandBias`: Calculate bias in sequencing depth towards either strand for each C position.  


## Interactive web app  

This package is also implemented in an interactive app built using [shiny](http://shiny.rstudio.com/). To run the app locally, first make sure a few extra dependencies are installed:

```R
install.packages(c("shiny", "ploty"))
```

Then just run the app by running the `app.R` script. This will open a browser window where you can interact with the data.
