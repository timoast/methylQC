[![Build Status](https://travis-ci.org/timoast/methylQC.svg?branch=master)](https://travis-ci.org/timoast/methylQC)


# methylQC
Quality control for methylome data

## Install

Install from github

```R
# install.packages("devtools")
devtools::install_github("timoast/methylqc")
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

### Data

* `methylome`: A small test dataset including sequencing data for the lambda genome and the first ~1 million cytosines of Arabidopsis chromosome 1.  
* `arabidopsis`: Cytosine information for the Arabidopsis genome.

## Interactive web app  

This package is also implemented in an interactive app built using [shiny](http://shiny.rstudio.com/), available here:  
* [code](https://github.com/timoast/methylQC)  
* [demo](https://timoast.shinyapps.io/shiny/)
