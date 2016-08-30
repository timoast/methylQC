# methylQC
Quality control for methylome data

## Install

Install from github

```
# install.packages("devtools")
install_github("timoast/methylqc")
```

## Package overview

### Functions

* `loadData`: Loads a BSseeker2 CGmap file into memory. This just calls `readr::read_tsv` and gives each column the correct names.  
* `plotBrowser`: Generate a genome browser view of sequencing coverage for a given genomic region.  
* `plotSurvival`: Plot the diminishing percentage of all cytosines with increasing levels of sequencing depth.  
* `methylomeStats`: Generate some summary statistics for the distribution of sequencing depth. This returns a dataframe with the quantiles, mean, and kurtosis (how "spiky" the coverage is).  


## Interactive web app  

This package is also implemented in an interactive app built using [shiny](http://www.shinyapps.io/). To run the app locally, first make sure a few extra dependencies are installed:

```
install.packages(c("shiny", "ploty", "shinythemes"))
```

Then just run the app by executing the `app.R` script. This will open a browser window where you can interact with the data.