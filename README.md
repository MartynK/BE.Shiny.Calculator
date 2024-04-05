
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Before we begin

- Do **not** use ’\_’ in the project’s name
- Do **not** copy the .git folder to the project folder!
- Do initialize your repo on GitHub
- Do clone your repo onto a folder on your local machine
- Do open a new ‘devel’ branch
- Do copy the template stuff into the new folder now on said branch
- Do update the ‘PROJECTNAME’ parts
- Do some actual work (optional)

# PROJECTNAME

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/PROJECTNAME)](https://CRAN.R-project.org/package=PROJECTNAME)
<!-- badges: end -->

The goal of PROJECTNAME is to … Have a home of my template for a new
project.

As a reminder for myself, at each version bump I need to update:

- README.Rmd file (duh)  
- devtools::build_readme()
- NEWS file (project name)  
- DESCRIPTION (depends etc.)  
- devtools::document() your project after defining new functions under
  /R  
- devtools::build_site()
- then devtools::install()  
- only then devtools::check()

Notes:

- usethis::create_project() is a great resource  
- usethis::create_tidy_package() is also great

## Installation

You can install the development version of PROJECTNAME from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("MartynK/MartysProjectTemplate")
```

You’d need to have R and RStudio installed on your computer for the full
experience. The *.html* outputs are (usually) available in the
*vignettes* and *docs* subfolders.

## Example

This is a basic example which shows you how to solve a common problem:

``` r
#library(PROJECTNAME)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
