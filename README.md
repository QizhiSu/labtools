
<!-- README.md is generated from README.Rmd. Please edit that file -->

# labtools

<!-- badges: start -->
<!-- badges: end -->

The goal of labtools is to to provide a set of tools to help facilitate
the handling of laboratory data. At the moment, there is only one
function to help read and align compound lists of different samples. It
is still under active development. In the future, more functions will be
incorporated based on the need of our lab.

## Installation

You can install the development version of labtools from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("QizhiSu/labtools")
```

## Read and combine Canvas data

To combine Canvas data, we first need to manually analyze GC x GC data
in Canvas, mark down peak of interest, and then export marked data in
.txt format. All .txt files should be put into a folder as the function
will read all .txt in the same folder and then combine them into a
single table by matching the chemical name. It will evaluate the
retention index and second dimensional retention time of a same compound
across them samples. If the differences are bigger than the defined
tolerance, it will tell you which samples have significantly different
RI or 2D RT, such that you can carefully check where the inconsistences
locate.

``` r
library(labtools)

data_path <- 'C:/data'
data <- read_canvas(data_path, 
                    ri_align_tolerance = 5,
                    rt_2d_tolerance = 0.05,
                    keep = 'area')

# to understand each argument, you can use the following code
?read_canvas
```
