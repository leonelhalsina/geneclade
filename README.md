# geneclade
Simulation platform for genetic diversity patterns in a clade radiation

### Install 

To install the R package LEMAD, make sure you have devtools installed and then type from R (Rstudio):

```
library(devtools)
remotes::install_github("leonelhalsina/geneclade")

library(withr)
library(devtools)
new_dir <- tempdir(new)
dir.create(new_dir) 
with_libpaths(new_dir,install_github("leonelhalsina/geneclade"))

```
