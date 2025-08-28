# geneclade
Simulation platform for genetic diversity patterns in a clade radiation

### Install 

To install the R package geneclade, make sure you have devtools installed and then type from R (Rstudio):

```
library(devtools)
.libPaths("/opt/software/uoa/spack-sw/linux-rhel8-x86_64/gcc-12.1.0/r-4.1.3-gixu37cbx7dcnwzya2lu7jimi4qd2gl6/rlib/R/library")

library(devtools)
.libPaths("/uoa/scratch/users/s06lh9/R/x86_64-pc-linux-gnu-library/4.1")

remotes::install_github("leonelhalsina/geneclade",ref="cleaning")
```
### Using geneclade

We have prepared a vignette (a sort of manual with chunks of code) that can
be called once you load the library:

```
library(geneclade)
devtools::install(build_vignettes = TRUE)
# or perhaps you need to build the vignette since package installation, so do:
remotes::remotes::install_github("leonelhalsina/geneclade",ref="cleaning",build_vignettes = TRUE)
browseVignettes("geneclade")
```
