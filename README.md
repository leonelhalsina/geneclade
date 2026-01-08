# geneclade
Simulation platform for genetic diversity patterns in a clade radiation

### Install 

To install the R package geneclade, make sure you have devtools installed and then type from R (Rstudio):

```
library(devtools)
remotes::install_github("leonelhalsina/geneclade",ref="cleaning")
```
### Using geneclade

We have prepared a vignette (a sort of manual with chunks of code) that can
be called once you load the library:

```
library(geneclade)
devtools::install(build_vignettes = TRUE)
# or perhaps you need to build the vignette since package installation, so do:
remotes::install_github("leonelhalsina/geneclade",ref="cleaning",build_vignettes = TRUE)
browseVignettes("geneclade")
```
