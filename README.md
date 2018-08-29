# appinetwork
An R package for construction and Analysis of Protein-Protein Interactions (PPI) NETWORKS for complexes and biological processes


## Installation

To install this package from git, you will need the devtools package :

    install.packages("devtools")
    library("devtools")
    
Now we can install our package from Github using the following line:

    devtools::install_github("melinagallopin/appinetwork")

You need the version R 3.2.0+. If you do not have the latest version of R installed, it should work as long as you install the dependencies first with
the following block of code:

    install.packages( pkgs = c("R.methodsS3","rPython","stringr"), dependencies = TRUE)

## Examples

This package only implements a function `interface()` which ... 
