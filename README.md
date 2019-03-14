# appinetwork
An R package for constructing and Analysing of Protein-Protein Interactions (PPI) NETWORKs for complexes and biological processes


## Installation

To install this package from git, you need the devtools package using the following command in R :

    install.packages("devtools")
    library("devtools")
    
Now you can install our package from Github using the following command in R :

    devtools::install_github("melinagallopin/appinetwork")

You need the version R 3.2.0+. If you do not have the latest version of R installed, it should work as long as you install the dependencies first with
the following commands in R:

    install.packages( pkgs = c("R.methodsS3","rPython","stringr"), dependencies = TRUE)

## Examples

This package implements a function `interface()` which ...  Put a link to the vignette, userguide ...

Currently, the final function `tfit` to cluster the network of second degree physical interactions does not work using the interface. We have to use the function  `tfit()` directly in command line. 

Try.