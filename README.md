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
    
You can copy paste the following R commands all at once :


    install.packages("devtools")
    library("devtools")
    install.packages( pkgs = c("R.methodsS3","rPython","stringr"), dependencies = TRUE)
    devtools::install_github("melinagallopin/appinetwork")

## Examples

This package implements a function `interface()` which enables to format PPI databses and performs two mains analysis : model complexes assembly intermediaries and cluster PPI networks. The step by step procedure is described here : [userguide](https://github.com/melinagallopin/data/blob/master/userguide.pdf)



## For more information : read the  [userguide](https://github.com/melinagallopin/data/blob/master/userguide.pdf)

