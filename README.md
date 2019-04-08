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

Currently, the final function `tfit` to cluster the network of second degree physical interactions does not work using the interface. We have to use the function  `tfit()` directly in command line. 


## Formatting database and PPI network analysis 

### Construct an ID correspondence file (theasaurus)

The internal function __theasaurus_window__ displays a dialog window that enables to select the name of the organism of interest and the type and the name of the Uniprot file, from which an ID correspondence file is built.

This function return a file of protein IDs correspondances with the isoforms number stored in a directory having the organism name. This file is used to format the Biogrid database and to build the Protein Interaction Network.

Organisms Names proposed by the dialog window are: Caenorhabditis elegans, Drosophila melanogaster, Escherichia coli, Homo sapiens, Mus musculus, Rattus norvegicus, Saccharomyces cerevisiae. If the organism of interest is not in this list, one can enter its name by selecting "Other" in the drop-down menu.
The uniprot file is a file (in .txt format) with all proteins informations of an organism. Two types of files can be downloaded from the website  [http://www.uniprot.org/](http://www.uniprot.org/): the file containing the proteome of the organism and the one containing all its ORFs. The file's type has to be chosen in the dialog window: "proteome"" or "all" (by default) .

### Format databases

#### Intact database formatting

This function displays a dialog window that enables to select the name of the organism of interest and the name of the original Intact file. Then it parses and formates this file in a unique format common to all databases. 

This function returns a file with parsed interactions data in a unique format common to all databases. Original Intact datafiles come from the website, [http://]().

Organisms Names proposed by the dialog window are: Caenorhabditis elegans, Drosophila melanogaster, Escherichia coli, Homo sapiens, Mus musculus, Rattus norvegicus, Saccharomyces cerevisiae. Their organism ID are : Caenorhabditis elegans : 6239, Drosophila melanogaster : 7227, Escherichia coli : 562, Homo sapiens : 9606, Mus musculus : 10090, Rattus norvegicus : 10116, Saccharomyces cerevisiae : 559292. This ID is usefull for the original Intact file name. If the organism of interest is not in this list, one can enter its name by selecting "Other" in the drop-down menu.

#### IrefIndex database formatting & parsing

Function irefindex_window

This function displays a dialog window that enables to select the name of the organism of interest, the maximum number of proteins in a complex and the name of the original IrefIndex file. Then it parses and formates this file in a unique format common to all databases. 


Original Irefindex datafiles come from the website: [http://irefindex.org/download/irefindex/data/archive/release_14.0/psi_mitab/MITAB2.6/](http://irefindex.org/download/irefindex/data/archive/release_14.0/psi_mitab/MITAB2.6/)

Organisms Names proposed by the dialog window are: Caenorhabditis elegans, Drosophila melanogaster, Escherichia coli, Homo sapiens, Mus musculus, Rattus norvegicus, Saccharomyces cerevisiae. If the organism of interest is not in this list, one can enter its name by selecting "Other" in the drop-down menu.

#### Biogrid database formatting

The internal function __biogrid__ formates a biogrid file to use it in build_network function.

This function returns the BIOGRID file in a unique format common to all databases. 

The Biogrid file is the file BIOGRID-ORGANISM*.tab2 downloaded from the website [https://downloads.thebiogrid.org/BioGRID/Release-Archive/](https://downloads.thebiogrid.org/BioGRID/Release-Archive/).
The uniprot file is a file (in .txt format) with all proteins informations of an organism. Two types of files can be downloaded from the website [http://www.uniprot.org/](http://www.uniprot.org/): the file containing the proteome of the organism and the one containing all its ORFs. The file's type has to be chosen in the dialog window: "proteome"" or "all" (by default). Organisms ID are the following : (Caenorhabditis elegans : 6239, Drosophila melanogaster : 7227, Escherichia coli : 562, Homo sapiens : 9606, Mus musculus : 10090, Rattus norvegicus : 10116, Saccharomyces cerevisiae : 559292).

Organisms Names proposed by the dialog window are: Caenorhabditis elegans, Drosophila melanogaster, Escherichia coli, Homo sapiens, Mus musculus, Rattus norvegicus, Saccharomyces cerevisiae. If the organism of interest is not in this list, one can enter its name by selecting "Other" in the drop-down menu.


### Modeling complexes assembly intermediaries

Uses a network and an input list of proteins in complexe to modelise intermediary subcomplexes of assembly.
To perform lmodeling complexes assembly intermediaries, the internal function needs a network and an input list. The function creates a .jpeg file representing a tree of subcomplexes modelised.

### Clustering with TFit method

Network clustering by TFit to find all interactants with a complexe or in a biological process. This clustering only applied on second degree networks.
