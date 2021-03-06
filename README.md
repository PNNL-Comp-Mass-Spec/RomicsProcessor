# RomicsProcessor

RomicsProcessor is an R package that can be used to analyze omics data.
The package provides a structured R object to store the data, allowing for reproducible
data analysis. The package also supports creating analytical pipelines from 
previously processed objects and applying these pipeline to other objects.
This allows for rapid development and reuse of bioinformatics methods.


## Installation

To install the package in R please first make sure devtools is installed

```
install.packages("devtools")

```

When devtools is installed, run the following command to install RomicsProcessor and its dependencies

```
devtools::install_github(“PNNL-Comp-Mass-Spec/RomicsProcessor”)

```

Alternately, download the 
[built package](https://github.com/PNNL-Comp-Mass-Spec/RomicsProcessor/blob/master/RomicsProcessor_1.0.0.tar.gz)
and install the package manually.

We recommend the installation of the Bioconductor package 'sva' if batch corrections are needed and of 'pmartR'.

## Example of use

The folder /Example contain an 
[R markdown file](https://github.com/PNNL-Comp-Mass-Spec/RomicsProcessor/blob/master/Example/Bacillus_cereus_media_experiment.Rmd)
.

This example consist in an analysis of the proteome of B. cereus grown in different media. 
The Raw LC-MS/MS data is publicly available, on [MassIVE](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=6d0ca42ca79244a49d66a80fd741ba28), FTP at the following address (ftp://massive.ucsd.edu/MSV000085696/).

It was analyzed using MaxQuant set with the 
[attached parameters](https://github.com/PNNL-Comp-Mass-Spec/RomicsProcessor/blob/master/Example/parameters.txt)
. 

The MaxQuant output 
[proteinGroups.txt](https://github.com/PNNL-Comp-Mass-Spec/RomicsProcessor/blob/master/Example/proteinGroups.txt) 
Can be used as data and the 
[metadata](https://github.com/PNNL-Comp-Mass-Spec/RomicsProcessor/blob/master/Example/metadata.csv) 
is also provided.


The different examples dataset are accessible directly inside the package 
[data folder](https://github.com/PNNL-Comp-Mass-Spec/RomicsProcessor/tree/master/data).

## Cite the code

To cite the package please use the following DOI:
[![DOI](https://zenodo.org/badge/206400976.svg)](https://zenodo.org/badge/latestdoi/206400976)

## Contacts

Written by @GeremyClair for the Department of Energy (PNNL, Richland, WA) \
E-mail: geremy.clair@pnnl.gov or proteomics@pnnl.gov \
Website: https://omics.pnl.gov/ or https://panomics.pnnl.gov/

## License

RomicsProcessor is licensed under the 2-Clause BSD License; 
you may not use this file except in compliance with the License.  You may obtain 
a copy of the License at https://opensource.org/licenses/BSD-2-Clause

Copyright 2019 Battelle Memorial Institute
