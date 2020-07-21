# Infinity Flow: High-throughput single-cell quantification of 100s of proteins using conventional flow cytometry and machine learning

This repository contains the software package ***infinityFlow*** designed to analyze massively-parallel cytometry experiments, such as BioLegend's *LEGENDScreen* or's BD *Cell surface screening panels*. 

You can learn more about the Infinity Flow approach using the following sources. If you use this software in your research, please cite **Infinity Flow: High-throughput single-cell quantification of 100s of proteins using conventional flow cytometry and machine learning, Becht, Tolstrup, Dutertre, Ginhoux, Newell, Gottardo and Headley, 2020**.

1. [Our preprint](https://www.biorxiv.org/content/10.1101/2020.06.17.152926v1)
1. The peer-reviewed article - [pending]
1. [An imaged digest of the paper](https://twitter.com/EtienneBecht/status/1274039148781826049)

The package is currently pending publication on Bioconductor. In the meantime you can install it from github:

```
if(!require(devtools)){
	install.packages("devtools")
}
if(!require(infinityFlow)){
	library(devtools)
	install_github("ebecht/infinityFlow")
}
```

After installing the package, you can learn about how to use it by running the following command in an R session

```
browseVignettes("infinityFlow")
```