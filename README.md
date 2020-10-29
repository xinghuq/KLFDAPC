
## Kernel Local Fisher Discriminant Analysis of Principal Components (KLFDAPC) for large genomic data



### Install packages

Install the most recent version of the KLFDAPC package using devtools:
`````{r}
library("devtools")

devtools::install_github("xinghuq/KLFDAPC")
``````
Alternatively, you can install from the source files, run the following commands in the shell:

```{shell}
R CMD build KLFDAPC
R CMD check --as-cran KLFDAPC_0.1.0.tar.gz
R CMD INSTALL KLFDAPC_0.1.0.tar.gz
```


### Dependencies

Before install or during installation, make sure the below dependences are installed.
``````{r}
requireNamespace("SNPRelate")

if (!requireNamespace("BiocManager", quietly=TRUE))

  install.packages("BiocManager",repos = "http://cran.us.r-project.org")
  
if (!requireNamespace("SNPRelate", quietly=TRUE))

  BiocManager::install("SNPRelate")
  
 if (!requireNamespace("DA", quietly=TRUE))
 
  devtools::install_github("xinghuq/DA")
``````

### Vignettes and tutorials

``````{r}

vignette("Population_structure_of_Covid")

vignette("Population_structure_of_RegMap")

vignette("Genome_scan_KLFDAPC")

``````



Welcome any [feedback](https://github.com/xinghuq/KLFDAPC/issues) and [pull request](https://github.com/xinghuq/KLFDAPC/pulls). 


## Citation

Qin. X. 2020. KLFDAPC: Kernel Local Fisher Discriminant Analysis of Principal Components (KLFDAPC) for large genomic data. R package version 0.2.0.

Qin, X., Gaggiotti, OE., Chiang, C. 2020. Kernel Local Fisher Discriminant Analysis of Principal Components (KLFDAPC) significantly improves the accuracy of predicting geographic origin of individuals. In submission.
