
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

The best practice for KLFDAPC is to use the below code instead of using the specific kernel function from kernelab (i.e., rbfdot kernel).

`````{r}

library(KLFDAPC)

# Open the GDS file
genofile <- SNPRelate::snpgdsOpen(snpgdsExampleFileName())
## obtaining pop code
pop_code <- read.gdsn(index.gdsn(genofile, "sample.annot/pop.group"))
pop_code <- read.gdsn(index.gdsn(genofile, path="sample.annot/pop.group"))
pop_code=factor(pop_code,levels=unique(pop_code))
## Doing PCA
pcadata <- SNPRelate::snpgdsPCA(genofile)
snpgdsClose(genofile)

## normalization function
normalize <- function(x) {
return ((x - min(x)) / (max(x) - min(x)))
}

pcanorm=apply(pcadata$eigenvect[,1:20], 2, normalize)

the Gaussian kernel
kmat <- kmatrixGauss(pcanorm,sigma=5)

klfdapc=KLFDA(kmat, y=pop_code, r=3, knn = 2)

``````


Welcome any [feedback](https://github.com/xinghuq/KLFDAPC/issues) and [pull request](https://github.com/xinghuq/KLFDAPC/pulls). 


## Citation

Qin. X. 2020. KLFDAPC: Kernel Local Fisher Discriminant Analysis of Principal Components (KLFDAPC) for large genomic data. R package version 0.2.0.

Qin, X., Chiang, C.W.K., and Gaggiotti, O.E. (2021). [Kernel Local Fisher Discriminant Analysis of Principal Components (KLFDAPC) significantly improves the accuracy of predicting geographic origin of individuals](https://www.biorxiv.org/content/10.1101/2021.05.15.444294v2.full). bioRxiv, 2021.2005.2015.444294.


