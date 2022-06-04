
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

## GDS format 

GDS is portable across platforms with hierarchical structure to store multiple scalable array-oriented data sets with metadata information. Details can be found [here](http://bioconductor.org/packages/release/bioc/vignettes/gdsfmt/inst/doc/gdsfmt.html).

You can convert PLINK/Oxford files to GDS using SNPRelate (Zheng et al. 2012).

[snpgdsBED2GDS	Conversion from PLINK BED to GDS](https://rdrr.io/bioc/SNPRelate/man/snpgdsBED2GDS.html).
[snpgdsGEN2GDS	Conversion from Oxford GEN format to GDS](https://rdrr.io/bioc/SNPRelate/man/snpgdsGEN2GDS.html).
[snpgdsPED2GDS	Conversion from PLINK PED to GDS](https://rdrr.io/bioc/SNPRelate/man/snpgdsPED2GDS.html)

Below, is an example on how to convert PLINK bed files to GDS, 

``````
library("SNPRelate")
# PLINK BED files
bed.fn ="/file_location/target.bed"
fam.fn = "/file_location/target.fam"
bim.fn = "/file_location/target.bim"

# convert
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "HapMap.gds")

# open
genofile <- snpgdsOpen("HapMap.gds")
genofile

# close
snpgdsClose(genofile)
``````

## Implementation of KLFDAPC

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
## Plotting the genetic structure

``````{r}
plot(klfdapc$Z[,1], klfdapc$Z[,2], col=as.integer(pop_code), xlab="RD 1", ylab="RD 2")
legend("topright", legend=levels(pop_code), pch="o", col=1:nlevels(pop_code))


``````

Welcome any [feedback](https://github.com/xinghuq/KLFDAPC/issues) and [pull request](https://github.com/xinghuq/KLFDAPC/pulls). 


## Citation

Qin. X. 2020. KLFDAPC: Kernel Local Fisher Discriminant Analysis of Principal Components (KLFDAPC) for large genomic data. R package version 0.2.0.

Qin, X., Chiang, C.W.K., and Gaggiotti, O.E. (2022). [KLFDAPC: A Supervised Machine Learning Approach for Spatial Genetic Structure Analysis](https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbac202/6596986). Briefings in Bioinformatics, bbac202, https://doi.org/10.1093/bib/bbac202.


