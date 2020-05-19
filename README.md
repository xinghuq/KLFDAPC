# KLFDAPC
# Kernel local Fisher discriminant analysis of principal components (KLFDAPC) for large genomic data
## Install packages

library("devtools")

devtools::install_github("xinghuq/KLFDAPC")

## Deependencies
requireNamespace("SNPRelate")

if (!requireNamespace("BiocManager", quietly=TRUE))

  install.packages("BiocManager",repos = "http://cran.us.r-project.org")
  
if (!requireNamespace("SNPRelate", quietly=TRUE))

  BiocManager::install("SNPRelate")
  
 if (!requireNamespace("DA", quietly=TRUE))
 
  devtools::install_github("xinghuq/DA")

### example
f <- system.file('extdata',package='KLFDAPC')

infile <- file.path(f, "2019-nCoV_total.gds")

##Note the labels below is random

y1=rep(1,times=1736)
y2=rep(2,times=2000)

y=rbind(as.matrix(y1),as.matrix(y2))

y=as.factor(y)
##### using gaussan kernel
##### This will take longer than PCA, denpending on the number of samples and n.pcs. We will not show the results here. Users can test on their own clusters

virus_klfdapc=KLFDAPC(infile,y,kernel=kernlab::rbfdot(sigma = 0.5),r=3,snp.id=NULL, maf=0.05, missing.rate=0.05,n.pc=10,tol=1e-30, num.thread=2,metric = "plain",prior = NULL)

showfile.gds(closeall=TRUE)

##### Plot the reduced features
