\name{KLFDAPC}
\alias{KLFDAPC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Kernel local Fisher discriminant analysis of principal components (KLFDAPC)
%%  ~~function to do ... ~~
}
\description{ This function performs Kernel local Fisher discriminant analysis of principal components (KLFDAPC) for genomic data
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
KLFDAPC(infile, y, n.pc,sample.id=NULL, snp.id=NULL,
                 autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
                 algorithm=c("exact", "randomized"),
                 eigen.cnt=ifelse(identical(algorithm, "randomized"), 16L, 32L),
                 num.thread=1L, bayesian=FALSE, need.genmat=FALSE,
                 genmat.only=FALSE, eigen.method=c("DSPEVX", "DSPEV"),
                 aux.dim=eigen.cnt*2L, iter.num=10L, verbose=TRUE,
                 kernel=kernlab::polydot(degree = 1, scale = 1, offset = 1), r=3,
                 tol=1e-30,prior=NULL, CV=FALSE,usekernel = TRUE,
                 fL = 0.5,metric = c('weighted', 'orthonormalized', 'plain'),
                 knn = 6, reg = 0.001, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genofile}{ The genomic files. This package uses GDS file to perform a high performance computation. GDS formatted data is designed for efficient random access to large data sets.The GDS file can be converted from various formats such as PLINK data format (BED) and vcf files via SNPRelate package.
%%     ~~Describe \code{genofile} here~~
}
  \item{y}{ The group labels
%%     ~~Describe \code{y} here~~
}
  \item{n.pc}{ Number of PCs to proceed the KLFDA.
%%     ~~Describe \code{n.pc} here~~
}
  \item{snp.id}{The sample IDs that the users want to select to perform the analysis, if NULL or missing, the default is to use all samples.
%%     ~~Describe \code{snp.id} here~~
}
  \item{num.thread}{the number of (CPU) cores used; if NA, detect the number of cores automatically
%%     ~~Describe \code{num.thread} here~~
}
  \item{kernel}{ The kernel techniques used in KLFDA
%%     ~~Describe \code{kernel} here~~
}
  \item{r}{The number of reduced features
%%     ~~Describe \code{r} here~~
}
  \item{tol}{ The tolerance to decide if a matrix is singular; it will reject variables and linear combinations of unit-variance variables whose variance is less than tol^2.
%%     ~~Describe \code{tol} here~~
}
  \item{prior}{ The prior of the groups
%%     ~~Describe \code{prior} here~~
}
  \item{CV}{If true, returns results (classes and posterior probabilities) for leave-one-out cross-validation. Note that if the prior is estimated, the proportions in the whole dataset are used.
%%     ~~Describe \code{CV} here~~
}
  \item{usekernel}{Whether to use kernel classifier.whether to use kernel classifier, if TRUE, pass to Naive Bayes classifier.
%%     ~~Describe \code{usekernel} here~~
}
  \item{fL}{If usekernel is TRUE, pass to the kernel function. See the kernel function from kernlab.
%%     ~~Describe \code{fL} here~~
}
  \item{metric}{The type of metric in the embedding space (default: 'weighted') 'weighted' - weighted eigenvectors 'orthonormalized' - orthonormalized 'plain' - raw eigenvectors
%%     ~~Describe \code{metric} here~~
}
  \item{knn}{The number of nearest neighbours
%%     ~~Describe \code{knn} here~~
}
  \item{reg}{The regularization parameter
%%     ~~Describe \code{reg} here~~
}
  \item{\dots}{Other parameters pass to klfda and snpgdsPCA.
%%     ~~Describe \code{\dots} here~~
}
}
\details{ The kernel local discriminant analysis becomes complex once  the kernel tricks are intruduced into the function. We also employ the kernel classifier to make inference based on the non-linear features.
%%  ~~ If necessary, more details than the description above ~~
}
\value{ The results return two parts, the results of PCA and the results of KLFDAPC including the classified classes and the posterior possibility of each class using different classifier.
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{KLFDAPC}{The results of KLFDAPC including the classified classes and the posterior possibility of each class using different classifier}
\item{PCA}{PCA results}
%% ...
}
\references{

Sugiyama, M (2007).Dimensionality reduction of multimodal labeled data by local Fisher discriminant analysis. Journal of Machine Learning Research, vol.8, 1027-1061.

Sugiyama, M (2006). Local Fisher discriminant analysis for supervised dimensionality reduction. In W. W. Cohen and A. Moore (Eds.), Proceedings of 23rd International Conference on Machine Learning (ICML2006), 905-912.

Original Matlab Implementation: http://www.ms.k.u-tokyo.ac.jp/software.html#LFDA

Tang, Y., & Li, W. (2019). lfda: Local Fisher Discriminant Analysis inR. Journal of Open Source Software, 4(39), 1572.

Moore, A. W. (2004). Naive Bayes Classifiers. In School of Computer Science. Carnegie Mellon University.

Pierre Enel (2020). Kernel Fisher Discriminant Analysis (https://www.github.com/p-enel/MatlabKFDA), GitHub. Retrieved March 30, 2020.

Karatzoglou, A., Smola, A., Hornik, K., & Zeileis, A. (2004). kernlab-an S4 package for kernel methods in R. Journal of statistical software, 11(9), 1-20.

Zheng, X., Levine, D., Shen, J., Gogarten, S. M., Laurie, C., & Weir, B. S. (2012). A high-performance computing toolset for relatedness and principal component analysis of SNP data. Bioinformatics, 28(24), 3326-3328.
}
\author{qinxinghu@gmail.com
%%  ~~who you are~~
}

\examples{
# example file
f <- system.file('extdata',package='KLFDAPC')
infile <- file.path(f, "2019-nCoV_total.gds")
##Note the labels below is random
y1=rep(1,times=1736)
y2=rep(2,times=2000)
y=rbind(as.matrix(y1),as.matrix(y2))
y=as.factor(y)
### using gaussan kernel
### This will take longer than PCA, denpending on the number of samples and n.pcs.
#We will not show the results here. Users can test on their own clusters
###virus_klfdapc=KLFDAPC(infile,y,kernel=kernlab::rbfdot(sigma = 0.5),r=3,snp.id=NULL, maf=0.05, missing.rate=0.05,n.pc=10,tol=1e-30, num.thread=2,metric = "plain",prior = NULL)
#showfile.gds(closeall=TRUE)

}
