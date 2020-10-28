\name{KLFDAPT}
\alias{KLFDAPT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Genome scan for detecting natural selection based on KLFDAPC
%%  ~~function to do ... ~~
}
\description{KLFDAPT implements genome scan based on KLFDAPC. The implementation is similar with PCA-based genome scan that reported by Patterson, N.,et al. 2006 and Luu, K., et al. 2017. The loci are weighted by the correlations between genetic variation and the genetic structure (first K reduced features).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
KLFDAPT(klfdapc, gdsfile, estimation="auto", snp.id = NULL, n.RD = NULL, num.thread = 1L, with.id = TRUE, outgds = NULL, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{klfdapc}{KLFDAPC matrix. The results returned from KLFDAPC analysis.
%%     ~~Describe \code{klfdapc} here~~
}
  \item{gdsfile}{ The genotype file, it is an object of class SNPGDSFileClass, a SNP GDS file
%%     ~~Describe \code{gdsfile} here~~
}
\item{estimation}{a character string specifying the robust estimator to be used. The choices are: "mcd" for the Fast MCD algorithm of Rousseeuw and Van Driessen, "weighted" for the Reweighted MCD, "donostah" for the Donoho-Stahel projection based estimator, "M" for the constrained M estimator provided by Rocke, "pairwiseQC" for the orthogonalized quadrant correlation pairwise estimator, and "pairwiseGK" for the Orthogonalized Gnanadesikan-Kettenring pairwise estimator. The default "auto" selects from "donostah", "mcd", and "pairwiseQC" with the goal of producing a good estimate in a reasonable amount of time.
%%     ~~Describe \code{verbose} here~~
}
  \item{snp.id}{a vector of snp id indicating the selected SNP id, if NULL, all SNPs are used
%%     ~~Describe \code{snp.id} here~~
}
  \item{n.RD}{a vector of integers indicating which reduced features to be used in KLFDAPC
%%     ~~Describe \code{eig.which} here~~
}
  \item{num.thread}{the number of (CPU) cores used; if NA, detect the number of cores automatically
%%     ~~Describe \code{num.thread} here~~
}
  \item{with.id}{if TRUE, the returned value with sample.id and sample.id
%%     ~~Describe \code{with.id} here~~
}
  \item{outgds}{NULL or a character of file name for exporting correlations to a GDS file, see details
%%     ~~Describe \code{outgds} here~~
}
  \item{verbose}{if TRUE, show information
%%     ~~Describe \code{verbose} here~~
}
}
\details{ For more details, please read our vignette (https://github.com/xinghuq/KLFDAPC).
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{cor}{The correlation coefficients between the genotypes and the population genetic structure}
 \item{q.values }{The p.values and q.values for the correlation coefficients}
%% ...
}
\references{
Patterson, N., Price, A. L., & Reich, D. (2006). Population structure and eigenanalysis. PLoS genet, 2(12), e190.

Luu, K., Bazin, E., & Blum, M. G. (2017). pcadapt: an R package to perform genome scans for selection based on principal component analysis. Molecular ecology resources, 17(1), 67-77.

A High-performance Computing Toolset for Relatedness and Principal Component Analysis of SNP Data. Xiuwen Zheng; David Levine; Jess Shen; Stephanie M. Gogarten; Cathy Laurie; Bruce S. Weir. Bioinformatics 2012; doi: 10.1093/bioinformatics/bts606.

R. A. Maronna and R. H. Zamar (2002) Robust estimates of location and dispersion of high-dimensional datasets. Technometrics 44 (4), 307–317.

Wang, J., Zamar, R., Marazzi, A., Yohai, V., Salibian-Barrera, M., Maronna, R., ... & Maechler, M. (2017). robust: Port of the S+“Robust Library”. R package version 0.4-18. Available online at: https://CRAN. R-project. org/package= robust.
%% ~put references to the literature/web site here ~
}
\author{
qinxinghu@gmail.com
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
###Hap_map data
### Doing KLFDAPC, default kernel=kernlab::polydot(degree = 1, scale = 1, offset = 1)
klfdapcfile <- KLFDAPC(snpgdsExampleFileName(), y=pop_code, snp.id=snpset.id,n.pc = 3, kernel=kernlab::polydot(degree = 1, scale = 1, offset = 1), r=5,num.thread=2)

# Make a data.frame
pop_str <- data.frame(sample.id = klfdapcfile$PCA$sample.id,
    pop = factor(pop_code)[match(klfdapcfile$PCA$sample.id, sample.id)],
    RD1 = klfdapcfile$KLFDAPC$Z[,1],    # the first reduced feature
    RD2 =  klfdapcfile$KLFDAPC$Z[,2],    # the second reduced feature
    stringsAsFactors = FALSE)
head(pop_str)

rownames(pop_str)=pop_str$sample.id
genofile <- snpgdsOpen(snpgdsExampleFileName())
chr <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
klfdapcmat=as.matrix(pop_str[,-c(1:2)])

### using the correlation to estimate the contribution of SNPs to the pop str (which is the same as PCA does )

klfdapcscanqvalue <- KLFDAPT(klfdapcmat, genofile, n.RD=1:2)

###Mmanhattan plot
 plot(-log10(klfdapcscanqvalue$q.values$q.values), ylim=c(0,4), xlab="", ylab="-log(q-value)",col=chr, pch="+")
   abline(h=2, col="blue")
 }


