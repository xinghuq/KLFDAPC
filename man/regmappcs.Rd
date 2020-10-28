\name{regmappcs}
\docType{data}
\alias{regmappcs}
\title{Principal components of RegMap data on A. thaliana}
\description{
  Data on 1,307 \emph{Arabidopsis thaliana} accessions collected around
  the world. Data include geographic co-ordinates and results of
  applying principal components analysis (PCA) to genotypes
  at 214,051 single nucleotide polymorphisms (SNPs). This data is obtained from Novembre, et al.2008.See PCAviz package.
}
\usage{data("regmappcs")}
\format{
  \code{regmappcs} is a data frame with 1,307 rows and 21 columns, in which
  each row represents an \emph{A. thaliana} accession included
  in the RegMap panel. The last 11 columns (labeled \code{"PC1"} through
  \code{"PC10"}) store the projection of the samples onto the first
  10 PCs.
}
\source{
  \emph{Arabidopsis thaliana} SNP genotypes were obtained from the
  website of Joy Bergelson's lab at the University of Chicago,
  \url{http://bergelson.uchicago.edu/?page_id=790}. We used the data from Novembre, et al.2018 PCAviz package.
}
\details{
  PCs were computed using \code{rpca} from the \pkg{rsvd} package after
  first converting the genotypes to a 1,307 x 214,051 binary matrix.
}
\references{
  Horton, M. W., Hancock, A. M., Huang, Y. S., Toomajian, C., Atwell,
  S., Auton, A., \emph{et al} (2012). Genome-wide patterns of genetic
  variation in worldwide Arabidopsis thaliana accessions from the RegMap
  panel. \emph{Nature Genetics} \bold{44}, 212–216.
  Novembre, J., R. Williams, H. Pourreza, Y. Wang, and P. Carbonetto.(2018) "PCAviz: Visualizing principal components analys." R Package.

}
\examples{

# View the first few data samples.
data(regmappcs)
head(regmappcs)

\dontrun{
# See the vignette for more examples.
vignette("regmappcs")
}
}
\keyword{datasets}
