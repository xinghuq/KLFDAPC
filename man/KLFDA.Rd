\name{KLFDA}
\alias{KLFDA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Kernel local Fisher discriminant analysis
%%  ~~function to do ... ~~
}
\description{This function performs Kernel local Fisher discriminant analysis. The function provided here allows users to carry out the KLFDA using other input datasets. We used the gaussan matrix as example. Users can input the
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
KLFDA(kdata, y, r, sigma = 0.5, metric = c("weighted", "orthonormalized", "plain"),tol=1e-5,knn = 6, reg = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kdata}{The input dataset (kernel matrix). The input data can be a genotype matrix, dataframe, species occurence matrix, or principal components. The dataset have to convert to a kernel matrix before feed into this function.
%%     ~~Describe \code{pcdata} here~~
}
  \item{y}{ The group lables
%%     ~~Describe \code{y} here~~
}
  \item{r}{ Number of reduced features
%%     ~~Describe \code{r} here~~
}
  \item{sigma}{ The sigma used to estimate gaussan matrix.
%%     ~~Describe \code{sigma} here~~
}
  \item{metric}{Type of metric in the embedding space (default: 'weighted') 'weighted' - weighted eigenvectors 'orthonormalized' - orthonormalized 'plain' - raw eigenvectors
%%     ~~Describe \code{metric} here~~
}
  \item{knn}{The number of nearest neighbours
%%     ~~Describe \code{knn} here~~
}
  \item{reg}{The regularization parameter
%%     ~~Describe \code{reg} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 \item{Z}{The reduced features}
\item{Tr}{The transformation matrix}
%% ...
}
\references{
Tang, Y., & Li, W. (2019). lfda: Local Fisher Discriminant Analysis inR. Journal of Open Source Software, 4(39), 1572.
%% ~put references to the literature/web site here ~
}

\seealso{
KLFDAPC
}
\examples{
kmat <- kmatrixGauss(iris[, -5],sigma=1)
zklfda=KLFDA(kmat, iris[, 5], r=3,tol=1e-5, metric = "plain")
print(zklfda$Z)
}

