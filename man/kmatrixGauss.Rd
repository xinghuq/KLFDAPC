\name{kmatrixGauss}
\alias{kmatrixGauss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimating Gaussian Kernel matrix
%%  ~~function to do ... ~~
}
\description{This function estimates Gaussian kernel computation for klfda, which maps the original data space to non-linear and higher dimensions. See the deatils of kmatrixGauss from lfda.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
kmatrixGauss(x, sigma = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Input data matrix or dataframe
%%     ~~Describe \code{x} here~~
}
  \item{sigma}{ The Gaussian kernel parameter
%%     ~~Describe \code{sigma} here~~
}
}
\details{ Return a n*n matrix
%%  ~~ If necessary, more details than the description above ~~
}
\value{Return a n*n matrix
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Tang, Y., & Li, W. (2019). lfda: Local Fisher Discriminant Analysis inR. Journal of Open Source Software, 4(39), 1572.
%% ~put references to the literature/web site here ~
}

