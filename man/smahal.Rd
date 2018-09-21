\name{smahal}
\alias{smahal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates robust, rank based, Mahalanobis distance matrix.
}
\description{
This is a variation of a common function in the literature of matching.  It takes
as input a numeric data matrix of \code{n} rows and outputs a \code{n} by \code{n}
matrix of pairwise \code{distance} of points.  The distances satisfy triangle 
inequality.
}
\usage{
smahal(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
	A numeric data matrix.
}
}
\details{
A call to this function will create a distance matrix which can be used
as an input of the main function \code{\link{makeblocks}}.  For more specialized
use see the examples of \code{\link{makeblocks}}.  

The robust Mahalanobis distance uses ranks of covariates rather than the covariates 
themselves, but the variances of the ranks are not adjusted for ties, so ties do not make a 
variable more important.

The returned matrix is the (entry wise) squared root of the rank based Mahalanobis distance matrix.

Any positive linear combination of outputs of multiple calls to this function
gives a valid distance matrix.
}
\value{
A \code{n} by \code{n} matrix of pairwise distances of the units.
}
\references{
Rosenbaum, P. R. (2010). Design of Observational Studies. New York: Springer.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
# See examples of 'makeblocks'.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distance }% use one of  RShowDoc("KEYWORDS")

