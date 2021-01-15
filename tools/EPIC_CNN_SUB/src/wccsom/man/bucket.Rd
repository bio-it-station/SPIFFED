\name{bucket}
\alias{bucket}
\alias{debucket}
\title{Variable averaging (bucketing) for data matrices}
\description{
Function \code{bucket} decreases the size (i.e., the number of columns)
of a data matrix by averaging variables. Function \code{debucket}
achieves the reverse by linear interpolation.
}
\usage{
bucket(x, factor)
debucket(x, nout)
}
\arguments{
  \item{x}{Data matrix: each variable corresponds with a column.}
  \item{factor}{Bucket factor: this number of variables will be
    averaged.}
  \item{nout}{Required number of variables after debucketing.}
}
\value{
  Returns a data matrix of the new dimensions.
}
\author{Ron Wehrens}

\examples{
data(cepha)
gr <- somgrid(3,3, "hexagonal")

X <- bucket(cepha$patterns, 4)
system.time(x <- wccsom(X, grid=gr, trwidth=5,
                        rlen=500, FineTune=FALSE))
## this is much faster than the following code:
\dontrun{
set.seed(7)
system.time(x <- wccsom(cepha$patterns, grid=gr, trwidth=20,
                        rlen=500, FineTune=FALSE))
}
}
\keyword{classif}
