\name{wcc}
\alias{wcc}
\alias{wac}
\alias{wacmat}
\title{Agreement between patterns including peak shifts}
\description{
  Weighted cross correlation and autocorrelation, as described in De
  Gelder et al. (2001), for assessing similarities in spectra-like data
  containing peak shifts. Euclidean distances are useless in this
  situation.
}
\usage{
wcc(pattern1, pattern2, trwdth, wghts, acors)
wac(pattern1, trwdth, wghts)
wacmat(patterns, trwdth, wghts, do.transpose = TRUE)
}
\arguments{
  \item{pattern1}{Pattern.}
  \item{pattern2}{Another pattern.}
  \item{patterns}{Pattern matrix: rows correspond with
    patterns.}
  \item{trwdth}{Triangle width, given in the number of data points.}
  \item{wghts}{Optional weights vector, will be calculated from triangle
  width if necessary. Sometimes it is more efficient to pre-calculate it
  and give it as an argument.}
  \item{acors}{Autocorrelation, also optional to speed up calculations.}
  \item{do.transpose}{Internally, columns should correspond with
    patterns, so normally one should leave this value to its default:
    TRUE. If a matrix is already in the correct format, one can avoid
    unnecessary double transpose operations.}
}
\value{
  Function \code{wcc} returns the WCC value, a similarity value between
  0 and 1. Functions \code{wac} and \code{wacmat} return weighted
  autocorrelations for one pattern and a matrix of patterns,
  respectively.
}
\references{R. de Gelder, R. Wehrens, and J.A. Hageman. 
  A generalized expression for the similarity spectra: application to
  powder diffraction pattern classification. 
  J. Comput. Chem., 22(3), 273-289, 2001.}
\author{Ron Wehrens}

\seealso{\code{\link{wccsom}}, \code{\link{wccxyf}}}
\examples{
data(cepha)
wac(cepha$patterns[1,], 20)
wacmat(t(cepha$patterns), 20)
wcc(cepha$patterns[1,], cepha$patterns[2,], 20)
}
\keyword{classif}
