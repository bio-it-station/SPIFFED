\name{wccxyf}
\alias{wccxyf}
\title{Supervised mapping of spectra with self-organising maps}
\description{
  Supervised self-organising maps for mapping high-dimensional spectra
  or patterns to 2D; instead of Euclidean distance, the weighted cross
  correlation (WCC) similarity measure is used. Modelled after the SOM
  function in package 'class'. wccxyf takes 'continous' patterns,
  i.e. datapoints are equidistant.

  At this point, no facilities are implemented for growing networks or
  k-means-like fine-tuning of the maps, such as in function \code{wccsom}.
}
\usage{
wccxyf(data, Y, grid=somgrid(), rlen = 100, alpha = c(0.05, 0.01),
       radius = quantile(nhbrdist, 0.67), xweight = 0.5, trwidth = 20,
       toroidal = FALSE, keep.data = TRUE)
}
\arguments{
  \item{data}{Spectra or patterns to be mapped: a matrix, with each row
    representing a compound.}
  \item{Y}{Property for each pattern, either a numerical vector or
    matrix, or a class matrix. In the latter case, the Tanimoto distance
    is used for Y; in all other cases (also for combinations of
    numerical and class properties) the Euclidean distance is used.} 
  \item{grid}{A grid for the representatives: see 'somgrid'.}
  \item{rlen}{the number of times the complete data set will be
    presented to the network.}
  \item{alpha}{a vector of two numbers indicating the amount of
    change. Default is to decline linearly from 0.05 to 0.01
    over \code{rlen} updates.}
  \item{radius}{the initial radius of the neighbourhood to be used for
    each update: the decrease is exponential over \code{rlen} updates in
    such a way that after one-third of the updates only the winning unit
    is updated. The default is to start with a
    value that covers 2/3 of all units.}
  \item{xweight}{weight of X matrix in determining the distances of
    objects to units.}
  \item{trwidth}{width of the triangle function used in the WCC
    measure, given in the number of data points.}
  \item{toroidal}{if TRUE, then the edges of the map are joined. Note
    that in a toroidal hexagonal map, the number of rows must be even.}
  \item{keep.data}{store training data and their mapping in the network.}
}
\value{
  an object of class '"wccsom"' with components
  \item{grid}{the grid, an object of class '"somgrid"'.}
  \item{changes}{vector of mean average deviations from code vectors}
  \item{codes}{a matrix of code vectors.}
  \item{trwdth}{the triangle width used for the WCC measure}
  \item{acors}{autocorrelations of the code vectors.}
  \item{toroidal}{setting of parameter 'toroidal'.}
  \item{FineTune}{setting of parameter 'FineTune'.}
  \item{unit.classif}{mapping of training data: a vector of unit
    numbers. Only if \code{keep.data} equals TRUE.}
  \item{wccs}{WCC values of all training data, compared to the best
    matching codebook vector. Only if \code{keep.data} equals TRUE.}
  \item{data.acors}{WAC values for training data. Only if
    \code{keep.data} equals TRUE.}
}
\references{FIXME: this page is a copy of wccsom, should be edited further}
\author{Ron Wehrens}

\seealso{\code{\link[class]{SOM}}, \code{\link{plot.wccsom}},
  \code{\link{wccsom}}, \code{\link{wcc}}}
\examples{
\dontrun{
data(degelder)
gr <- somgrid(5, 5, "hexagonal")
set.seed(7)
x <- wccxyf(degelder$patterns, degelder$properties[,"cell.vol"],
            grid=gr, trwidth=20, rlen=100)
}
}
\keyword{classif}
