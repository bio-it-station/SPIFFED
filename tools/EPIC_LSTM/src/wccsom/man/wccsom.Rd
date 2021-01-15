\name{wccsom}
\alias{wccsom}
\title{Mapping spectra with self-organising maps}
\description{
  Self-organising maps for mapping high-dimensional spectra or patterns
  to 2D; instead of Euclidean distance, the weighted cross correlation
  (WCC) similarity measure is used. Modelled after the SOM function in
  package 'class'. wccsom takes 'continous' patterns, i.e. datapoints
  are equidistant.
}
\usage{
wccsom(data, grid=somgrid(), rlen = 100, alpha = c(0.05, 0.01),
       radius = quantile(nhbrdist, 0.7), init, nhbrdist, trwidth = 20,
       toroidal = FALSE, FineTune = TRUE, keep.data = TRUE)
}
\arguments{
  \item{data}{Spectra or patterns to be mapped: a matrix, with each row
    representing a compound.}
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
  \item{init}{the initial representatives, represented as a matrix.  If
    missing, chosen (without replacement) randomly from 'data'.}
  \item{nhbrdist}{optionally, the distance matrix for the units.}
  \item{trwidth}{width of the triangle function used in the WCC
    measure, given in the number of data points.}
  \item{toroidal}{if TRUE, then the edges of the map are joined. Note
    that in a toroidal hexagonal map, the number of rows must be even.}
  \item{FineTune}{apply kmeans for fine-tuning the codebook vectors.}
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
\references{R. Wehrens, W.J. Melssen, L.M.C. Buydens and R. de Gelder.
Representing Structural Databases in a Self-Organising Map.
Acta Cryst. B61, 548-557, 2005.}
\author{Ron Wehrens}

\seealso{\code{\link[class]{SOM}}, \code{\link{plot.wccsom}},
  \code{\link{wccxyf}}, \code{\link{wcc}}}
\examples{
data(cepha)
gr <- somgrid(3,3, "hexagonal")
set.seed(7)
x <- wccsom(cepha$patterns, grid=gr, trwidth=20, rlen=100)
}
\keyword{classif}
