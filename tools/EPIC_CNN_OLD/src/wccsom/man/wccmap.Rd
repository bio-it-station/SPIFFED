\name{wccmap}
\alias{wccmap}
\title{Map one powder pattern to a trained Kohonen map}
\description{
  Function calculates the agreement - as measured by WCC - of a powder
  pattern to all units in a network (SOM or XYF).
}
\usage{
wccmap(x, obj)
}
\arguments{
  \item{x}{Kohonen network, either from wccsom or wccxyf}
  \item{obj}{Vector: new powder pattern}
}
\value{
  Returns a vector of length equal to the number of units in the
  network, containing all WCC values, i.e. similarities of the new
  pattern to every unit.
}
\references{R. Wehrens, W.J. Melssen, L.M.C. Buydens and R. de Gelder.
Representing Structural Databases in a Self-Organising Map.
Acta Cryst. B61, 548-557, 2005.}
\author{Ron Wehrens}
\seealso{ \code{\link{wccsom}}, \code{\link{wccxyf}}, \code{\link{wccassign}}}
\examples{
data(cepha)
gr <- somgrid(3,3, "hexagonal")
set.seed(7)
x <- wccsom(cepha$patterns, grid=gr, trwidth=20, rlen=100)

wccs1 <- wccmap(x, cepha$patterns[1,])
par(mfrow=c(1,2))
plot(x, "property", property = wccs1, main="Unit similarities to object 1")
plot(x, "property", property = wccs1,
     main="Unit similarities to object 1", zlim=c(0.96, 1))
}
\keyword{classif}
