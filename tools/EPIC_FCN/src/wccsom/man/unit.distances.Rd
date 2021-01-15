\name{unit.distances}
\alias{unit.distances}
\title{Function to calculate distances between units in a SOM}
\description{
  Function calculates Euclidean distances between units in a SOM; if
  argument '"toroidal"' is TRUE, the edges of the map are considered to
  be joined so that the overal shape of the map is a torus. The
  distances are calculated correspondingly.
}
\usage{
unit.distances(grid, toroidal)
}
\arguments{
  \item{grid}{A \code{somgrid} object.}
  \item{toroidal}{For toroidal maps, equal to TRUE. Default is FALSE.}
}
\value{Returns a distance matrix.}
\author{Ron Wehrens}
\seealso{\code{\link{wccsom}}, \code{\link{wccxyf}}}
\examples{
gr <- somgrid(3, 3, "hexagonal")
x <- list(grid = gr)
class(x) <- "wccsom" 

par(mfrow = c(1,2))
unit.dists <- unit.distances(gr, toroidal = FALSE)
plot(x, type = "property", property = unit.dists[1,],
     main = "Distances to unit 1", zlim = c(0,2.75), contin = TRUE)
unit.dists <- unit.distances(gr, toroidal = TRUE)
plot(x, type = "property", property = unit.dists[1,],
     main = "Toroidal distances to unit 1", zlim = c(0,2.75), contin = TRUE)
}
\keyword{classif}
