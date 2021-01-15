\name{wccassign}
\alias{wccassign}
\title{Assign patterns to nodes in a SOM network by WCC value}
\description{
  KNN assignment of patterns to units in a Kohonen map, with maximal WCC
  as the criterion.
}
\usage{
wccassign(x, data)
}
\arguments{
  \item{x}{Trained Kohonen map}
  \item{data}{Data matrix}
}
\value{
  Returns a list with components:
  \item{classif}{Unit numbers to which rows in the data matrix are assigned}
  \item{wccs}{wcc value of rows in the data matrix and the units to
    which they are assigned.}
}
\author{Ron Wehrens}
\seealso{\code{\link{wccsom}}}
\examples{
data(cepha)
gr <- somgrid(3,3, "hexagonal")
set.seed(7)
x <- wccsom(cepha$patterns, grid = gr, trwidth = 20, rlen = 50,
            FineTune = FALSE, keep.data = FALSE)
sombins <- wccassign(x, cepha$patterns)
}
\keyword{classif}
