\name{summary.wccsom}
\alias{summary.wccsom}
\title{Summarizing wccsom objects}
\description{
Summary for objects of class '"wccsom"'
}
\usage{
\method{summary}{wccsom}(object, type = c("unit", "object", "smoothness", "quality"),
               nr, labels, data = object$data,
               classif = object$unit.classif,
               wccs = object$wccs, properties = NULL, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Trained Kohonen map.}
  \item{type}{One of "unit", "object", "smoothness" or "quality".}
  \item{nr}{Number of the unit or object that is to be summarised
    (only needed for \code{type} "unit" and "object").}
  \item{labels}{Labels for all objects (not needed for \code{type}
    equaling "smoothness" and "quality"). Defaults to integer numbers.}
  \item{data}{Data matrix (not needed for \code{type} equaling
    "smoothness"). Default is to take the training
    objects (if available).}
  \item{classif}{Mapping of all objects (not needed for
    \code{type} equaling "smoothness"). Default is to take the mapping
    of the training objects (if available). Alternatively, a mapping
    object (output of \code{wccassign} may be given for this argument.}
  \item{wccs}{WCC values of all mapped objects (not needed for
    \code{type} equaling "smoothness"). Default is to take the training
    objects (if available).}
  \item{properties}{Other properties of the crystals that should be
    shown in the summary, as e.g. reduced cell parameters.}
  \item{\dots}{Further arguments. Currently ignored.}
}
\details{
  Several types of summary are calculated and printed to the screen.
  If \code{type} equals "unit", a summary is given of all objects mapped
  to that unit. In case \code{type} equals "object", the summary
  indicates the unit to which the object is mapped and shows other
  objects mapped to it. For a \code{type} of "smoothness", the function
  returns the ratio of WCC values between neighbouring and
  non-neighbouring units. Values larger than 1 show that neighbouring
  units are "more alike" than non-neighbouring units. Type "quality"
  does something similar: it compares all data objects and takes the
  ratio of WCCs between objects mapped to the same unit and objects
  mapped to different units.
}
\author{Ron Wehrens}
\seealso{\code{\link{wccsom}}}
\examples{
data(cepha)
gr <- somgrid(3,3, "hexagonal")
set.seed(7)
x <- wccsom(cepha$patterns, grid=gr, trwidth=20, rlen=100)

summary(x, "unit", nr=1, labels=paste(cepha$class.names, 1:20, sep=""))
summary(x, "object", nr=1, labels=paste(cepha$class.names, 1:20, sep=""))
summary(x, type = "smoothness")
summary(x, "quality")
}
\keyword{classif}
