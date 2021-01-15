\name{plot.wccsom}
\alias{plot.wccsom}
\alias{plot.wccchanges}
\alias{plot.wcccounts}
\alias{plot.wcccodes}
\alias{plot.wccpred}
\alias{plot.wccprop}
\alias{plot.wccmapping}
\alias{plot.wccquality}
\alias{plot.heatkey}
\title{Plot wccsom object}
\description{
  Plot self-organising map, obtained from function wccsom. Several types
  of plots are supported.
}
\usage{
\method{plot}{wccsom}(x, type = c("changes", "codes", "counts", "mapping",
                        "prediction", "property", "quality"),
            classif = NULL, labels = NULL, pchs = NULL, main = NULL,
            palette.name = heat.colors, ncolors, unit.colors,
            unit.bgcol = NULL, zlim = NULL, property = NULL,
            heatkey=TRUE, contin, \dots)
}
\arguments{
  \item{x}{SOM object.}
  \item{type}{Type of plot. (Wow!)}
  \item{classif}{classification object, as returned by \code{wccassign}.
    Only needed if \code{type} equals '"mapping"', '"quality"',
    and '"counts"'. Also a vector of class numbers may be given. If the
    network was trained with \code{keep.data} equal to TRUE, then the
    object already contains this information.}
  \item{labels}{(optional) labels to plot when \code{type} equals
    '"mapping"'.}
  \item{pchs}{(optional) plotting symbols to use when \code{type} equals
    '"mapping"'.}
  \item{main}{title of the plot.}
  \item{palette.name}{colors to use as unit background for '"codes"',
    '"counts"', '"prediction"', '"property"' and '"quality"' plotting
    types.}
  \item{ncolors}{number of colors to use in the palette.}
  \item{unit.colors}{explicit definition of the colors for the units in
    the '"codes"', '"counts"', '"prediction"', '"property"' and
    '"quality"' plotting types.}
  \item{unit.bgcol}{background color for units if no other color
    information is to be plotted. For plotting type is '"codes"', the
    default is '"transparent"'; in other cases the default is '"gray"'.}
  \item{zlim}{Optional range for color coding of unit backgrounds.}
  \item{property}{Values to use if the '"property"' plotting type. Can be used
    for colouring units in general.}
  \item{heatkey}{Whether or not to generate a heatkey at the left side
    of the plot in the '"quality"', '"counts"', and '"property"'
    plotting types.}
  \item{contin}{Whether the heatkey should show a range of values
    (TRUE) or a series of discrete values. The function tries to make a
    good guess; in case of strange-looking results it may pay to
    explicitly provide a value for this argument.}
  \item{\dots}{Other graphical parameters, e.g. colours of labels in the
    '"mapping"' plotting type.}
}
\details{Several different types of plots are supported:
  \describe{
    \item{'"changes"'}{Shows the mean change in similarity to the best
      matching codebook vector for each epoch. Since codebook vectors
      become more similar to the data that are mapped to them, the
      changes should always be positive. Upon convergence, the changes
      should be very small.}
    \item{'"codes"'}{Shows the codebook vectors.}
    \item{'"counts"'}{Shows the number of objects mapped to the
      individual units. Empty units are depicted in gray.}
    \item{'"mapping"'}{Show where a set of objects is mapped. It
      needs a '"labels"' argument: a string name for each object.}
    \item{'"prediction"'}{Shows predictions for units; if no '"property"'
      argument is given, the function will for supervised maps the
      predictions for every unit; for unsupervised maps, where this
      information is not available, it gives an error.}
    \item{'"property"'}{Plot a map with the units coloured according to
      a specific property. The standard application is to precompute the
      similarities of an object to all units in the map, and plot these
      with a colour key. Also other quantities may be
      used to colour the units: see the example of distances below.}
    \item{'"quality"'}{Shows the units coloured according to the mean
      agreement (WCC) of mapped objects to the unit vector. A colour key is
      plotted on the left. The variation in the WCCs of the mapped
      objects is indicated by the blue line: if it is pointing downwards it
      indicates low variation, if pointing upwards large variation.}
  }
}
\value{If \code{type} equals '"property"', the wcc values for all units
  are returned.}
\references{R. Wehrens, W.J. Melssen, L.M.C. Buydens and R. de Gelder.
Representing Structural Databases in a Self-Organising Map.
Acta Cryst. B61, 548-557, 2005.}
\author{Ron Wehrens}
\seealso{\code{\link{wccsom}}, \code{\link{wccxyf}},
  \code{\link{unit.distances}}, \code{\link{predict.wccsom}}}
\examples{
data(cepha)
gr <- somgrid(3, 3, "hexagonal")
set.seed(7)
x <- wccsom(cepha$patterns, grid = gr, trwidth = 20, rlen = 100)

par(mfrow = c(3,2))
plot(x, type = "codes", main = "Codebook vectors")
plot(x, type = "changes", main = "Convergence")
plot(x, type = "counts", main = "Mapping counts")
plot(x, type = "quality", main = "Mapping quality")
plot(x, type = "mapping", main = "Mapping", 
     labels = cepha$class.names, col = as.integer(factor(cepha$class.names)))
plot(x, type = "mapping", main = "Mapping",
     pchs = as.integer(factor(cepha$class.names)),
     col = as.integer(factor(cepha$class.names)))

par(mfrow=c(1,1))
obj1.wccs <- wccmap(x, cepha$patterns[1,])
plot(x, type = "property", property = obj1.wccs)
}
\keyword{classif}
