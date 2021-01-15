\name{predict.wccsom}
\alias{predict.wccsom}
\title{Predict properties from self-organising maps}
\description{
  Function to predict property values for every unit in a supervised or
  unsupervised SOM. These, in turn, are used to provide predictions for
  individual objects.
}
\usage{
\method{predict}{wccsom}(object, newdata,
               trainX, trainY, unit.predictions, \dots)
}
\arguments{
  \item{object}{Trained SOM.}
  \item{newdata}{If new objects are supplied (in the form of a matrix),
    they are mapped to the SOM; predictions for the new data are the
    predicted values associated with the units to which they are mapped.
    In order to calculate these the maps should be either supervised, or
    the training data should be stored in the map, or these data should
    be provided through arguments '"trainX"' and '"trainY"', or finally,
    the unit predictions can be explicitly given ('"unit.predictions"').}
  \item{trainX}{Training data, only used when they have not been stored
    in the trained map.}
  \item{trainY}{Dependent values for the training data.}
  \item{unit.predictions}{Alternatively, one can provide predictions for
    every unit.}
  \item{\dots}{Not used.}
}
\details{
  For supervised SOMs, predictions per unit are available after
  training. For unsupervised SOMs, these predictions can be obtained
  from the average values of the properties of training set objects
  mapping to specific units. New objects that are mapped to the SOM will
  receive the predicted value of the unit to which they are mapped.
}
\value{
  The function returns a list with components
  \item{unit.predictions}{Property predictions per unit of the map.}
  \item{predictions}{Property predictions for the new data.}
}
\references{ ~put references to the literature/web site here ~ }
\author{R. Wehrens}
\seealso{\code{\link{wccsom}}, \code{\link{wccxyf}},
  \code{\link{plot.wccsom}}}
\examples{
\dontrun{
data(degelder)
gr <- somgrid(5, 5, "hexagonal")
set.seed(7)
x <- wccxyf(degelder$patterns, degelder$properties[,"cell.vol"],
            grid=gr, trwidth=20, rlen=100)
plot(x, "predict")

predicted.volumes <- predict(x)
plot(degelder$properties[,"cell.vol"], predicted.volumes$predictions,
     xlab="Cell volume", ylab="Predicted cell volume")
abline(0,1, col="gray")
}
}
\keyword{classif}
