\name{cepha}
\alias{cepha}
\docType{data}
\title{Cephalosporin data set}
\description{
  X-ray powder patterns of 20 cephalosporin / antibiotic complexes.
}
\usage{data(cepha)}
\format{
  This yields a list with three components: the first component,
  '"patterns"', is a matrix of 20 rows and 425 variables, containing
  the powder patterns; the second component is '"class.names"', and
  gives information on the class of the crystal structure. The final
  component, "thetas", contains the 2theta values at which intensities
  have been measured.
}
\references{R. de Gelder, R. Wehrens, and J.A. Hageman.
J. Comput. Chem., 22(3), 273-289, 2001.}
\examples{
data(cepha)
plot(cepha$thetas, cepha$patterns[1,], type="l", xlab="2theta", ylab="Counts")
matplot(cepha$thetas, t(cepha$patterns), type="l",xlab="2theta", ylab="Counts",
        col=as.integer(factor(cepha$class.names))+1, lty=1)
}
\keyword{datasets}
