\name{dataset2}
\alias{dataset2}
\docType{data}
\title{
Real bivariate data set 
}
\description{
This example represents the measurements on the brinell hardness (X_1) and tensile strength (X_2) 
discussed by Wang and Chen (1998).
}
\usage{data(dataset2)}
\format{
  A data frame with 25 observations on the following 2 variables.
  \describe{
    \item{\code{X1}}{ brinell hardness }
    \item{\code{X2}}{ tensile strength}
  }
}

\references{
Wang F, Chen J (1998). "Capability index using principal components analysis." Quality Engineering, 11, 21-27.
}
\examples{
data("dataset2")
## maybe str(dataset2) ; plot(dataset2) ...
}
\keyword{datasets}
