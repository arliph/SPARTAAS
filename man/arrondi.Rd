\name{arrondi}
\alias{arrondi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Returns the rounded value.
}
\description{
Always returns the upper value if the next digit is 5.
}
\usage{
arrondi(x, acc)
}
\arguments{
  \item{x}{
The number to round.
}
  \item{acc}{
Accuracy (number of digits). A negative number of digits means rounding to the power of ten, e.g. arrondi(x, digits = -2) will round to the nearest hundred.}
}
\value{
  \item{res }{Value or vector of values rounded.}
}

\author{
A. COULON

L. BELLANGER

P. HUSI
}

\examples{
  library(SPARTAAS)

  x1 <- c(15,25,35,45,55)
  round(x1,-1)
  arrondi(x1,-1)

  x2 <- c(-15,-25,-35,-45,-55)
  round(x2,-1)
  arrondi(x2, -1)

  x3 <- 1.125
  round(x3,2)
  arrondi(x3, 2)

  x4 <- seq(-0.55,0.55,0.1)
  data.frame(
    val = x4,
    round = round(x4,1),
    arrondi = arrondi(x4, 1),
    equal = (arrondi(x4, 1) == round(x4,1))
  )

}

