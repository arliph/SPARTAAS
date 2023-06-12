\name{datacerardat}
\alias{datacerardat}
\docType{data}
\title{
Data set of archaeology.
}
\description{
Data sets on Tours pottery.
Contains three objects:
- Pottery sherd contingency table (data.frame 280 sets and 391 technical groups)
- The date (if known) of the sets (vector)
- Column index for supplementary  sets (vector)
}
\usage{data("datacerardat")}
\format{
  A list with 3 objects.
  \describe{
    \item{\code{df}}{Ceramic contingency table. a integer data.frame}
    \item{\code{date}}{The date (if known) of the sets. a numeric or NA vector}
    \item{\code{col.sup}}{Column index for supplementary  sets. a numeric vector}
  }
}
\examples{
data(datacerardat)
datacerardat
## maybe str(datacerardat);
}
\keyword{datasets}
