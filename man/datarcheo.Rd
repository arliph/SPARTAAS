\name{datarcheo}
\alias{datarcheo}
\docType{data}
\title{
Data set of archeology
}
\description{
Latitude, longitude, absolute difference between two dates and the name of the archaeological site. The data concern the dating of archaeological contexts in West-central France based on a large collection of medieval pottery finds. Two original statistical models are developed to estimate context dates using pottery. The  absolute difference is calculated for each context. The data are based on a collective research on medieval pottery directed by P. Husi (\href{https://citeres.univ-tours.fr/spip45d1.html?article493}{\code{"La céramique médiévale dans la vallée de la Loire moyenne"}}).
}
\author{
P. Husi IR CNRS, UMR CITERES-LAT, CNRS/Tours University, France :
}
\usage{data("datarcheo")}
\format{
  List of three objects with 240 observations on the following 4 variables.
  \describe{
\strong{$coord (data.frame):}
    \item{\code{longitude}}{a numeric vector: geographical coordinate}
    \item{\code{latitude}}{a numeric vector: geographical coordinate}
    \strong{$var (vector):}
    \item{\code{regionalized_var}}{a numeric vector: difference between two dating model}
    \strong{$label (vector):}
    \item{\code{noms}}{a character vector(factor): name of archeological site}
  }
}

\examples{
  library(SPARTAAS)
  data(datarcheo)
  str(datarcheo)
  head(datarcheo$coord)
  str(datarcheo$var)
  levels(datarcheo$label)
}
