\name{datarcheo}
\alias{datarcheo}
\docType{data}
\title{
Data set of archeology
}
\description{
Latitude, longitude, Absolute difference between two dating and the name of the archeological site. Data concern  archaeological contexts datation in West Central France based on a large collection of mediaval pottery finds. Two original statistical models are developed to estimate context dates using pottery. The  absolute difference is calculated for each context. Data are based on a collective research about medieval ceramic directed by P. Husi (\href{http://citeres.univ-tours.fr/spip.php?article493}{\code{"La céramique médiévale dans la vallée de la Loire moyenne"}}) and from archeo-statistical modelling project \href{http://citeres.univ-tours.fr/spip.php?article506#CeraR}{\code{CeraR}} (L. Bellanger & Ph. Husi).
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
