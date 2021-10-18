\name{dataceram}
\alias{dataceram}
\docType{data}
\title{
Data set of archeology
}
\description{
This important dataset come from the Collective Research Project (CRP) on medieval and modern ceramics in the Middle Loire Basin. This project is a long-term one, since it began in 1996 and has already been the subject of two books on the subject (Husi dir. 2003 and 2013) as well as an online logicist publication (Husi dir. 2022).
}
\author{
Husi P. (dir.) – La céramique médiévale et moderne du bassin de la Loire moyenne, chrono-typologie et transformation des aires culturelles dans la longue durée (6e —19e s.), XXe Supplément à la Revue Archéologique du Centre de la France, FERACF, Tours, 2022.
}
\usage{data("dataceram")}
\format{
  List of three objects with 226 observations.
  \describe{
\item{$contingency}{(data.frame) Contingency table of the quantities of 183 types of pottery shards in the 226 sets.}
\item{$timerange}{(data.frame) The first column corresponds to the identifier (sets), the second and the third to the lower and upper bound of temporal range of the observations.}
\item{$geographic_area}{Vector of the geographical area of each observation.}
  }
}

\examples{
  library(SPARTAAS)
  data(dataceram)
  str(dataceram)
  str(dataceram$contingency)
  head(dataceram$timerange)
  head(dataceram$geographic_area)
}
