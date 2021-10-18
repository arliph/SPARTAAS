\name{datacancer}
\alias{datacancer}
\docType{data}
\title{
Data set of cancerology.
}
\description{
Longitude, latitude, number of thyroid cancers. The data concern two departments of France (Loire Atlantique and Vendee) between 1998 and 2012. For reasons of confidentiality they are simulated data.
}
\author{
M. Karakachoff (IR CHU - l'institut du Thorax INSERM UMR 1087 - CNRS UMR 6291) Nantes, France

F. Molinie (resp. Loire-Atlantique-Vendee cancer registry -
\href{https://www.registre-cancers-44-85.fr/}{\code{registre-des-cancers}}) France
}

\usage{data("datacancer")}
\format{
  List of two objects with 359 observations on the following 3 variables.

  \strong{$coord (data.frame):}
  \describe{
    \item{\code{longitude}}{a numeric vector: geographical coordinate}
    \item{\code{latitude}}{a numeric vector: geographical coordinate}
  }
  \strong{$var (vector):}
  \describe{
    \item{\code{var}}{a numeric vector: number of thyroid cancers (simulated)}
  }
}


\examples{
  library(SPARTAAS)
  data(datacancer)
  str(datacancer)
  head(datacancer$coord)
  str(datacancer$var)
}

