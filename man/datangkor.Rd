\name{datangkor}
\alias{datangkor}
\docType{data}
\title{
Data set of archeology
}
\description{
The archaeological data come from excavations carried out at Angkor Thom (Cambodia), the capital of the Khmer Empire between the 9th and 15th centuries (Gaucher, 2004). The dataset consists of the pottery assemblages (quantities of different types of pottery shards contained in the sets - ..$contingency) and the stratigrpahy of the sets from 3 disconnected archaeological sites (..$stratigraphy).
}
\author{
Gaucher, J. (2004). Angkor Thom, une utopie réalisée ? Structuration de l’espace et modèle indien d’urbanisme dans le Cambodge ancien. Arts Asiatiques, Volume 59, pp. 58-86.
}
\usage{data("datangkor")}
\format{
  List of two objects with 17 observations.
  \describe{
\item{$contingency}{(data.frame) Contingency table of the quantities of 12 types of pottery shards in the 17 sets.}
\item{$stratigraphy}{(data.frame) Saves the stratigraphic network. First column corresponds to the nodes (sets) and the second to the edges by listing the nodes connected to it.}
  }
}

\examples{
  library(SPARTAAS)
  data(datangkor)
  str(datangkor)
  str(datangkor$contingency)
  str(datangkor$stratigraphy)
}
