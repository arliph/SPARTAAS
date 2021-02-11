\name{adjacency}
\alias{adjacency}
\title{
Dissimilarity matrix base on connectivity information.
}
\description{
From the data of a network, we build a contiguity matrix. Based on this matrix, we generate a dissmilarity matrix.
The matrix contains only 0 or 1, 1 if there is no relationship and 0 if there is a relationship.
The network object is a two-column data frame. The first column contains the elements of the network and the second column contains a list of all other elements related to it. The list is a character string consisting of the names of the elements separated by commas (see example).
}
\usage{
adjacency(network)
}
\arguments{
  \item{network}{
data frame with 2 columns. The first one contains all the elements (nodes) and the second one a string with all the elements related to it (links).
}
}
\value{
\item{D }{Normalized connection dissimilarity matrix}
}
\author{
A. COULON

L. BELLANGER

P. HUSI
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--  or do  help(data=index)  for the standard data sets.
library(SPARTAAS)

## network stratigraphic data (Network)
network <- data.frame(
  nodes = c("AI09","AI08","AI07","AI06","AI05","AI04","AI03",
  "AI02","AI01","AO05","AO04","AO03","AO02","AO01","APQR03","APQR02","APQR01"),
  edges = c("AI08,AI06","AI07","AI04","AI05","AI01","AI03","AI02","","","AO04","AO03",
  "AO02,AO01","","","APQR02","APQR01","")
)

dissimilarity <- adjacency(network)
dissimilarity
}
