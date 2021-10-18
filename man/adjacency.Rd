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
\item{D }{Dissimilarity matrix base on adjacency.}
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
data(datangkor)

## network stratigraphic data (Network)
network <- datangkor$stratigraphy

dissimilarity <- adjacency(network)
dissimilarity
}
