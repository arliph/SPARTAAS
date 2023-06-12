\name{CAdist}
\alias{CAdist}
\title{
Distance matrix based on correspondence analysis results}
\description{
Perform a correspondence analysis on a contingency table and then return the distance matrix of the coordinates (you can choose the number of axes to use to build the distance matrix with the nCP parameter).
}
\usage{
CAdist(df, nPC = NULL, graph = TRUE)
}
\arguments{
  \item{df}{
Data.frame, matrix or table with the data for the correspondence analysis.}
  \item{nPC}{
Number of principal components to be retained for the construction of the distance matrix. Must be between 1 and the minimum of ncol - 1 and nrow - 1. Could also be "max".}
  \item{graph}{
Logical parameter for plotting the Correspondence Analysis (Axis1, Axis2).}
}

\value{
\item{D }{The distance matrix}
}
\author{
A. COULON

L. BELLANGER

P. HUSI
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
library(SPARTAAS)
data(datangkor)

## contingency table
cont <- datangkor$contingency

distance <- CAdist(cont, nPC = "max")
distance

## run without printing the plot
distance <- CAdist(cont, nPC = "max", graph=FALSE)
}
