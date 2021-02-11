\name{CAdist}
\alias{CAdist}
\title{
Distance matrix based on correspondence analysis results}
\description{
Run an Correspondences analysis on a contingence table then return the matrix distance of the coordinates (you can choose the number of axes to use to build the distance matrix with the nCP parameter).
}
\usage{
CAdist(df, nPC = NULL, graph = TRUE)
}
\arguments{
  \item{df}{
Data.frame, matrix or table with the data for the correspondence analysis}
  \item{nPC}{
Number of principal components to be retained for the construction of the distance matrix. Must be between 1 and the minimum of ncol - 1 and nrow - 1. Could also be "max".}
  \item{graph}{
Logical parameter for plot the Correspondences Analysis (axis1, axis2)}
}

\value{
\item{dist }{The distance matrix}
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

## contingency table
cont <- data.frame(
  Cat10 = c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
  Cat20 = c(4,8,0,0,0,0,0,0,0,6,0,0,0,0,0,0,0),
  Cat30 = c(18,24,986,254,55,181,43,140,154,177,66,1,24,15,0,31,37),
  Cat40.50 = c(17,121,874,248,88,413,91,212,272,507,187,40,332,174,17,288,224),
  Cat60 = c(0,0,1,0,0,4,4,3,0,3,0,0,0,0,0,0,0),
  Cat70 = c(3,1,69,54,10,72,7,33,74,36,16,4,40,5,0,17,13),
  Cat80 = c(4,0,10,0,12,38,2,11,38,26,25,1,18,4,0,25,7),
  Cat100.101 = c(23,4,26,51,31,111,36,47,123,231,106,21,128,77,10,151,114),
  Cat102 = c(0,1,2,2,4,4,13,14,6,6,0,0,12,5,1,17,64),
  Cat110.111.113 = c(0,0,22,1,17,21,12,20,30,82,15,22,94,78,18,108,8),
  Cat120.121 = c(0,0,0,0,0,0,0,0,0,0,66,0,58,9,0,116,184),
  Cat122 = c(0,0,0,0,0,0,0,0,0,0,14,0,34,5,0,134,281),
  row.names = c("AI01","AI02","AI03","AI04","AO03","AI05","AO01","AI07","AI08",
  "AO02","AI06","AO04","APQR01","APQR02","AO05","APQR03","AI09")
)

distance <- CAdist(cont, nPC = "max")
distance

## run without printing the plot
distance <- CAdist(cont, nPC = "max", graph=FALSE)
}
