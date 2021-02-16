\name{mapclust_app}
\alias{mapclust_app}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Shiny application for ClustCutMap method
}
\description{
This function performs a divisive hierarchical clustering on a regionalized variable using spatial patches algorithm  (Woillez et al. 2007). It is a hierarchical top-down clustering with geographical constraint. It is possible to cut the tree by clicking on the dendrogramm at the level you want. Results include description of the clusters and graphics. When cutting the dendrogram you can look at the two plots (WSSPlot and AveSilPlot) which shows the relative good quality of the partitions. The first one shows the Within Sum of Square (WSS) for each partition and you can use the Elbow approach in order to select one partition. The second plot show the average silhouette width. This index is between -1 and 1. The closer it is to 1 the more the partition is good. See the fonction help page of \code{silhouette{cluster}}
}
\usage{
mapclust_app()
}

\details{
Run the function without argument to open the application. You can import your own data in the Uploading tab. The calculation starts in the "ClustCutMap" tab. Select data, press run and choose the partition.
}

\author{
A. Coulon
}

\examples{
#open application
#library(SPARTAAS)
#if(interactive()){
  #mapclust_app()
#}
}

