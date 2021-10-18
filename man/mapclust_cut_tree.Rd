\name{mapclust_cut_tree}
\alias{mapclust_cut_tree}
\encoding{utf-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to cut the dendrogram for a new height (distance limit) or a new number of clusters and map the new partition
}
\description{
The function returns the new map, one dendrogram with the cutting line, the silhouette width and the new vector cluster. You must run \link[SPARTAAS]{mapclust} beforehand to have a \code{mapclust_cl} object that can be used by mapclust_cut_tree.

}
\usage{
mapclust_cut_tree(classification, nb_grp = NA, dlim = NA)
}
\arguments{
  \item{classification}{
The object return by mapclust  \code{Class: mapclust_cl}
}
  \item{nb_grp}{
The number of groups you want on the partition. Must be integer. (don't use dlim in this case)
}
  \item{dlim}{
The value of dlim where you want to cut the dendrogram. You can enter wathever you want (numeric value) and the function select for you the nearest lower dlim with the same partition. (don't use nb_grp in this case)
}
}

\value{
the function returns a list.

\strong{Plot:}
\item{dendrogram}{The global dendrogram (hclust object)}
\item{dendrogram_ggplot}{The global dendrogram (ggplot2 object)}
\item{cuttree}{The cut dendrogram}
\item{map}{The map of the selected partition}
\item{AveSilPlot}{The average silhouette width plot (for each partition)}
\item{WSSPlot}{The Within Sum of Square plot (for each partition)}
\item{silhouette}{The silhouette plot of the selected partition}
\strong{Value:}
\item{X}{The x-coordinates data you used}
\item{Y}{The y-coordinates data you used}
\item{var}{The regionalized variable.s data you used}
\item{label}{The label vector you used}
\item{density}{The estimate density based on var. Equal to var if you used a unidimensionnal density variable}
\item{cluster}{The vector of cluster of the selected partition}
\item{Plabel}{Logical parameter for activate or not the print of labels on the dendrogram}
\item{fullhist}{The composition cluster for each observation}
\item{hist}{The composition cluster without duplicates (matches to leaf of the dendrogram)}
\item{dlim}{The vector of the different limit distances}
\item{cutdlim}{The select dlim for the cut of the current partition}
\item{DiMatrix}{The matrix of Weighted Euclidiean distances}
\item{silhouetteData}{The silhouette data of the selected partition}
\item{AveSilData}{The average silhouette value for each partition}
\item{Moran}{The Moran index for each groups for each partitions}
\item{lonlat}{Logical parameter if your coordinates are in longitude latitude format or not.}
}

\author{
A. COULON
L. BELLANGER
P. HUSI
}

\examples{
\donttest{
## loading data
library(SPARTAAS)
data(datarcheo)

##First you need to run the mapclust function.
#This function allow you to obtain one partition
# object <- mapclust( coord = ..., var = ..., label = ...)
OldClassif <- mapclust(datarcheo$coord, datarcheo$var, datarcheo$label, n = 4)

##In order to cut the dendrogram for another dlim
NewCut <- mapclust_cut_tree(classification = OldClassif, dlim = 0.37)
##In order to cut the dendrogram for another number of cluster
NewCut2 <- mapclust_cut_tree(classification = OldClassif, nb_grp = 4)

#See evaluation using Silhouette width by running:
NewCut$silhouette
#If the plot is empty try to increase the height of the window (full screen).

#See summary of the data by running:
summary(NewCut$silhouetteData)
}
}

