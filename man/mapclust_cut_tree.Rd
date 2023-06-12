\name{mapclust_cut_tree}
\alias{mapclust_cut_tree}
\encoding{utf-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to cut the dendrogram for a new height (distance limit) or a new number of clusters and map the new partition
}
\description{
The function returns the new map, one dendrogram with the cut line, the silhouette width and the new vector cluster. You must call \link[SPARTAAS]{mapclust} beforehand to get a \code{mapclust_cl} object that can be used by mapclust_cut_tree.

}
\usage{
mapclust_cut_tree(classification, nb_grp = NA, dlim = NA)
}
\arguments{
  \item{classification}{
The object return by mapclust  \code{Class: mapclust_cl}
}
  \item{nb_grp}{
The number of groups you want to have on the partition. Must be an integer. (don't use dlim in this case)
}
  \item{dlim}{
The value of dlim at which you wish to cut the dendrogram. You can enter any value (numerical) and the function will select the nearest lower dlim with the same partition. (do not use nb_grp in this case).
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
\item{X}{The x-coordinate data you have used}
\item{Y}{The y-coordinate data you have used}
\item{var}{The regionalised variable data you have used}
\item{label}{The label vector you have used}
\item{density}{The estimated density based on var. Equal to var if you used a unidimensional density variable}
\item{cluster}{The cluster vector of the selected partition}
\item{Plabel}{Logical parameter to activate or not the printing of labels on the dendrogram}
\item{fullhist}{The compositional cluster for each observation}
\item{hist}{The compositional cluster without duplicates (corresponds to the split on the dendogram)}
\item{dlim}{The vector of the different limit distances}
\item{cutdlim}{The select dlim for the cut of the current partition}
\item{DiMatrix}{The weighted euclidean distance matrix}
\item{silhouetteData}{The silhouette data of the selected partition}
\item{AveSilData}{The average silhouette value for each partition}
\item{Moran}{The Moran index for each group for each partition}
\item{lonlat}{Logical parameter, whether your coordinates are in latitude-longitude format or not}
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

