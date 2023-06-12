\name{mapclust}
\alias{mapclust}
\encoding{utf-8}
\title{
Divise hierarchical Clustering using Spatialpatches algorithm.
}

\description{

This function performs a divisive hierarchical clustering on a regionalised variable using the spatial patches algorithm (Woillez et al. 2007; Woillez, Rivoirard and Petitgas 2009). It is a top-down hierarchical clustering with a geographical constraint. It is possible to cut the tree by clicking on the dendrogram at the desired level. The results include a description of the clusters and graphics. When slicing the dendrogram, you can look at the two plots (WSSPlot and AveSilPlot) that show the relatively good quality of the partitions. The first shows the Within Sum of Square (WSS) for each partition and you can use the elbow approach to select a partition. The second graph shows the average silhouette width. This index is between -1 and 1. The closer it is to 1, the better the partition. See \link[cluster]{silhouette}.

If you want to cut the dendogram to a different dlim or number of clusters, you can do so without re-running  \code{mapclust() with \link[SPARTAAS]{mapclust_cut_tree}}.
}
\usage{
mapclust(
    coord,
    var,
    label = "Nolabel",
    iter = 20,
    Plabel = TRUE,
    lonlat = TRUE,
    positive_var = FALSE,
    n = NULL
)
}

\arguments{
  \item{coord}{
The x (longitude) and y (latitude) coordinates of the data.frame or matrix dimension 2.
}
  \item{var}{
The regionalised variable(s) of interest
}
  \item{label}{
(optional) The names of the samples or an identifier. Must be a factor.
}
  \item{iter}{
The number of iterations. The number of different dlim you want to test (must be greater than 10).
}
  \item{Plabel}{
Logical parameter to activate or not the printing of labels on the dendrogram.
}
  \item{lonlat}{
Logical parameter to activate or not the cartography in lonlat system with leaflet (base map).
}
  \item{positive_var}{
logical parameter that defines whether your variable of interest is positive or not.
}
  \item{n}{
Number of groups. If NULL, you can select the number of groups by clicking on the dendrogram.
}
}
\details{
Dlim is the selected minimum distance from the sample to the patch centre: to identify patches (units are those of coordinates). The dlim is automatically initialised and does not need to be set by the user. The minimum data is a data frame or matrix with at least 3 columns.
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
\item{var}{The regionalized variable(s) data you have used}
\item{label}{The label vector you have used}
\item{density}{The estimated density based on var. Equal to var if you used a unidimensional density variable.}
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

\references{
Bellanger L., Coulon A. and Husi P. (2021) Determination of cultural areas based on medieval pottery using an original divisive hierarchical clustering method with geographical constraint (MapClust), Journal of Archaeological Science, Volume 132 \doi{10.1016/j.jas.2021.105431}.

Bellanger L., Husi P., Laghzali Y. (2015). Spatial statistic analysis of dating using pottery: an aid to the characterization of cultural areas in West Central France. In : Traviglia A. ed., Across Space and Time, Proceedings of the 41th International Conference on Computer Applications and Quantitative Methods in Archaeology (CAA-2013), Perth (Australie), Amsterdam University Press : 276-282.

Woillez M., Poulard J.C., Rivoirard J., Petitgas P., Bez N. (2007). Indices for capturing spatial patterns and their evolution in time,
with application to European hake (Merluccius merluccius) in the
Bay of Biscay. ICES J. Mar. Sci. 64, 537-550.

Woillez M., Rivoirard J. and Petitgas P. (2009) Notes on survey-based spatial indicators for monitoring fish populations, Aquatic Living Resources, 22 :155-164.
}

\examples{
\donttest{
###################################
## loading data
 library(SPARTAAS)
 data(datarcheo)
 data(datacancer)

###################################
### Example: 1
## Function "mapclust"
# object <- mapclust( coord = ..., var = ..., label = ...)

classification <- mapclust(datarcheo$coord, datarcheo$var, datarcheo$label, n=4)

#Global dendrogram
 classification$dendrogram
#Cut dendrogram
 classification$cuttree
#silhouette of selected partition
 classification$silhouette


#You can cut the dendrogram for another dlim
 NewCut <- mapclust_cut_tree(classification, dlim=0.30)

#See evaluation using Silhouette width by running:
 NewCut$silhouette
 #If the plot is empty try to increase the height of the window (full screen)

#See summary of the data by running:
 summary(NewCut$silhouetteData)


###################################
## kmeans comparison
# pepare data (only geographical data)
 datakmeans <- datarcheo$coord

#kmeans
 number_cluster <- 4
 cl <- kmeans(datakmeans, number_cluster)
 plot(datakmeans, col = cl$cluster)


}
}
