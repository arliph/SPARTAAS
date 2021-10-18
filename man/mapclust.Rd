\name{mapclust}
\alias{mapclust}
\encoding{utf-8}
\title{
Divise hierarchical Clustering using Spatialpatches algorithm.
}

\description{

\href{https://spartaas.frama.io/r-package/mapclust.html}{\code{Get started with MapClust}}

This function performs a divisive hierarchical clustering on a regionalized variable using spatial patches algorithm (Woillez and al. 2007; Woillez, Rivoirard and Petitgas 2009). It is a hierarchical top-down clustering with geographical constraint. It is possible to cut the tree by clicking on the dendrogram at the level you want. Results include description of the clusters and graphics. When cutting the dendrogram you can look at the two plots (WSSPlot and AveSilPlot) which shows the relative good quality of the partitions. The first one shows the Within Sum of Square (WSS) for each partition and you can use the Elbow approach in order to select one partition. The second plot show the average silhouette width. This index is between -1 and 1. The closer it is to 1 the more the partition is good. See the fonction help page of \code{silhouette{cluster}}


See \link[SPARTAAS]{mapclust_cut_tree} if you want cut the dendrogram to another dlim or number of clusters without re run \code{mapclust()}.
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
The x-coordinates (longitude) and the y-coordinates (latitude) data.frame or matrix dimension 2.
}
  \item{var}{
The regionalized variable(s) of interest
}
  \item{label}{
(optional) The names of the samples or an id. Must be factor.
}
  \item{iter}{
The number of iterations. The number of different dlim you want to test. (must be higher than 10)
}
  \item{Plabel}{
Logical parameter for activate or not the print of labels on the dendrogram
}
  \item{lonlat}{
Logical parameter for activate or not the cartography in lonlat system with Leaflet (basemap)
}
  \item{positive_var}{
logical parameter defining whether your variable of interest is positive or not
}
  \item{n}{
Number of groups. If NULL you can select the number of groups by clicking on the dendrogram
}
}
\details{
Dlim is the select minimum distance from sample to patch centre: to identify patches (units are those of coordinates). The dlim is automatically initialized and the user does not need to use it. The minimal data are data frame or matrix with at least 3 columns.
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
\item{var}{The regionalized variable(s) data you used}
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

\references{
Bellanger L., Coulon A. and Husi P. (2021) Determination of cultural areas based on medieval pottery using an original divisive hierarchical clustering method with geographical constraint (MapClust), Journal of Archaeological Science, Volume 132 \href{https://doi.org/10.1016/j.jas.2021.105431}{doi:10.1016/j.jas.2021.105431}.

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
