\name{seriograph}
\alias{seriograph}
\title{
Plot seriograph (B. DESACHY).
}

\description{
  Visualization of contingency data over time. \strong{Rows} must be individuals (archaeological site,...) and \strong{columns} must be categories (type,...).
}

\usage{
seriograph(cont, order, insert, show, permute, col_weight)
}

\arguments{
  \item{cont}{
    Contingency table or hclustcompro object. Note: Your contingency table must have the rows sorted in chronological order. (the order parameter allows you to change the order of the rows if necessary)
}
  \item{order}{
     Vector for change the order of the rows (use row's names or cluster names if cont is a hclustcompro object, as a character vector). The oldest one at the bottom. Missing names will not be ploted. You can remove row by simply remove the name in the vector.
}
  \item{show}{
     The element to plot. This should be (an unambiguous abbreviation of) one of "both", "EPPM" or "frequency".
}
  \item{permute}{
     Logical for permute columns in order to show seriation.
}
  \item{col_weight}{
     Logical for activiate or not the coloration of the last column: weight.
}
  \item{insert}{
     Vector with the position after where you want insert one or more Hiatus. Could be a list with two vector: position and label to print instead of hiatus. (see last examples)
}
}
\details{
  \strong{Seriograph}\cr\cr
  We chose the seriograph (B. DESACHY). This tool makes it possible to highlight artisanal evolutions over time as well as to understand commercial relations thanks to imported potteries. In this representation, the oldest element is at the bottom. The percentages of each pottery category are displayed. The percentages are calculated independently for each class. The percentage display allows you to compare the different classes but does not provide information on the number of individuals per class. To fill this gap, the proportion of each class of their workforce is displayed on the seriograph (weight).

  We can generalized this representation for other contingency data or with \link[SPARTAAS]{hclustcompro} object.

  The visualization of positive deviations from the average percentage allows us to observe a series that results from changes in techniques and materials dedicated to ceramic manufacturing over time.

  \strong{Positive deviation from the average percentage (EPPM in French)}\cr\cr
  The average percentage is calculated for each ceramic category (columns) on the total number of accounts (all classes combined). From the average percentage we recover for each category and for each rows the difference between the percentage of the category in the class with the average percentage. The EPPM corresponds to the notion of independence deviation (between rows and columns, between categories and time classes) in a chi-square test approach. Although this approach is fundamental in statistical analysis, independence deviations are here purely indicative and are not associated with a p_value that could determine the significance of deviations.

    \strong{Weight}\cr\cr
  Weight is the number of observations divided by the total number of observations. It indicates for each row the percentage of the data set used to calculate the frequencies of the elements (row).

  \strong{Permutation}\cr\cr
  order argument:\cr
  The rows of the contingency table are initially in the order of appearance (from top to bottom). It must be possible to re-order the classes in a temporal way (You can also order as you want your contingency table).

  permute argument:\cr
  In addition, it is possible to swap ceramic categories (contingency table columns) in order to highlight a serialization phenomenon. Matrix permutation uses an algorithm called "reciprocal averages". Each line is assigned a rank ranging from 1 to n the number of lines. A barycentre is calculated for each column by weighting according to the row rank. Finally, the columns are reorganized by sorting them by their barycentre.

  \strong{Insert}\cr\cr
  It's possible to insert a row in the seriograph in order to represent a archeological hiatus or other temporal discontinuities.
}
\value{
The function returns a list (class: seriograph).
\item{seriograph }{The seriograph plot}
\item{dendrogram}{If cont is a hclustcompro object return the dendrogram with the period order as label}
\item{contingency }{Data frame of the contingencies data group by cluster}
\item{frequency }{Data frame of the frequencies data group by cluster}
\item{ecart }{Data frame of the gap data group by cluster}
}

\author{
Desachy Bruno. Le sériographe EPPM : un outil informatisé de sériation graphique pour tableaux de comptages. In: Revue archéologique de Picardie, n°3-4, 2004. Céramiques domestiques et terres cuites architecturales. Actes des journées d'étude d'Amiens (2001-2002-2003) pp. 39-56. DOI : https://doi.org/10.3406/pica.2004.2396.

A. COULON

L. BELLANGER

P. HUSI
}

\examples{
\donttest{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
library(SPARTAAS)

## network stratigraphic data (Network)
network <- data.frame(
  nodes = c("AI09","AI08","AI07","AI06","AI05","AI04","AI03",
  "AI02","AI01","AO05","AO04","AO03","AO02","AO01","APQR03","APQR02","APQR01"),
  edges = c("AI08,AI06","AI07","AI04","AI05","AI01","AI03","AI02","","","AO04","AO03",
  "AO02,AO01","","","APQR02","APQR01","")
)
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
## default
seriograph(cont)

seriograph(cont,show = "EPPM")
seriograph(cont,show = "frequency")

## clustering <- hclustcompro(cont,network,alpha=0.7,k=7) # number of cluster 7
seriograph(clustering)

## change order with cluster name (letters on dendrogram) to sort them in a chronological order
seriograph(clustering,order = c("C","F","A","G","E","B","D"))

## Don't allow permutation of columns
seriograph(clustering,order = c("C","F","A","G","E","B","D"),permute = FALSE)

## Don't allow coloration
seriograph(clustering,order = c("C","F","A","G","E","B","D"),col_weight = FALSE)

## insert Hiatus (position, 1 -> after first row from bottom: oldest)
seriograph(clustering,order = c("C","F","A","G","E","B","D"),insert = 2)
seriograph(clustering,order = c("C","F","A","G","E","B","D"),insert = c(2,3))

## insert costum label element
insert <- list(
  position = c(2,3),
  label = c("Hiatus.100years","Missing data")
)
seriograph(clustering,order = c("C","F","A","G","E","B","D"),insert = insert)
}
}

