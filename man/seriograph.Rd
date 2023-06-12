\name{seriograph}
\alias{seriograph}
\title{
Plot seriograph (B. DESACHY).
}

\description{
  Visualization of contingency data over time. \strong{Rows} must be individuals (archaeological site,...) and \strong{columns} must be categories (type,...).
}

\usage{
seriograph(cont, order, insert, show, permute)
}

\arguments{
  \item{cont}{
    Contingency table or hclustcompro object. Note: Your contingency table must have the rows sorted in chronological order (the order parameter allows you to change the order of the rows if necessary).
}
  \item{order}{
     Vector to change the order of the rows (use row names or cluster names if cont is a hclustcompro object, as a character vector). The oldest (at the bottom) must be at the end of the vector. Missing names are not plotted. You can remove a row by simply removing the name in the vector.
}
  \item{show}{
     The element to be plotted. This should be (a unique abbreviation of) one of 'both', 'EPPM' or 'frequency'.
}
  \item{permute}{
     Logical for permute columns in order to show seriation.
}

  \item{insert}{
     Vector with the position after which you want to insert one or more hiatuses. Could be a list with two vectors: position and label to be printed instead of hiatus (see last examples).
}
}
\details{
  \strong{Seriograph}\cr\cr
  We have chosen the serigraph (DESACHY 2004). This tool makes it possible to highlight the evolution of ceramics over time as well as to understand the commercial relations thanks to the imported ceramics. The percentages of each category of ceramics per set are displayed. The percentages are calculated independently for each set (row). The display of the percentages allows comparison of the different sets but does not provide information on the differences in numbers. To fill this gap, the proportion of the numbers in each class is displayed on the seriograph (weight column).

  We can generalized this representation for other contingency data or with \link[SPARTAAS]{hclustcompro} object.

  The visualization of positive deviations from the average percentage allows us to observe a series that results from changes in techniques and materials dedicated to ceramic manufacturing over time.

  In order to facilitate the exploitation of the data tables, we propose here a computerised graphic processing tool (EPPM serigraph - for Ecart Positif aux Pourcentages Moyens - positive deviation from the average percentage), which does not require specialised statistical skills and is adapted to the case of stratified sites, where the study of the evolution of artefacts can be based on the relative chronology provided by the excavation.

The treatment consists firstly of transforming this table of counts into a table of percentages, the total number in each set (each row) being reduced to 100; these are the proportions, or frequencies, of the types in the sets are thus compared.

The display of positive deviations from the mean percentages (EPPM) shows in black on a lighter background the percentage shares that are higher than the mean percentage of the variable, so as to highlight the most significant part of the values in the table.This display is simply adapted to the seriograph: when a percentage is greater than the average percentage of the type, the excess share (called here EPPM: positive deviation from the average percentage) is shown in black, centred around the axis of the type, on the grey background of the percentage bar.

The table is then transformed into a graphic matrix where these percentages are expressed, for each type, by horizontal bars centred on the column axis. When the rows are ordered chronologically, the silhouette formed by the superposition of these frequency bars bars makes it possible to visualise the evolution over time of the type concerned.

The display of the percentages allows comparison of the different sets but does not provide information on the differences in numbers. To fill this gap, the proportion of the numbers in each class is displayed on the seriograph (weight column).

The processing technique applies to sets whose chronological order is not known; the lines of the graph are to be reorganised so as to obtain boat-shaped silhouettes following the hypothesis of a chronological evolution corresponding to the seriation model.

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

\references{
Desachy B. (2004). Le sériographe EPPM : un outil informatisé de sériation graphique pour tableaux de comptages. In: Revue archéologique de Picardie, n°3-4, 2004. Céramiques domestiques et terres cuites architecturales. Actes des journées d'étude d'Amiens (2001-2002-2003) pp. 39-56 \doi{10.3406/pica.2004.2396}
}
\author{
B. DESACHY

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
data(datangkor)

## network stratigraphic data (Network)
network <- datangkor$stratigraphy

## contingency table
cont <- datangkor$contingency

## default
seriograph(cont)

seriograph(cont,show = "EPPM")
seriograph(cont,show = "frequency")

## Don't allow permutation of columns
seriograph(cont, permute = FALSE)

## insert Hiatus (position, 1 -> after first row from bottom: oldest)
seriograph(cont,insert = 2)
seriograph(cont,insert = c(2,3))

## insert costum label element
insert <- list(
  position = c(2,3),
  label = c("Hiatus.100years","Missing data")
)
seriograph(cont,insert = insert)

## change order with cluster name (letters on dendrogram) to sort them in a chronological order
seriograph(cont,order=c("AI03","AI09","AI01","AI02","AI04","APQR01","AO05",
"AO03","AI05","AO01","APQR02","AI07","AI08","AO02","AI06","AO04","APQR03"))
## by omitting the row names, you delete the corresponding rows
seriograph(cont,order=c("AI02","AI08","APQR03","AI09"))
}
}

