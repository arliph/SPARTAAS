\name{timerange}
\alias{timerange}

%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Plot the timerange of observations sorted by cluster.
}
\description{
Vizualisation of cluster timerange.}
\usage{
timerange(
  data,
  cluster = NULL,
  add = NULL,
  density = NULL,
  color = NULL,
  reorder = FALSE
)
}

%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{data}{
data.frame (id, inf, sup) for each observation}
  \item{cluster}{
vector number cluster of observations}
  \item{add}{
data.frame of information to display on hover.}
  \item{density}{
vector of the density for each observation.}
  \item{color}{
vector of the color. if you want the colors to correspond to something other than clusters. Character vector of the same length as the number of observations}
  \item{reorder}{
Logical for reorder the color. If TRUE, the first color is for the leftmost cluster on the graph. If FALSE, the first color is for cluster number 1, wherever it is..}
}

\value{
The function returns a list.
  \item{plot }{The timerange plot.}
}

\author{
A. COULON

B. MARTINEAU

L. BELLANGER

P. HUSI
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
library(SPARTAAS)
data <- data.frame(
  id = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),
  Lower_bound = c(400,401,401,350,500,460,482,432,399,
    489,750,740,704,700,758,789,802,755,750,820),
  Upper_bound = c(550,689,755,700,780,700,700,699,650,
    850,1100,1100,1010,889,999,999,1050,1002,1000,1100)
)

cluster = c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2)

add <- data.frame(
  Site = c("Angers","Angers","Angers","Angers","Angers",
    "Angers","Angers","Angers","Angers","Angers",
    "Blois","Blois","Blois","Blois","Blois",
    "Blois","Blois","Blois","Blois","Blois")
)

timerange(data, cluster, add)

## with sub group (cluster 1 is sub divided in 2: 1.1 and 1.2)
cluster_with_sub = c(1.1,1.1,1.1,1.1,1.1,1.2,1.2,1.2,1.2,1.2,2,2,2,2,2,2,2,2,2,2)

timerange(data, cluster_with_sub, add)

## with density
density <- c(32,34,35,19,9,25,19,29,28,18,10,13,9,10,9,6,3,7,7,1)
timerange(data=data, cluster=cluster, density=density)

}

