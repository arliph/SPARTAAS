\name{overlap}
\alias{overlap}
\title{
Temporal overlap index
}
\description{
The overlap index is the ratio between internal overlap and total overlap over time. We define the total overlap limit as follows: the minimum of the lower limits of the pair of individuals and the maximum of the upper limits. We define the limit of the internal overlap as follows: the maximum of the lower limits and the minimum of the upper limits.
}

\details{
The lower and upper limits must be number.

The dissimilarity between time ranges is calculated from the ratio of the Overlap of the time ranges (Distance in the case where the two ranges are disjoint) by the cumulative extent of the two.

As the ratio is bounded between -1 and 1, we add 1 to make it positive and normalise it so that it is between 0 and 1.

This overlap index must then be transformed into a dissimilarity index between sets. We do 1 - ratio to do this. It is equal to 0 if the two time ranges are identical and 1 when they are infinitely distant.
}
\usage{
overlap(temporal)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{temporal}{
A data frame with tree columns: the name of the element, the lower limit and the upper limit.
}
}

\value{

\item{D }{The dissimilarity matrix base on the overlap index.}

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
data(dataceram)
dist <- overlap(dataceram$timerange)

}
