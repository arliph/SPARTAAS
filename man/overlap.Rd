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

\item{overlap_matrix }{The dissimilarity matrix with for each couple i,j the value of the overlap index.}

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

}
