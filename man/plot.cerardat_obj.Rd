\name{plot.cerardat_obj}
\alias{plot.cerardat_obj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot cerardat model
}
\description{
plot cerardat model
}
\usage{
\method{plot}{cerardat_obj}(x, which = NULL,
                  col1 = rgb(0.93,0.23,0.23,0.5), col2 = "black",
                  xlim=NULL, ylim=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a \link[SPARTAAS]{cerardat} output.
}
  \item{which}{
  Vector containing the plots you want (row number on the contingence table: numeric).
}
  \item{col1}{
  Color of the Event curve.
}
  \item{col2}{
  Color of the Accumulation curve.
}
  \item{xlim}{
  Two numeric values, specifying the left limit and the right limit of the scale.
}
  \item{ylim}{
  Two numeric values, specifying the lower limit and the upper limit of the scale.
}
  \item{...}{
  other parameters to be passed through to plotting functions.
}
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

## The function is currently defined as

}
