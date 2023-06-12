\name{cerardat_estim_nf}
\alias{cerardat_estim_nf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Number of axes to keep
}
\description{
Estimate the correct number of axes to keep in the regression model.
}
\usage{
cerardat_estim_nf(df, row.sup, date)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
The data (data.frame) is a contingency table with the observations in rows and the technical groups in columns.
}
  \item{row.sup}{
Index of supplementary rows in df (vector).
}
  \item{date}{
The dates of each observation or NA (vector).
}
}
\value{


\item{nf }{Number of axes to keep (minimal PRESS value)}
\item{MSE }{plot of the Mean Squared Error.}
\item{PRESS }{plot of the PRediction Error Sum Of Squares.}
\item{adj.R_sq }{plot of the Coefficient of determination R².}
\item{data }{data frame of MSE, PRESS and R_sq values.}

}
\references{
Bellanger L. and Husi P. (2012) Statistical tool for dating and interpreting archaeological contexts using pottery. Journal of Archaeological Science, Elsevier, 39 (4), pp.777-790. doi:10.1016/j.jas.2011.06.031.
}
\author{
A. COULON

L. BELLANGER

P. HUSI
}


\examples{
library(SPARTAAS)
data(datacerardat)
res = cerardat_estim_nf(datacerardat$df, datacerardat$row.sup, datacerardat$date)

#Number of axes to keep (minimal PRESS value)
res$nf

#the plots
res$MSE
res$PRESS
res$adj.R_sq
}

