\name{hclust}
\alias{hclust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hierarchical Clustering
}
\description{
Overload of \link[stats]{hclust} for dealing with two dissimilarities matrices.
Hierarchical cluster analysis on a set of dissimilarities and methods for analyzing it.
}
\usage{
hclust(d, method = "complete", members = NULL, d2 = NULL, alpha = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
a dissimilarity structure as produced by dist.
}
  \item{method}{
the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
}
  \item{members}{
NULL or a vector with length size of d. See the ‘Details’ section.
}
  \item{d2}{
a second dissimilarity structure as produced by dist.
}
  \item{alpha}{
The mixing parameter in order to generate the D_alpha matrix on which the classical hclust method is applied. Formula: D_alpha = alpha * d + (1-alpha) * d2.
}
}
\details{
Data fusion (alpha optimal value parameter see \link[SPARTAAS]{hclustcompro_select_alpha}. It is necessary to define the appropriate proportion for each data source. This is the first sensitive point of the method that the user must consider. A tool is provided to guide his decision.
}
\value{
\link[stats]{hclust}
}

\author{
The hclust function is based on Fortran code contributed to STATLIB by F. Murtagh.

A. COULON

L. BELLANGER

P. HUSI
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (d, method = "complete", members = NULL, d2 = NULL,
    alpha = NULL)
{
    if (!is.null(d2)) {
        if (!length(d) == length(d2)) {
            stop("d and d2 have not the same size.")
        }
        if (is.null(alpha)) {
            sa <- hclustcompro_select_alpha(d, d2, method = method,
                resampling = FALSE)
            alpha <- sa$alpha[1]
        }
        alpha <- as.numeric(alpha)
        if (!(alpha > 0 & alpha < 1)) {
            warning("Alpha must be between 0 and 1.")
            sa <- hclustcompro_select_alpha(d, d2, method = method,
                resampling = FALSE)
            alpha <- sa$alpha[1]
        }
        d <- dist(alpha * d + (1 - alpha) * d2)
    }
    stats::hclust(d, method, members)
  }
}

