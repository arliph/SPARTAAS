\name{hclustcompro}
\alias{hclustcompro}
\alias{perioclust}
\title{hclustcompro}
\usage{
hclustcompro(
  D1,
  D2,
  alpha="EstimateAlphaForMe",
  k=NULL,
  title="notitle",
  method="ward.D2",
  suppl_plot=TRUE
)
}
\description{
Compromised Hierarchical bottom-up clustering method.
The method use two sources of informations. The merging of the two data sources is done by a parameter (alpha) which allows to weight each source.
Formula: D_alpha = alpha * D1 + (1-alpha) * D2\cr
}

\arguments{
  \item{D1}{
First dissimilarity matrix (square matrix) or distance matrix. Could be a contingency table (see \link[SPARTAAS]{CAdist}). A factorial correspondences analysis is carried out and the distances are used (Chi-square Metric).
}
  \item{D2}{
Second dissimilarity matrix (square matrix), same size than D1, or distance matrix.
}
  \item{alpha}{
The mixing parameter in order to generate the D_alpha matrix. Formula: D_alpha = alpha * D1 + (1-alpha) * D2
}

  \item{k}{
The number of cluster you want.
}

  \item{title}{
The title to display on the dendrogram plot.
}

  \item{method}{
The agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
}

  \item{suppl_plot}{
Logical define whether to display additional plots (WSS and averages sil plot).
}
}

\details{
\strong{CAH}\cr
Data fusion (alpha optimal value parameter see \link[SPARTAAS]{hclustcompro_select_alpha}). It is necessary to define the appropriate proportion for each data source. This is the first sensitive point of the method that the user must consider. A tool is provided to guide his decision.

\strong{Cut dendrogram}\cr
The division into classes, and subclasses, is the second crucial point. It must be done on the basis of knowledge of the study area and some decision support tools such as the cluster silhouette or the calculation of intra cluster variability (WSS: Within Sum of Square). You can use \link[SPARTAAS]{hclustcompro_subdivide} in order tu sub-divide a cluster into sub-cluster.
}

\value{
The function returns a list (class: hclustcompro_cl).
\item{D1 }{First dissimilarity matrix (square matrix)}
\item{D2 }{Second dissimilarity matrix (square matrix)}
\item{D_alpha }{The matrix use in the CAH result of the mixing of the two matrix (D1 and D2)}
\item{alpha }{Alpha}
\item{tree }{An object of class hclust which describes the tree produced by the clustering process (see \link[stats]{hclust})}
\item{cluster }{The vector of cluster of the selected partition}
\item{cutree }{Plot of the cut dendrogram}
\item{call }{How you call the function}
\item{cont }{Original contingency data (if D1 is a contingency table)}
}

\author{
The hclust function is based on Fortran code contributed to STATLIB by F. Murtagh.

A. COULON

L. BELLANGER

P. HUSI
}

\examples{
library(SPARTAAS)
data(datangkor)

#network stratigraphic data (Network)
network <- datangkor$stratigraphy

#contingency table
cont <- datangkor$contingency

#obtain the dissimilarities matrices
distance <- CAdist(cont, nPC = 11)
constraint <- adjacency(network)

#You can also run hclustcompro with the dist matrix directly
hclustcompro(D1 = distance, D2 = constraint, alpha = 0.7, k = 4)
}
