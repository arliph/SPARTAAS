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
  method="ward.D2"
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

#network stratigraphic data (Network)
network <- data.frame(
  nodes = c("AI09","AI08","AI07","AI06","AI05","AI04","AI03",
  "AI02","AI01","AO05","AO04","AO03","AO02","AO01","APQR03","APQR02","APQR01"),
  edges = c("AI08,AI06","AI07","AI04","AI05","AI01","AI03","AI02","","","AO04","AO03",
  "AO02,AO01","","","APQR02","APQR01","")
)
#contingency table
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
#obtain the dissimilarities matrices
distance <- CAdist(cont, nPC = 11)
constraint <- adjacency(network)

#You can also run hclustcompro with the dist matrix directly
hclustcompro(D1 = distance, D2 = constraint, alpha = 0.7, k = 4)


}
