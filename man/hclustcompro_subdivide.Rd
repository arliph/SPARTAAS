\name{hclustcompro_subdivide}
\alias{hclustcompro_subdivide}
\title{
subdivide a cluster after runing hclustcompro.
}

\description{
Allow the user to subdivide one cluster into sub clusters. This function work only with 'hclustcompro_cl' object return by hclustcompro function.
}

\usage{
hclustcompro_subdivide(hclustcompro_cl,cluster,nb_class)
}

\arguments{
  \item{hclustcompro_cl }{A hclustcompro_cl object}
  \item{cluster }{The number of the cluster. Numbered from left to right on the dendrogram (1, 2, ...)}
  \item{nb_class }{The number of sub-clusters you want}
}

\value{
  \item{hclustcompro_cl }{A new hclustcompro_cl object updated see \link[SPARTAAS]{hclustcompro}}
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
data(datangkor)

#network stratigraphic data (Network)
network <- datangkor$stratigraphy

#contingency table
cont <- datangkor$contingency

#obtain the dissimilarities matrices
distance <- CAdist(cont, nPC = 11)
constraint <- adjacency(network)

#You can also run hclustcompro with the dist matrix directly
clustering <- hclustcompro(D1 = distance, D2 = constraint, alpha = 0.7, k = 7) #number of cluster 7
clustering <- hclustcompro_subdivide(clustering,cluster = 5, nb_class = 2)

#subdivide more than one cluster
clustering2 <- hclustcompro(D1 = distance, D2 = constraint,0.7,k=7) #number of cluster 7
clustering2 <- hclustcompro_subdivide(clustering2,cluster = c(5,7), nb_class = c(2,2))
}
