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
clustering <- hclustcompro(D1 = distance, D2 = constraint, alpha = 0.7, k = 7) #number of cluster 7
clustering <- hclustcompro_subdivide(clustering,cluster = 5, nb_class = 2)

#subdivide more than one cluster
clustering2 <- hclustcompro(D1 = distance, D2 = constraint,0.7,k=7) #number of cluster 7
clustering2 <- hclustcompro_subdivide(clustering2,cluster = c(5,7), nb_class = c(2,2))
}
