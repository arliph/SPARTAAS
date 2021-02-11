\name{hclustcompro_app}
\alias{hclustcompro_app}
\title{
Launch the shiny application.
}
\description{
see \link[SPARTAAS]{hclustcompro}, \link[SPARTAAS]{hclustcompro_select_alpha}, \link[SPARTAAS]{seriograph}.
You can also check the wiki on the application.}
\usage{
hclustcompro_app()
}

\value{No return value}

\author{
A. COULON

L. BELLANGER

P. HUSI
}



\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
\donttest{
library(SPARTAAS)
if(interactive()){
  hclustcompro_app()
}
}

}

