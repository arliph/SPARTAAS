\name{extract_results}
\alias{extract_results}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generates all plots (in jpeg format) in a subfolder
}
\description{
Generates all plots (in jpeg format) in a subfolder (relative path from the working directory).
}
\usage{
extract_results(cerardat, width=480, height=480, path="figures",
                  col1 = rgb(0.93,0.23,0.23,0.5), col2 = "black",
                  xlim=NULL, ylim=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cerardat}{
a \link[SPARTAAS]{cerardat} output.
}
  \item{width}{
  the width of the graphics region in pixels.
}
  \item{height}{
  the width of the graphics region in pixels.
}
  \item{path}{
  directory where plots are exported (relative path from the working directory). Automatically generates two sub-folders 'ref' for the reference plot and 'sup' for the supplementary  plot.
}
  \item{col1}{
  Color of the the Event curve (color name, hex code, \link[grDevices]{rgb} function \{grDevices\} for transparency).
}
  \item{col2}{
  Color of the the Accumulation curve (color name, hex code, \link[grDevices]{rgb} function \{grDevices\} for transparency).
}
  \item{xlim}{
  Two numeric values, specifying the left limit and the right limit of the scale.
}
  \item{ylim}{
  Two numeric values, specifying the lower limit and the upper limit of the scale.
}
}


\author{
A. COULON

L. BELLANGER

P. HUSI
}

\examples{
\donttest{
data("datacerardat")

resultat = cerardat(df = datacerardat$df,
           row.sup = datacerardat$row.sup,
           date = datacerardat$date,
           nf = NULL,
           confidence = 0.95,
           graph = TRUE
        )
#Generates all plots (in jpeg format) in a folder 'figures'
#You will find the jpeg in 'ref' and 'sup' subfolder
#extract_results(resultat,width=480, height=480, path="figures",
#     col1=grDevices::rgb(0.93,0.23,0.23,0.5),
#     col2="black",
#     xlim=NULL,
#     ylim=c(0,0.03)
#    )

#You can extract the plots and find them in the directory :
paste0(getwd(),"/figures")
}
}
