\name{cerardat}
\alias{cerardat}

\title{
cerardat
}
\description{
The methodology is based on a statistical and visual approach using two estimated density curves to date each archaeological context. The statistical procedure required two steps, each leading to the construction of a density curve. The first allowed us to estimate a date corresponding to the terminus post quem of the context, a cursor reflecting an event dated in calendar time. This statistical tool allows the archaeologist to easily visualise and analyse chronological patterns.
}
\usage{
cerardat(df, row.sup, date, nf = NULL, confidence = 0.95, graph = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
The data (data.frame) is a contingency table with the observations in the rows and the technical groups in the columns.
}
  \item{row.sup}{
Index of supplementary rows in df (vector).
}
  \item{date}{
The dates of each observation or NA (vector).
}
  \item{nf}{
an integer representing the number of axes retained in the correspondence analysis. If NULL, it is automatically chosen to keep a number corresponding to at least 60\% of the inertia.
}
  \item{confidence}{
The desired confidence interval (0.95 for 95\%).
}
  \item{graph}{
logical to display the plots or not.
}


}
\details{
The corpus data is a contingency table with the observations in the rows and the technical groups in the columns. There are two types of observations: the reference corpus observations and the supplementary observations. The supplementary rows (observations) are identified by the argument `row.sup`.

\strong{step 1: modelling events dated in calendar time (dateEv)}\cr
This step involves estimating the date of an event recorded in the ground (an archaeological context for the archaeologist) from the pottery assemblage of which it is composed, by fitting a regression model that relates a known date in calendar time, such as the date of issue of a coin, to its pottery profile. The reference corpus used to fit the regression model. We then used the previously fitted model to calculate a predicted value for contexts not included in the reference corpus, sometimes stratigraphically separated or poorly documented, with a 95\% confidence interval for the predicted date.

A correspondence analysis (\link[FactoMineR]{CA}) was carried out to summarize the information in the reference corpus data. We then kept only the first factorial axes. In this way, our contingency table becomes a reduced size table, an incomplete reconstruction of the data. This principle is used in many factor analysis techniques to reduce the number of explanatory variables in the linear regression model.

After estimating the beta parameters of the model using the classical results of multiple regression analysis and checking that the model fits the data correctly, we can deduce the estimated date of an observation and also predict the date of another observation that has no coins and is therefore not dated.

\strong{step 2: from event time (dateEv) to accumulation time (dateAc)}\cr
We used the results of the first step and the properties of the CA to obtain an estimate of the date of each fabric. We could then define the archaeological time represented as dateAc, in other words the accumulation time of a context, as the weighted sum of the fabric dates; the weights being the proportions of MINVC of each fabric in the context. Assuming that the random variables dateEvj are independent, the distribution of the accumulation time of each context can be approximated by the Gaussian mixture. In this way, for each context, we obtained a plurimodal density curve representing the estimated law of accumulation time based on the mixture of normal densities (dates of each tissue). By definition, the area under the density curve has a value of 1 (i.e. 100\%).

\strong{date}\cr
In order to estimate a date for the context, it is essential to refer to objects that have been dated by another source, such as coins. These contexts were selected on a very strict basis for their chronostratigraphic reliability, level of domestic occupation or enclosures with long urban stratigraphic sequences, thereby minimising any bias associated with the disparity between the date of the coin and that of the context.

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{prediction }{Estimated date for archaeological context (event: dateEV and accumulation: dateAC) with confidence interval.}
\item{date_gt }{Estimated date for technical groups with confidence interval. (use for dateAC)}
\item{lm }{Linear model on the components of the correspondance analysis.}

\item{predict_obj_row }{date prediction of archaeological contexts (rows) using \link[stats]{predict.lm}.}
\item{predict_obj_col }{date prediction of technical groups (columns) using \link[stats]{predict.lm}.}
\item{cont_gt }{Contingency table of the reference corpus.}

\item{statistical.summary }{Statistical summary of the model:\cr
Adjusted R-squared\cr
R-squared\cr
sigma (Residual standard error)\cr
The Shapiro-Wilks test is used to verify the normality of the residuals.\cr
The Durbin-Watson test checks for first order autocorrelation.\cr
The Breusch-Pagan test checks for heteroscedasticity.
}

\item{obs_ca_eval }{Quality of row representation in the correspondence analysis.}

\item{check_ref }{Plot of estimated dates (with confidence interval) and real dates of reference data. Only when the real date is known.}
\item{check_sup }{Plot of estimated dates (with confidence interval) and real dates of supplementary data. Only when the real date is known.}

\item{Shapiro_Wilks }{Summary of the Shapiro-Wilks test. see \link[stats]{shapiro.test}.}
\item{Durbin_Watson }{Summary of the Durbin-Watson test. see \link[lmtest]{dwtest}.}
\item{Breusch_Pagan }{Summary of the Breusch-Pagan test. see \link[lmtest]{bptest}.}



%% ...
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
\donttest{
data("datacerardat")

resultat = cerardat(df = datacerardat$df,
           row.sup = datacerardat$row.sup,
           date = datacerardat$date,
           nf = NULL,
           confidence = 0.95,
           graph = TRUE
        )

resultat
#The Shapiro-Wilks test is used to verify the normality of the residuals.
#The Durbin-Watson test checks for first order autocorrelation.
#The Breusch-Pagan test checks for heteroscedasticity.



#See the first plot
plot(resultat,
     which = 1,
     col1=rgb(0.93,0.23,0.23,0.5),
     col2="black",
     xlim=NULL,
     ylim=c(0,0.03)
    )

#See the first ten plots
#plot(resultat,
#     which = 1:10,
#     col1=rgb(0.93,0.23,0.23,0.5),
#     col2="black",
#     xlim=NULL,
#     ylim=c(0,0.03)
#    )

#See all plots
#plot(resultat,
#     which = NULL,
#     col1=rgb(0.93,0.23,0.23,0.5),
#     col2="black",
#     xlim=NULL,
#     ylim=c(0,0.03)
#    )

#You can extract the plots and find them in the directory :
paste0(getwd(),"/figures")
#With the 'extract_results' function
#extract_results(resultat,width=480, height=480, path="figures")
}
}
