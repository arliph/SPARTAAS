library(SPARTAAS)
suppressMessages(library(shinyjs))
library(shinyjqui)
suppressMessages(library(shinydashboard))
library(colorspace)
suppressMessages(library(dplyr))
library(tidyr)
library(ggplot2)
library(shinyWidgets)
library(ade4)
library(shinythemes)
library(shinycssloaders)
library(plotly)
library(scatterD3)
library(caret)


ui <- fluidPage(useShinyjs(),withMathJax(),style="padding-top: 150px;",theme = shinytheme("flatly"),
    tags$style(HTML('#sort,#data,#subdivise,#tri,#sub1,#sub2,
         #import2,#import3,#import3,#Finishimport,#back1,#back2{color:white;background-color:#337ab7;border-color:#337ab7}
         #back{color:white;background-color:#f0ad4e}
         #tabsetperso{background-color:#e5e5e5;}
         #help{text-align;right;}
         #tabsetperso{margin-bottom:20px;}
         .btn-default {
           background-color:#5da3a8;
           border-color:#5da3a8;
         }
         .dropdown-menu > li > a {
           background-color:#fff;
           border-color:#dce4ec;
         }
         #cerardat_plot_ref > div {display:flex;flex-wrap:wrap;width:94vw;}
         #cerardat_plot_ref > div > div.shiny-plot-output {width:450px !important;}
         #cerardat_plot_sup > div {display:flex;flex-wrap:wrap;width:94vw;}
         #cerardat_plot_sup > div > div.shiny-plot-output {width:450px !important;}
         #lmplot > div {display:flex;flex-wrap:wrap;width:65vw;}
         #lmplot > div > div.shiny-plot-output {width:450px !important;}
         .disabled{color:grey}
         .disabled:hover{color:grey;background-color:#e5e5e5 !important}
    ')),
    tags$head(
      tags$link(rel = "icon", type = "image/gif", href = "https://spartaas.gitpages.huma-num.fr/r-package/img/lambda.png"),
      tags$title("cerardat")
    ),
  #header panel
  absolutePanel(class = "panel panel-default",
    style="z-index: 2000;padding: 8px; background: #ecf0f1; opacity: 1;border-bottom: 1px solid #2c3e50;",
    top = 0, left = 0, right = 0,
    fixed = TRUE,
    h2("cerardat: Ceramic estimation date"),
    div(span(strong("SPARTAAS | cerardat")))
  ),
  #footer
  absolutePanel(style="z-index: 2000;padding: 0px; border-bottom: 0px solid #CCC; background: #fff;opacity: 1;",
                bottom = 0, left = 0,
                fixed = TRUE,
                div(span(strong("SPARTAAS")),span("[Bellanger,Coulon,Husi]",style ="font-size:14px;"))
  ),
  #go top panel
  absolutePanel(class = "panel-hover",
    style="opacity:0.8;z-index: 1000;padding: 8px; border-bottom: 1px solid #CCC;
    background: #FFFFEE;width:30px;margin-right:10px;",
    bottom = 20, right = 0,
    fixed = TRUE,
    HTML("<a href=\"#top\"><span class=\"glyphicon glyphicon-chevron-up\"></span></a>")
  ),
  #-----------------------------------#
  #      NavBar
  #-----------------------------------#
  navbarPage("cerardat",id="navbar",
    tabPanel(value = "1",span(icon("home", lib = "glyphicon"),strong("Home")),style="max-width:1200px;",
             br(),
             sidebarLayout(
               #-----------------------------------#
               sidebarPanel(style="",
                            h1("R Package:"),
                            br(),
                            HTML("<p>This method is part of the <a target='_blank' href='https://spartaas.gitpages.huma-num.fr/r-package/index.html'>SPARTAAS</a> package.</p>
                                                        <p>If you are interested you can install our R package available on the <a target='_blank' href='https://cran.r-project.org/package=SPARTAAS'>CRAN</a> and on <a target='_blank' href='https://github.com/arliph/SPARTAAS'>GitHub</a>.</p>

                                                        ")
               ),
               #-----------------------------------#
               mainPanel(
                 HTML("
                      <h3>Introduction</h3>
                      <p>Analyzing chronological patterns is one of the major issues in archaeology. How can the date of a specific
context be estimated? Is it possible to identify residual and intrusive material in it at the same time?
Numerous statistical methodological approaches have been developed and implemented to estimate
dates but have less often addressed the issue of socio-economic area or the functional interpretation of
contexts. This article deals with the construction and analysis of two different probability estimate
density curves of context dates using pottery. By contrasting the two curves we can define the boundaries
of the socio-economic area and make a chrono-functional interpretation of a context. This statistical
tool allows the archaeologist to visualize and analyze chronological patterns easily.</p>
<p>Dating an archaeological context is always delicate. It is
important to keep in mind what is to be dated, either a one-off
event occurring at a specific moment, or a process built up over
time. This study proposes datings based on the statistical modelling
of pottery data, in contrast to “traditional” datings which are
frequently based on an intuitive comparison of pottery assemblages.
Of course, this does not in any way deny the importance of
the pottery expert’s knowledge in the dating process; on the
contrary, it enables this knowledge to be better integrated through
a systemic procedure, which is essential when there is a very large
corpus of data to be analyzed.</p>
<p>The methodology is based on a statistical and visual approach
using two estimated density curves to date each archaeological
context. Two steps were required in the statistical procedure, each
leading to the construction of a density curve. The first enabled us
to estimate a date corresponding to the <i>terminus post quem</i> of the
context, a cursor reflecting an event dated in calendar time.
The second step, based on the results of the first,
allows the chronological profile of the context to be estimated,
giving a picture that is closer to archaeological time, in other words
the rate of accumulation.</p>
                      ")
               )
             )
    ),
    tabPanel(value = "2",span(icon("chart-area",lib="font-awesome"),strong("cerardat")),style="max-width:1200px;",
              tabsetPanel(id = "tabsetperso",type="tabs",
                      tabPanel(id="tab11","Overview",
                               sidebarPanel(width = 5,
                                            fluidRow(
                                              h4("Import your data:"),
                                              column(10,
                                                     actionButton("data",icon("import", lib = "glyphicon"))
                                              ),
                                              column(1,
                                                     # Help
                                                     dropdownButton(
                                                       h3("data format"),
                                                       hr(style="border-color: #222222;"),
                                                       HTML("<p>The data are in the form of two datasets.</p>
                                                       <p>1/ Reference data: count table plus a first column with dates, if known.</p>
                                                       <p>2/ Supplementary data: count table (same number of columns) plus a column for dates, if known.</p>
                                                       <p>Example:</p>
                                                       <img style='width:280px;' src='data_exemple.png' alt='data'>"),
                                                       circle = TRUE, status = "danger", icon = icon("question"), width = "300px",
                                                       size = "sm",right=FALSE,
                                                       tooltip = tooltipOptions(title = "Help", placement = "top")
                                                     )
                                                     # Help
                                              )
                                            ),



                                            div(style="height:20px;"),
                                            HTML("
<p>In settings:</p>
<ul>
  <li>You can check the results of the Correspondance Analysis</li>
  <li><b>You must choose the number of components</b> to keep in the linear model</li>
  <li>You can compare the predicted dates with the actual dates and define the confidence interval you want.</li>
</ul>
<p>And then you can see the results in the result and visualization tabs.</p>
                                                 ")
                               ),
                               mainPanel(style = 'max-width:1000px;',width = 7,
                                         HTML("
<p>For each archaeological context, we used pottery to construct
and analyze two different probability density curves of estimated
dates, one based on the dating of events in calendar time (hereafter
called “event time”), and one based on duration/time-span (hereafter
called “accumulation time”). The chain of reasoning for a given
archaeological context can be summarized as follows:</p>
<ul>
  <li>Obtain two density curves representing the dating of an
archaeological context.<ul>
    <li>Step 1: The first is a Gaussian curve derived from a linear
regression model. This allows the mean date of issue of
a coin to be estimated from the pottery profile. From an
archaeological point of view, this estimation of a terminus
post quem, with all the biases involved in the use of coins for
dating, provides a fundamental basis for chronology building
in calendar time;</li>
    <li>Step 2: The second curve is a mixture of Gaussians derived
from the previous model. The date of an archaeological
context is estimated using the weighted average of the
estimated dates of the pottery fabrics found within it.
Assuming that the fabrics are statistically independent, the
associated probability density of the date can be estimated
as a weighted sum of Gaussian distributions. From an
archaeological point of view, this dating estimation provides
a closer representation of accumulation time recorded in the
soil (Olivier, 2001; Wirtz and Olivier, 2003). At best,
depending on the quality of the archaeological context, it
can be interpreted as a formation process reflecting the
duration or succession of events on the scale of archaeological
time, and at worst, as imprecise dating due to
contamination of the context by residual or intrusive
material.</li>
  </ul></li>
  <li>Compare the density curves, in order to:<ul>
    <li>Validate the method from a chronological perspective,
exploring the duration or intensity of occupation for each
context;</li>
    <li>Identify the boundaries of socio-economic entities within
the broader area (supplemantary data), using a chronological model based on the reference site;</li>
    <li>Obtain a clearer understanding of chrono-functional issues
through a better interpretation of the type of archaeological context.</li>
  </ul></li>
</ul>

                                              ")
                                )
                      ),
             tabPanel(id="tab12","Settings",
                      tabsetPanel(id = "tabsetperso1",type="tabs",
                                  tabPanel(id="set1","Correspondance Analysis (CA)",
                                           sidebarPanel(
                                              HTML("
<p>A correspondence analysis (CA) was carried out to summarize
the information in the reference corpus data; the data matrix
consisted of the archaeological contexts and fabrics quantified
using MINVC whether or not they had
been dated with coins.</p>
<p>In order to estimate a date for the context, it is essential to refer
to objects dated by another source, in this instance, nonresidual
coins. These contexts were selected on a very strict
basis for their chrono-stratigraphic reliability, level of domestic
occupation, or enclosures with long urban stratigraphic
sequences, thereby minimizing any bias linked to disparity
between the date of the coin and that of the context.</p>
<p>estimating the date of an event recorded in
the ground (an archaeological context for the archaeologist) based
on the pottery assemblage of which it is comprised is achieved
by fitting a regression model linking a known date in calendar time,
in this instance the date of issue of a coin, to its pottery profile.</p>
<p>We used the results of the first step and the properties of the CA
to obtain an estimated value of the date of each fabric. We could
then define the archaeological time shown as dateAc, in other
words the accumulation time of a context, as the weighted sum of
fabric datings; the weights are the proportions of MINVC of each
fabric in the context.</p>

                                                   ")
                                           ),
                                           mainPanel(

                                             ####################################
                                             ##        Select Component        ##
                                             ####################################
                                             br(),br(),
                                             h3("Correspondence analysis"),
                                             hr(style="border-color: #222222;"),
                                             h4("Eigenvalues"),

                                             withSpinner(jqui_resizable(plotlyOutput("eigenvalue")),color='#18bc9c'),br(),

                                             numericInput("axe1", "Axe X", 1, min = 1, max = NA, step = 1),
                                             numericInput("axe2", "Axe Y", 2, min = 1, max = NA, step = 1),
                                             withSpinner(jqui_resizable(scatterD3Output("CA")),color='#18bc9c'),

                                             br()
                                           )
                                           ),
                                  tabPanel(id="set2","Linear model (lm)",
                                           sidebarPanel(
                                             verbatimTextOutput("cerar"),
                                             sliderInput("k", label = "Number of Components", min = 1, max = 20, step = 1, value = 2),
                                             HTML("
<p>We then kept only the first factorial axes, accounting for a part of the total
variance. In this way, becomes an incomplete
reconstitution of the data. This principle is used in many factor
analysis techniques, providing away of reducing the number of
explanatory variables in the linear regression model(the other components
may reasonably be ignored!) (Matzner-Løber, 2010).</p>
                                                  "),
                                             hr(style = "border-top: 1px solid #acacac;"),

                                             h4("statistical testing of residual hypotheses"),
                                             span("normality of residuals:"),
                                             verbatimTextOutput("SW"),
                                             span(textOutput("sw"), style="color:red"),
                                             span("no autocorrelation of residuals:"),
                                             verbatimTextOutput("DW"),
                                             span(textOutput("dw"), style="color:red"),
                                             span("heteroskedasticity of the variance of the residuals:"),
                                             verbatimTextOutput("BP"),
                                             span(textOutput("bp"), style="color:red")

                                           ),
                                           mainPanel(
                                             h4("Number of components to be kept in the linear model"),
                                             selectInput("quid", "display:",
                                                         c("both" = "both",
                                                           "PRESS" = "PRESS",
                                                           "MSE" = "MSE")),
                                             withSpinner(jqui_resizable(plotOutput("RMSE")),color='#18bc9c'),br(),
                                             withSpinner(jqui_resizable(plotOutput("R_sq")),color='#18bc9c'),br(),

                                             br(),
                                             uiOutput("formula_dateEv"),
                                             br(),
                                             verbatimTextOutput("lm"),
                                             uiOutput("lmplot")
                                           )
                                  ),
                                  tabPanel(id="set2","Control",
                                           sidebarPanel(
                                             numericInput("confidence", "Confidence:", 0.95, min = 0.00, max = 0.99, step= 0.01)#,
                                             #verbatimTextOutput("checkbox"),
                                             #checkboxGroupInput("col.refbox", "Choose which data are the reference:",
                                             #                   choiceNames = list(""),
                                             #                   choiceValues = list(0)
                                             #)
                                           ),
                                           mainPanel(
                                             h4("Comparison of estimated dates (black) with actual dates in ref data (red)"),
                                             HTML("<p>It is the dates that are used to build the model. It is advisable to choose reliable sets. You can see the effect of adding or removing components in the model on the date estimates.</p>"),
                                             withSpinner(jqui_resizable(plotlyOutput("evaldateref")),color='#18bc9c'),
                                             h4("Comparison of estimated dates (black) with actual dates in sup data (red)"),
                                             HTML("<p>These dates are not used in the construction of the model. It is advisable to add sets that are either far from the study area or whose pottery quality is not as reliable as the one above. These dates can also be used to check if the model is working.</p>"),
                                             withSpinner(jqui_resizable(plotlyOutput("evaldatesup")),color='#18bc9c'),

                                             br()
                                           )
                                  )
                      )

             ),
             tabPanel(id="tab13","Results",
                      tabsetPanel(id = "tabsetperso2",type="tabs",
                                  tabPanel(id="datesite","Dating archaeological contexts (dateEV and dateAC)",
                                           tableOutput('tablesite')
                                  ),
                                  tabPanel(id="dategt","Dating technical group",
                                           tableOutput('tableGT')
                                  )
                      )
             ),
             tabPanel(id="tab14","Visualization",
                      downloadButton("downloadData","Download plots as png"),br(),
                      tabsetPanel(id = "tabsetperso2",type="tabs",
                        tabPanel(id="col.ref","ref data",
                                 ####################################
                                 ##         cerardat plot          ##
                                 ####################################
                                 uiOutput("cerardat_plot_ref")
                                 ),
                        tabPanel(id="col.sub","sup data",
                                 uiOutput("cerardat_plot_sup")
                        )
                      )
             )
          )
    ),
    tabPanel(value = "3",span(icon("envelope", lib = "glyphicon"),strong("Contact")),style="max-width:1200px;",
      HTML("<style>
        .carte{
        border-left-width: 4px;
        border-left-style: solid;
        border-color: #2ac0a2;
        padding-left: 10px;
        }
        </style>
        <h3>Authors:</h3>
        <hr>
        <div class=\"carte\"><h4>L. Bellanger</h4><h5>mail: &lt;lise.bellanger@univ-nantes.fr&gt;</h5></div>
        <div class=\"carte\"><h4>P. Husi</h4><h5>mail: &lt;philippe.husi@univ-tours.fr&gt;</h5></div>
        <div class=\"carte\"><h4>A. Coulon</h4></div>
        <h3>Maintainer:</h3>
        <hr>
        <div class=\"carte\"><h4>A. Coulon</h4><h5>mail: &lt;arthur-coulon@outlook.fr&gt;</h5></div>
        ")
    )
  )
)
