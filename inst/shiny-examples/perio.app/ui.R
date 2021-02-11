library(SPARTAAS)
suppressMessages(library(shinyjs))
library(shinyjqui)
suppressMessages(library(shinydashboard))
library(colorspace)
suppressMessages(library(dplyr))
library(tidyr)
library(ggplot2)
library(shinyWidgets)
library(FactoMineR)
library(scatterD3)
library(dplyr)
library(fpc)
suppressMessages(library(plotly))
library(cluster) #for silhouette
library(stringr) #sub_extract

ui <- fluidPage(useShinyjs(),style="padding-top: 130px;",
    tags$style(HTML('#sort{color:white;background-color:#337ab7}')),
    tags$style(HTML('#data{color:white;background-color:#337ab7}')),
    tags$style(HTML('#subdivise{color:white;background-color:#337ab7}')),

    tags$style(HTML('#tri{color:white;background-color:#337ab7}')),
    tags$style(HTML('#sub1{color:white;background-color:#337ab7}')),
    tags$style(HTML('#sub2{color:white;background-color:#337ab7}')),
    tags$style(HTML('#import2{color:white;background-color:#337ab7}')),
    tags$style(HTML('#import3{color:white;background-color:#337ab7}')),
    tags$style(HTML('#Finishimport{color:white;background-color:#337ab7}')),

    tags$style(HTML('#back{color:white;background-color:#f0ad4e}')),
    tags$style(HTML('#back1{color:white;background-color:#337ab7}')),
    tags$style(HTML('#back2{color:white;background-color:#337ab7}')),
    tags$style(HTML('#tabsetperso{background-color:#e5e5e5;}')),
    tags$style(HTML('#help{text-align;right;}')),

    tags$head(
      tags$title("PerioApp!")
    ),
  #header panel
  absolutePanel(class = "panel panel-default",
    style="z-index: 2000;padding: 8px; border-bottom: 1px solid #CCC; background: #fff;opacity: 1;",
    top = 0, left = 0, right = 0,
    fixed = TRUE,
    h2("Perioclust: A chronological hierarchical agglomerative clustering method"),
    div(span(strong("SPARTAAS | Perioclust")),span("v0.9.19-5 ",style ="font-size:14px;"))
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
    background: #FFFFEE;width:30px;",
    bottom = 20, right = 0,
    fixed = TRUE,
    HTML("<a href=\"#top\"><span class=\"glyphicon glyphicon-chevron-up\"></span></a>")
  ),
  #-----------------------------------#
  #    NavBar
  #-----------------------------------#
  navbarPage("PerioClust",id="navbar",
    tabPanel(value = "1",span(icon("home", lib = "glyphicon"),strong("Home")),style="max-width:1200px;",
      imageOutput("spartaas")
    ),
    tabPanel(value = "2",span(icon("chart-area", lib = "font-awesome"),strong("Perioclust")),style="max-width:1200px;",
              tabsetPanel(id = "tabsetperso",type="tabs",
                      tabPanel(id="tab11","Overview",
                               sidebarPanel(style = 'max-width:300px;',
                                            h4("Import your data:"),actionButton("data",
                                                                                 icon("import", lib = "glyphicon")),
                                            h4("Refresh app:"),actionButton("refresh",
                                                                            icon("refresh", lib = "glyphicon"))
                               ),
                               mainPanel(style = 'max-width:1000px;',
                                         h4("Introduction"),
                                         hr(style="border-color: #222222;"),
                                         HTML("<p>PerioClust is a hierarchical ascending classification (HAC)
                                         algorithm under constraint.It has been developed to respond to problems
                                              related to chronology in Archaeology but its use can be extended to
                                              other fields (e.g. Ecology, Health).</p><p>This semi-supervised learning
                                              algorithm takes into account two types of information that can be subject
                                              to observational errors. The first source of information may be
                                              qualitative or quantitative, but the second must reflect the temporal
                                              structure of the data (e.g. stratigraphic information or time periods
                                              estimations). </p><p>The objective is to create a convex linear
                                              combination
                                              between these two sources of information. However, since errors have
                                              no reason to be of the same size, the choice of the mixing parameter
                                              (alpha) is important. This allows the two matrices to be merged in a
                                              balanced way.</p><br><br>"),
                                         h4("guide lines"),
                                         hr(style="border-color: #222222;"),
                                         HTML("<ul>
                                         <strong>overview tab (here)</strong>
                                         <li>Import your data <br>(you can skip this step and
                                         test perioclust with the default dataset)</li>
                                         <strong>Clustering settings tab</strong>
                                         <li>Choose the correct number of axes in the correspondance analysis</li>
                                         <li>Select the alpha</li>
                                         <li>Choose the number of groups <br>optimal numbers of groups button help
                                         you with
                                         two evaluations scores (e.g. WSS: Within Sum of Square and average
                                         silhouette score)</li>
                                         <strong>Results visualization tab</strong>
                                         <li>Sort the periods on the seriograph in temporal order
                                         (Sort periods button)</li>
                                         </ul><br><br>
                                              "),
                                         h4("First data source (pottery assemblages)"),
                                         hr(style="border-color: #222222;"),
                                         tableOutput("data1"),
                                         h4("Second data source (Stratigraphic network or timerange)"),
                                         hr(style="border-color: #222222;"),
                                         tableOutput("data2")
                                )
                      ),
             tabPanel(id="tab12","Clustering settings",
               sidebarPanel(style = 'max-width:300px;',
                            h4("Change number of axis in Correspondance Analysis:"),actionButton("axesCA", "Change"),

                            div(style="height:10px;"),
                            div(style="height:720px;"),
                            selectInput("method", label = "Aggregation method:",
                                        choices = list(
                                          "ward.D2" = 2,
                                          "complete" = 4,
                                          "single" = 3,
                                          "average (= UPGMA)" = 5,
                                          "mcquitty (= WPGMA)" = 6 #,"median (= WPGMC)" = 7,"centroid (= UPGMC)" = 8
                                        ), selected = 2),


                            actionButton("subdivise",span(icon("cut",lib="font-awesome"),"Subdivide")),
                            div(style="height:10px;"),
                            actionButton("reset",span(icon("undo",lib="font-awesome"),"Reset the dendrogram"))
                            ),
               mainPanel(style = 'max-width:1000px;',

                         ####################################
                         #         SelectAlpha PLOT         #
                         ####################################
                         br(),br(),
                         h3("Select Alpha"),
                         hr(style="border-color: #222222;"),
                         column(11,
                                dropdownButton(
                                  downloadButton("alpha.pdf","Save as pdf"),
                                  hr(style="border-color: #222222;"),
                                  downloadButton("alpha.png","save as png"),
                                  numericInput("width1", "width", 900, min = 0, max = NA, step = 1, width = NULL),
                                  numericInput("height1", "height", 500, min = 0, max = NA, step = 1, width = NULL),
                                  circle = TRUE, status = "primary",
                                  icon = icon("floppy-disk", lib = "glyphicon"), width = "150px",
                                  size = "sm",
                                  tooltip = tooltipOptions(title = "Download !")
                                )
                         ),
                         column(1,
                                # Help !
                                dropdownButton(
                                  h3("Select Alpha !"),
                                  hr(style="border-color: #222222;"),
                                  p("A criterion for choosing alpha IN [0;1] must be determined by balancing the
                                    weights between the two information sources in the final classification."),
                                  p("The CorCrit_alpha criterium represents the difference in absolute value
                                    between two cophenetic correlation (Cophenetic correlation is defined as
                                    the correlation between two distances matrices. It is calculated by considering
                                    the half distances matrices as vectors. It measures of how faithfully a dendrogram
                                    preserves the pairwise distances between the original unmodeled data points).
                                    The first correlation is associated with the comparison between D1 and
                                    ultrametric distances from the HAC with alpha fixed; while the second compares
                                    D2 and ultrametric distances from the HAC with alpha fixed. Then, in order
                                    to compromise between the information provided by D1 and D2, we decided to
                                    estimate alpha with hat(alpha) as the minimum."),
                                  circle = TRUE, status = "danger", icon = icon("question"), width = "500px",
                                  size = "sm",right=TRUE,
                                  tooltip = tooltipOptions(title = "Help !", placement = "top")
                                )
                                # Help !
                         ),br(),br(),

                         jqui_resizable(plotOutput("selectAlpha")),br(),

                         column(6,
                                sliderInput("alpha", label = "Alpha", min = 0, max = 1, step = 0.01, value = 0)
                         ),
                         column(6,
                                sliderInput("k",label ="Number of class",min = 2,max = 20,step = 1, value = 2)
                         ),
                         br(),
                         HTML("<p><img style=\"width:200px;\" src=\"GS/formule.PNG\"></p>"),
                         hr(style="border-color: #c2c2c2;"),
                         column(7,
                          h4("Result of correspondence analysis:"),
                          dropdownButton(
                            h3("Correspondence analysis"),
                            scatterD3Output("CA"),
                            circle = TRUE, status = "warning", icon = icon("info"), width = "600px", size = "sm",
                            right=TRUE
                          )
                         ),
                         column(5,
                           h4("Detail of resampling process:"),
                           dropdownButton(
                             h3("Resampling process"),
                             plotlyOutput("detail"),
                             circle = TRUE,status = "warning",icon = icon("info"),width = "800px",size = "sm",
                             right = TRUE
                           )
                         ),br(),br(),br(),br(),hr(style="border-color: #c2c2c2;"),br(),br(),


                         ####################################
                         #         Dendrogram PLOT          #
                         ####################################
                         h3("Dendrogram"),
                         hr(style="border-color: #222222;"),
                         column(11,
                                dropdownButton(
                                  downloadButton("dendrogram.pdf","Save as pdf"),
                                  hr(style="border-color: #222222;"),
                                  downloadButton("dendrogram.png","save as png"),
                                  numericInput("width2", "width", 900, min = 0, max = NA, step = 1, width = NULL),
                                  numericInput("height2", "height", 500, min = 0, max = NA, step = 1, width = NULL),
                                  circle = TRUE, status = "primary", icon = icon("floppy-disk", lib = "glyphicon"),
                                  width = "150px",
                                  size = "sm",
                                  tooltip = tooltipOptions(title = "Download !")
                                )
                         ),
                         column(1
                         ),br(),br(),
                         shinyjqui::jqui_resizable(plotOutput("dendrogramme")),br(),br(),
                         hr(style="border-color: #c2c2c2;"),
                         column(12,
                                h4("Optimal number of clusters:"),
                                dropdownButton(
                                  column(6,
                                         h3("Within Sum of Square !"),
                                         plotOutput("WSS"),
                                         column(6,
                                                downloadButton("wss.pdf","Save as pdf")
                                         ),
                                         column(6,
                                                downloadButton("wss.png","save as png")
                                         )
                                  ),
                                  column(6,
                                         h3("Average silhouette !"),
                                         plotOutput("silhouette"),
                                         column(6,
                                                downloadButton("silhouette.pdf","Save as pdf")
                                         ),
                                         column(6,
                                                downloadButton("silhouette.png","save as png")
                                         )
                                  ),
                                  circle = TRUE, status = "warning", icon = icon("info"),
                                  width = "800px", size = "sm",right = TRUE,
                                  tooltip = tooltipOptions(title = "Partition evaluation !")
                                )

                         ),br(),br(),br(),
                         hr(style="border-color: #c2c2c2;"),
                         br(),
                         column(11,
                                dropdownButton(
                                  downloadButton("silplot.pdf","Save as pdf"),
                                  hr(style="border-color: #222222;"),
                                  downloadButton("silplot.png","save as png"),
                                  numericInput("width4", "width", 900, min = 0, max = NA, step = 1, width = NULL),
                                  numericInput("height4", "height", 500, min = 0, max = NA, step = 1, width = NULL),
                                  circle = TRUE, status = "primary", icon = icon("floppy-disk",
                                                                      lib = "glyphicon"), width = "150px",
                                  size = "sm",
                                  tooltip = tooltipOptions(title = "Download !")
                                )
                         ),
                         column(1,
                                # Help !
                                dropdownButton(
                                  h3("Silhouette plot !"),
                                  hr(style="border-color: #222222;"),
                                  p("Silhouette refers to a method of interpretation and validation of
                                    consistency within clusters of data."),
                                  p("The silhouette value is a measure of how similar an object is to
                                    its own cluster (cohesion) compared to other clusters (separation).
                                    The silhouette ranges from −1 to +1, where a high value indicates
                                    that the object is well matched to its own cluster and poorly matched
                                    to neighboring clusters. If most objects have a high value, then the
                                    clustering configuration is appropriate. If many points have a low or
                                    negative value, then the clustering configuration may have too many or
                                    too few clusters. "),
                                  circle = TRUE, status = "danger", icon = icon("question"), width = "500px",
                                  size = "sm",right=TRUE,
                                  tooltip = tooltipOptions(title = "Help !", placement = "top")
                                )
                                # Help !
                         ),br(),br(),
                         jqui_resizable(plotOutput("SilPlot"))
                         )
             ),
             tabPanel(id="tab13","Results visualization",
                      div(id="nodata",
                          HTML("<br><p><strong>Warnings:</strong> You must go to the clustering settings
                               tab to set up the classification</p>")
                          ),
                      div(id="toggle",
                          sidebarPanel(style = 'max-width:300px;',
                                       div(style="height:200px;"),
                                       prettyToggle(
                                         inputId = "permute",
                                         label_on = "Seriation",
                                         label_off = "Seriation",
                                         icon_on = icon("check"),
                                         icon_off = icon("remove"),
                                         value = TRUE
                                       ),
                                       prettyToggle(
                                         inputId = "col_weight",
                                         label_on = "Weight color indicator",
                                         label_off = "Weight color indicator",
                                         icon_on = icon("check"),
                                         icon_off = icon("remove"),
                                         value = TRUE
                                       ),
                                       selectInput("show", label = "Show", choices = list("both" = 1,
                                                                                          "frequencies" = 2,
                                                                                          "EPPM" = 3),
                                                   selected = 1),
                                       actionButton("sort",span(icon("sort",lib="font-awesome"),"Sort periods"))
                          ),
                          mainPanel(style = 'max-width:1000px;',
                                    ####################################
                                    #         Seriograph PLOT          #
                                    ####################################
                                    br(),br(),
                                    h3("Seriograph"),
                                    hr(style="border-color: #222222;"),
                                    column(11,
                                           dropdownButton(
                                             downloadButton("seriograph.pdf","Save as pdf"),
                                             hr(style="border-color: #222222;"),
                                             downloadButton("seriograph.png","save as png"),
                                             numericInput("width3", "width", 900, min = 0, max = NA,
                                                          step = 1, width = NULL),
                                             numericInput("height3", "height", 500, min = 0, max = NA,
                                                          step = 1, width = NULL),
                                             circle = TRUE, status = "primary", icon = icon("floppy-disk",
                                                                    lib = "glyphicon"), width = "150px",
                                             size = "sm",
                                             tooltip = tooltipOptions(title = "Download !")
                                           )
                                    ),
                                    column(1,
                                           # Help !
                                           dropdownButton(
                                             h3("The seriograph (B. DESACHY) !"),
                                             hr(style="border-color: #222222;"),
                                             p("This tool makes it possible to highlight artisanal
                                               evolutions over time as well as to understand commercial
                                               relations thanks to imported ceramics. The percentages of
                                               each ceramic category are displayed. The percentages are
                                               calculated independently for each class. The percentage display
                                               allows you to compare the different classes but does not provide
                                               information on the differences in terms of enrolment. To fill
                                               this gap,
                                               the proportion of each class of their workforce is displayed on
                                               the seriograph (weight column). We have generalized this
                                               representation
                                               for other contingency data with Perioclust object."),
                                             circle = TRUE, status = "danger", icon = icon("question"),
                                             width = "500px",
                                             size = "sm",right=TRUE,
                                             tooltip = tooltipOptions(title = "Help !", placement = "top")
                                           )
                                           # Help !
                                    ),br(),br(),

                                    jqui_resizable(plotOutput("EPPM")),
                                    ####################################
                                    #          timerange PLOT          #
                                    ####################################
                                    conditionalPanel(
                                      condition = "1 == 0",
                                      checkboxInput("test","test",value = TRUE)
                                    ),
                                    conditionalPanel(
                                      condition = "input.datatype && input.test",
                                      br(),br(),
                                      h3("Timerange by cluster"),
                                      hr(style="border-color: #222222;"),
                                      column(11

                                      ),
                                      column(1,
                                             # Help !
                                             dropdownButton(
                                               h3("Timerange clust !"),
                                               hr(style="border-color: #222222;"),
                                               p("Each observation has a timerange estimation. We display the
                                                 timeranges grouping by cluster"),
                                               circle = TRUE, status = "danger", icon = icon("question"),
                                               width = "300px",
                                               size = "sm",right=TRUE,
                                               tooltip = tooltipOptions(title = "Help !", placement = "top")
                                             )
                                             # Help !
                                      ),br(),br(),
                                      jqui_resizable(plotlyOutput("timerange"))

                                    ),

                                    ####################################
                                    #         proportion PLOT          #
                                    ####################################
                                    br(),br(),
                                    h3("Top 4 (in proportion) for each cluster"),
                                    hr(style="border-color: #222222;"),

                                    jqui_resizable(plotlyOutput("prop_plot"))

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
        border-color: coral;
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
    ),
    #-----------------------------------#
    #    Wiki
    #-----------------------------------#
    tabPanel(value = "4",span(icon("question-sign", lib = "glyphicon"),strong("Wiki")),style="max-width:1200px;",
             HTML("<style>body {text-align: justify}</style><h1>Get started with the Perioclust application</h1>
                                      <h2>Table of content</h2>
                                      <div id=\"TOC\">
                                      <ul>
                                      <li><a href=\"#introduction\">Introduction</a>
                                        <ul>
                                          <li><a href=\"#perioclust\">Perioclust</a></li>
                                          <li><a href=\"#CAdist\">Construction of D1</a></li>
                                          <li><a href=\"#adjacency\">Construction of D2 for spatial data</a></li>
                                          <li><a href=\"#overlap\">Construction of D2 for temporal data</a></li>
                                          <li><a href=\"#corcrit\">Correlation Criterion</a></li>
                                          <li><a href=\"#resamp\">Resampling strategy</a></li>
                                        </ul>
                                      </li>
                                      <li><a href=\"#data1\">The inputs</a>
                                        <ul>
                                          <li><a href=\"#alpha\">Alpha</a></li>
                                          <li><a href=\"#k\">Number of class</a></li>
                                          <li><a href=\"#aggregation\">Aggregation method</a></li>
                                          <li><a href=\"#subdivide\">Subdivide</a></li>
                                          <li><a href=\"#seriation\">Seriation</a></li>
                                          <li><a href=\"#weight\">Weight color indicator</a></li>
                                          <li><a href=\"#show\">Vizualisation</a></li>
                                          <li><a href=\"#sort1\">Sort periods</a></li>
                                          <li><a href=\"#import\">Import your data</a></li>
                                          <li><a href=\"#csv\">CSV Format and write.table</a></li>
                                        </ul>
                                      </li>
                                      <li><a href=\"#outputs\">The outputs</a>
                                        <ul>
                                          <li><a href=\"#eval\">Evaluation Plot</a>
                                            <ul>
                                              <li><a href=\"#wssplot\">Within Sum of Square Plot</a></li>
                                              <li><a href=\"#avesilplot\">Average Silhouette Plot</a></li>
                                            </ul>
                                          </li>
                                          <li><a href=\"#selectalpha\">Select alpha</a></li>
                                          <li><a href=\"#dendrogram\">Dendrogram</a></li>
                                          <li><a href=\"#seriograph\">Seriograph</a></li>
                                          <li><a href=\"#timerangeclust\">Timerange.clust</a></li>
                                          <li><a href=\"#proportionplot\">Proportion plot</a></li>

                                          <li><a href=\"#CA\">Correspondences analysis plan 1-2</a></li>
                                          <li><a href=\"#resampling\">Resampling detail plot</a></li>

                                        </ul>
                                      </li>

                                      <li><a href=\"#references\">References</a></li>
                                      </ul>
                                      </div>

<hr>
<h2 id=\"introduction\">Introduction</h2>
<hr>
<p>Hierarchical bottom-up classification method with constraints. The method use two sources of informations.</p>

<h3 id=\"perioclust\">Perioclust</h3>
<p>The merging of the two data sources is done by a parameter (alpha) which allows to weight each source.</p>
<p><img style=\"width:200px;\" src=\"GS/formule.PNG\"></p>

<p>The first one is a contingency table. Rows must be individuals (archaeological site, ...)
and columns must be categories (type,...).

The second concerns network relationships that can be network type or time range.

In the network context the matrix contain only 0 or 1 (1 if there is not a relationship
and 0 if there is a relationship). In the time range context the matrix contains time overlap index.

The network object is a data frame with two columns. The first column contains the network elements
(same number as the number of lines of D1) and the second column,for temporal network, contains a
list of all the other elements connected to it. The list is a string composed of the names of the
elements separated by a comma.
</p><p><img src=\"GS/network.PNG\"></p><p>

The data frame for time range data contains the same first column. The second column contains the
lower temporal bound and the third column contains the upper temporal bound.</p>
<p><img src=\"GS/timerange.PNG\"></p>

<h3 id=\"CAdist\">Construction of D1</h3>
<p>We run an Correspondences analysis on a contingence table then return the distance matrix of
the coordinates in the plane 1-2.</p>
<p>You can choose the number of principal component to kept in the construction of D1.</p>
<p><img src=\"GS/ncp.PNG\"></p>
<p><img src=\"GS/ncpslider.PNG\"></p>

<h3 id=\"adjacency\">Construction of D2 for binary temporal data</h3>
<p>From a spatial network dataset we construct a adjacency matrix. Base on this matrix we generate
a dissmilarity matrix.
The matrix contain only 0 or 1, 1 if there is not a relationship and 0 if there is a relationship.</p>

<h3 id=\"overlap\">Construction of D2 for temporal data</h3>
<p>The overlap index is the ratio between internal overlap and total overlap over time. We define the
total overlap limit as follows: the minimum of the lower limits of the pair of individuals and the
maximum of the upper limits. We define the limit of the internal overlap as follows: the maximum of
the lower limits and the minimum of the upper limits.</p>
<p><img style=\"width:500px;\" src=\"GS/overlap.PNG\"></p>

<h3 id=\"corcrit\">Correlation Criterion</h3>
<p>A criterion for choosing alpha IN [0;1] must be determined by balancing the weights between the two
information sources in the final classification. To obtain alpha, we define the following criterion:
</p><p><img src=\"GS/formule2.PNG\"></p>
The CorCrit_alpha criterium in (1) represents the difference in absolute value between two cophenetic
correlation (Cophenetic correlation is defined as the correlation between two distances matrices.
It is calculated by considering the half distances matrices as vectors. It measures of how faithfully a
dendrogram preserves the pairwise distances between the original unmodeled data points).
The first correlation is associated with the comparison between D1 and ultrametric distances from the HAC
with alpha fixed; while the second compares D2 and ultrametric distances from the HAC with alpha fixed.
Then, in order to compromise between the information provided by D1 and D2, we decided to estimate alpha
with hat(alpha) such that:
</p><p><img src=\"GS/formule3.PNG\"></p>

<h3 id=\"resamp\">Resampling strategy</h3>
<p>To do this, a set of \"clones\" is created for each observation i. A clone c of observation i is a
copy of observation i for which the spatial (or temporal) relationships to others have been modified;
but it has the same values as observation i in N.</p>

<p>A HAC is then carried out. For each observation i of the real dataset, we generate a set of j clones
by checking all the possible (n-1) clones.</p>

<p>Intuitively, by varying α between 0 and 1, we will be able to identify when the clone and the initial
observation will be separated on the dendrogram. This moment will correspond to the value of α above which
the weight given to information on the connections between observations contained in D2 has too much impact
on the results compared to that of D1. In practice, to obtain the value α that minimizes the Correlation criterion,
we use a one-dimensional optimization procedure, combination of golden section search and successive parabolic
interpolation, implemented in the R function optimize.</p>

<hr>
<h2 id=\"data1\">The inputs</h2>
<hr>
<h3 id=\"alpha\">Alpha</h3>
<p>The slider input allow you to choose the value of the mixing parameter. You select the weight to give at the
first dissimilarity matrix (alpha) and the second (1-alpha). By default the value if estimate using the Correlation
Criterion define before.</p>
<p><img src=\"GS/alpha.PNG\"></p>

<h3 id=\"k\">Number of class</h3>
<p>The slider input defines the number of class you want to create. In order to help the decision you can look
the WSS plot and the silhouette information.</p>
<p><img src=\"GS/k.PNG\"></p>

<h3 id=\"aggregation\">Aggregation method</h3>
<p>The hierarchical bottom-up classification can use several method of aggregation during the construction of
the dendrogram. You can choose one of the methods except for the methods with risk of inversion
(if you know what you do, you can use package version for that).</p>
<p><img src=\"GS/aggregation.PNG\"></p>

<h4>Ward.D2</h2>
<p>Ward's method is a criterion applied in hierarchical cluster analysis. Ward's minimum variance method
is a special case of the objective function approach originally presented by Joe H. Ward, Jr.
Ward suggested a general agglomerative hierarchical clustering procedure, where the criterion for choosing
the pair of clusters to merge at each step is based on the optimal value of an objective function.</p>

<h4>complete</h2>
<p>At each step, the two clusters separated by the shortest distance are combined. In complete-linkage clustering,
the link between two clusters contains all element pairs, and the distance between clusters equals the distance
between those two elements (one in each cluster) that are farthest away from each other. The shortest of these
links that remains at any step causes the fusion of the two clusters whose elements are involved.</p>

<h4>single</h2>
<p>At each step combining two clusters that contain the closest pair of elements not yet belonging to the same
cluster as each other. the distance between two clusters is determined by a single element pair, namely those
two elements (one in each cluster) that are closest to each other. The shortest of these links that remains at
any step causes the fusion of the two clusters whose elements are involved.</p>

<h4>average (UPGMA)</h2>
<p>UPGMA (unweighted pair group method with arithmetic mean) is a simple agglomerative (bottom-up)
hierarchical clustering method. The method is generally attributed to Sokal and Michener.
The UPGMA method is similar to its weighted variant, the WPGMA method. Note that the unweighted term
indicates that all distances contribute equally to each average that is computed and does not refer to the math
by which it is achieved.
<br><br>
At each step, the nearest two clusters are combined into a higher-level cluster. The distance between any
two clusters A and B, each of size (i.e., cardinality) | A | and | B |, is taken to be the average of all
distances d ( x , y ) between pairs of objects x in A and y in B, that is, the mean distance between elements
of each cluster.</p>

<h4>mcquitty (WPGMA)</h2>
<p>WPGMA (Weighted Pair Group Method with Arithmetic Mean) is a simple agglomerative (bottom-up)
hierarchical clustering method, generally attributed to Sokal and Michener.
The WPGMA method is similar to its unweighted variant, the UPGMA method.
<br><br>
At each step, the nearest two clusters, say A and B, are combined into a higher-level cluster A ∪ B. Then,
its distance to another cluster C, d ( A U B , C ), is simply the arithmetic mean of the distances d ( A , C )
and d ( B , C ).</p>

<h3 id=\"subdivide\">Subdivide</h3>
<p>You can subdivide a class. clic the button select the class you want subdivide and select the number of subclasses.
You can reset the dendrogram with the reset button next to it</p>
<p><img src=\"GS/subdivide.PNG\"></p>

<h3 id=\"seriation\">Seriation</h3>
<p>This input allows you to activate or not the seriation of the columns. Matrix permutation uses an algorithm
called \"reciprocal averages\". Each line is assigned a rank ranging from 1 to n the number of lines.
A barycentre is calculated for each column by weighting according to the row rank. Finally,
the columns are reorganized by sorting them by their barycentre.</p>
<p><img src=\"GS/seriation.PNG\"></p>

<h3 id=\"weight\">Weight color indicator</h3>
<p>This input allows you to activate or not the coloration of the weight column in order to highlight
the confidence related to the quantity of data.</p>
<p><img src=\"GS/weight.PNG\"></p>

<h3 id=\"show\">Visualization</h3>
<p>This input let you choose which element to plot. There are tree options : plot
the Positive deviation from the average percentage (EPPM in French), plot the frequency or plot the both.
The average percentage is calculated for each category (columns) on the total number of accounts
(all classes combined).
From the average percentage we recover for each category and for each rows the difference between
the percentage of the
category in the class with the average percentage. The EPPM corresponds to the notion of independence deviation
(between rows and columns, between categories and time classes) in a chi-square test approach.
Although this approach is fundamental in statistical analysis, independence deviations are here purely indicative
and are not associated with a p_value that could determine the significance of deviations.</p>
<p><img src=\"GS/show.PNG\"></p>

<h3 id=\"sort1\">Sort periods</h3>
<p>The rows are initially in the order of appearance on the dendrogram. It must be possible to re-order
the classes in a temporal way. In the interface you can drag and drop the class (green) to change the order.
You can add a Hiatus between two periods and also remove a period by drag and drop this period in remove section.</p>
<p><img src=\"GS/sort.PNG\"></p>

<h3 id=\"import\">Import your data</h3>
<p>You can import your data. There are two step. You have to upload a first csv for the first source of information.
The data must be a contingency table.</p>
<p><img src=\"GS/import1.PNG\"></p>
<p>The second step consist on the upload of the second source. You can choose between two data type
(Network or Time range). The network object is a data frame with two columns. The first column contains
the elements (same number as the number of lines of the table of contingency) and the second column,
for temporal network (contains a list of all the other elements connected to it).
The list is a string composed of the names of the elements separated by a comma.
</p><p><img src=\"GS/network_ui.PNG\"></p><p>

The data frame for time range data contains the same first column.
The second column contains the lower temporal bound and the third column contains the upper temporal bound.
</p><p><img src=\"GS/timerange_ui.PNG\"></p><p>
<p><img src=\"GS/import3.PNG\"></p>
<p>The settings allow you to import diffrents data frame organization (Header, separator of column, ...).</p>
<h4 id=\"header\">Header</h4>
<p>Yes or not option. Do you have headers on your colunms ?</p>
<h4 id=\"rownames\">Rownames</h4>
<p>Yes or not option. Do you have rownames on your rows ?</p>
<h4 id=\"separator\">Separator</h4>
<p>Choose the character use to separate the colunms.</p>
<h4 id=\"quote\">Quote</h4>
<p>Choose the quote use to strings.</p>
<h4 id=\"dec\">Decimal</h4>
<p>Choose the character use to indicate decimal.</p>

<h4 id=\"csv\">CSV Format and write.table</h4>
<p>It is a data.frame with colunms separate by semicolon \";\".</p>
<p>The input format for importing data is the \".csv\" format but also supports the\".txt\" format as a csv file.</p>
<p>In R you can export your data frame into a csv file using write.csv2 or write.table.
In a csv you can choose a character to separate the columns.
In the same way, you can define the character to indicate the decimal point.</p>
<code>
write.table(data,file=\"path/to/name_file.csv\",sep=\";\",dec=\".\",row.names=FALSE,quote=FALSE)
</code>
<p>In Excel you can save as CSV format in order to import your data frame.</p>
<p>The import interface allows you to setup this values with the \"header\",
\"decimal\", \"separator\" and \"quote\" option.</p>


<hr>
<h2 id=\"output\">The outputs</h2>
<hr>

<h3 id=\"eval\">Evaluation Plot</h2>

<p>It is essential to be able to evaluate the different partitions of the hierarchy in order to identify
the one(s) that is (are) most relevant. The number k of classes of the partition to be retained is based
in our case on several indicators calculated for different values of K : the total within-class sum of square
(<code>WSS</code>) and the global average of the silhouette widths.</p>
<p>During the execution of the Perioclust method you have to make a choice. You have to cut the dendrogram.
This operation select the partition. In order to compare all the possibilities you can see the evaluation plot
(WSSPlot and AveSilPlot). This two plot evaluate the relative good quality of the partition.</p>
<h3 id=\"wssplot\">Within Sum of Square Plot (WSSPlot)</h3>
<p>This is the plot of within-groups sum of squares against number of clusters.
The Within Sum of Square decrease when the number of cluster increase.
In this plot the best partition is when add one or more clusters don’t decrease the WSS value.
It’s call the Elbow method.</p>
<h4>Example:</h4>
<p><img src=\"GS/WSS.PNG\"></p>
<p>On this graph we start by looking at the value for the lowest number of classes: 2.
if I add a third class we see that the WSS value will decrease (from 2.4 to 1.8). If I add another class
I will decrease this value again (from 1.8 to 1.4). We looking for the moment where ading a class is
therefore not interesting. In this case, we can keep a partition with 8 or 9 classes</p>
<h3 id=\"avesilplot\">Average silhouette Plot</h3>
<p>This graph shows the average silhouette width of each partition (ROUSSEEUW 1987).
The silhouette width is a limited index between -1 and 1, which is calculated for each observation.
The closer the value is to 1, the better the observation is classified. We look for the average value
for a partition closest to 1.</p>
<h4>Example:</h4>
<p><img src=\"GS/AveSil.PNG\"></p>
<p>On this graph we look for the maximum value. The best evaluation corresponds to the division into 9 classes.
Looking at the second best partition we identify the one with 10 classes.</p>

<h3 id=\"selectalpha\">Select alpha</h3>
<p>This plot show the Correlation criterion for each alpha. Thanks to resampling process we can add a 95% confidence
interval.</p>
<p><img src=\"GS/selectalpha.PNG\"></p>

<h3 id=\"dendrogram\">Dendrogram</h3>
<p>The dendrogram is the result of the classification process. You can see the classes and subclasses.</p>
<p><img src=\"GS/dendrogram.PNG\"></p>

<h3 id=\"seriograph\">Seriograph</h3>
<p>We chose the seriograph (B. DESACHY). This tool makes it possible to highlight artisanal
evolutions over time as well as to understand commercial relations thanks to imported ceramics.
The percentages of each ceramic category are displayed. The percentages are calculated independently for each class.
The percentage display allows you to compare the different classes but does not provide information on the
differences in terms of enrolment. To fill this gap, the proportion of each class of their workforce is displayed
on the seriograph (weight column).

We have generalized this representation for other contingency data with Perioclust object.</p>
<p><img src=\"GS/seriograph.PNG\"></p>

<h3 id=\"timerangeclust\">Timerange.clust</h3>
Show the diffrent time range of observations for each cluster. (you have to import data whith time range)
<p><img src=\"GS/timerangeclust.PNG\"></p>

<h3 id=\"proportionplot\">Proportion plot</h3>
Show the top 4 of the variable for each cluster. (you have to import data whith time range)
<p><img src=\"GS/prop_plot.PNG\"></p>

<h3 id=\"CA\">Correspondences analysis plan 1-2</h3>
<p>The plane 1-2 of the correspondences analysis use in the construction of the first dissimilarity matrix D1</p>
<p><img src=\"GS/CA.PNG\"></p>

<h3 id=\"resampling\">Resampling detail plot</h3>
<p>Base on a re-sampling process, we generate clone and we check for which alpha the clone and the original
object are separeted on the dendrogram. The function show each set of clone curve.</p>
<p><img src=\"GS/resampling.PNG\"></p>


<hr>
<h2 id=\"references\">References</h2>
<hr>
<p><br>
  Rousseeuw, P.J. (1987). “Silhouettes: a graphical aid to the interpretation and validation of cluster analysis”.
  In : J. Comput. Appl. Math. 20, 53–65. <a href=\"DOI:https://doi.org/10.1016/0377-0427(87)90125-7\"
  class=\"uri\">DOI:https://doi.org/10.1016/0377-0427(87)90125-7</a>.
  <br><br>Ward, J. H., Jr. (1963). \"Hierarchical Grouping to Optimize an Objective Function\".
  Journal of the American Statistical Association, 58, 236–244.
  <br><br>Sokal R. and Michener C. (1958). \"A statistical method for evaluating systematic
  relationships\". University of Kansas Science Bulletin. 38: 1409–1438.
  <br><br>Desachy B. (2004). \"Le sériographe EPPM : un outil informatisé de sériation graphique pour les tableaux
  de comptages\". Revue Archéologique de picardie. 3/4 39-56.
</p>
<br><br><br><br>

")
    )
  )
)
