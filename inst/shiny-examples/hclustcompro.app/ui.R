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
library(shinythemes)

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
    ')),
    tags$head(
      tags$link(rel = "icon", type = "image/gif", href = "https://spartaas.gitpages.huma-num.fr/r-package/img/lambda.png"),
      tags$title("hclustcompro")
    ),
  #header panel
  absolutePanel(class = "panel panel-default",
    style="z-index: 2000;padding: 8px; background: #ecf0f1; opacity: 1;border-bottom: 1px solid #2c3e50;",
    top = 0, left = 0, right = 0,
    fixed = TRUE,
    h2("hclustcompro: Hierarchical Agglomerative Clustering method by compromise"),
    div(span(strong("SPARTAAS | hclustcompro")))
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
  #    NavBar
  #-----------------------------------#
  navbarPage("hclustcompro",id="navbar",
    tabPanel(value = "1",span(icon("home", lib = "glyphicon"),strong("Home")),style="max-width:1200px;",
             br(),
             sidebarLayout(
               #-----------------------------------#
               sidebarPanel(style="",
                            h1("R Package:"),
                            br(),
                            HTML("<p>This method is part of the <a target='_blank' href='https://spartaas.gitpages.huma-num.fr/r-package/index.html'>SPARTAAS</a> package.</p>
                                                        <p>If you are interested you can install our R package available on the <a href='https://cran.r-project.org/package=SPARTAAS'>CRAN</a> and on <a href='https://github.com/arliph/SPARTAAS'>GitHub</a>.</p>

                                                        ")

               ),
               #-----------------------------------#
               mainPanel(
                 HTML("
                                              <h3>Introduction</h3>
                                              <p>hclustcompro is a hierarchical clustering that takes into account two sources of information associated with the same objects. This method, called CAH by compromise (hclustcompro), allows a compromise between the hierarchies obtained from each source considered separately. A convex combination of the dissimilarities associated with each of the sources is used to modify the dissimilarity measure in the classical CAH algorithm. The choice of the mixing parameter is the key point of the method. We propose an objective function to be minimised based on the absolute difference of the correlations between the initial dissimilarities and the kinetic distances, as well as a resampling procedure to ensure the robustness of the choice of the mixing parameter.</p>
                                              ")
               )
             )
    ),
    tabPanel(value = "2",span(icon("chart-area", lib = "font-awesome"),strong("hclustcompro")),style="max-width:1200px;",
              tabsetPanel(id = "tabsetperso",type="tabs",
                      tabPanel(id="tab11","Overview",
                               sidebarPanel(width = 5,
                                            h4("Import your data:"),actionButton("data",
                                                                                 icon("import", lib = "glyphicon")),
                                            h4("Refresh app:"),actionButton("refresh",
                                                                            icon("refresh", lib = "glyphicon")),
                                            hr(style="border-color: #222222;"),
                                            h4("guide lines"),
                                            HTML("<p>You have three subtabs to run hclustcompro:</p>"),
                                            actionLink(inputId = "over", label = "Overview (here)"),
                                            HTML("<ul>
                                         <li>Import your data <br>(you can skip this step and
                                         test hclustcompro with the default dataset)</li>
                                         </ul>"),
                                            actionLink(inputId = "clustsett", label = "Clustering settings"),
                                          HTML("<ul>
                                         <li>Select the correct number of axes in the correspondence analysis</li>
                                         <li>Select the alpha value</li>
                                         <li>Select the number of groups <br>The optimal number of groups button helps you with two evaluation scores (e.g. WSS: within sum of squares and average silhouette score)</li>
                                         </ul>"),
                                          actionLink(inputId = "resviz", label = "Results visualization"),
                                          HTML("<ul>
                                         <li>Sort the clusters on the seriograph in chronological order (Sort clusters button)</li>
                                         </ul><br><br>
                                              ")

                               ),
                               mainPanel(style = 'max-width:1000px;',width = 7,
                                         h4("The data"),
                                         HTML("<p>hclustcompro uses two types of input data. The first is a contingency table (counting ceramics or other furniture). The second can be either stratigraphic relationships or dates (time range) for the same individuals as the first source.</p>
                                              <p>See example below.</p>"),
                                         actionLink(inputId = "wiki", label = "Go to wiki for more detail"),
                                         h4("First data source (contingency table)"),
                                         hr(style="border-color: #222222;"),
                                         tableOutput("data1f"),
                                         h4("Second data source (Stratigraphic connection or Time range data)"),
                                         hr(style="border-color: #222222;"),
                                         h5("Stratigraphic connection:"),
                                         tableOutput("data2f"),
                                         h5("Timerange data:"),
                                         tableOutput("data3f")
                                )
                      ),
             tabPanel(id="tab12","Clustering settings",
               sidebarPanel(
                            h4("Change number of axes in correspondence analysis:"),actionButton("axesCA", "Change"),
                            hr(style="border-color: #c2c2c2;"),
                            column(7,
                                   h4("Result of correspondence analysis:"),
                                   dropdownButton(
                                     h3("Correspondence analysis"),
                                     numericInput("axe1", "Axe X", 1, min = 1, max = NA, step = 1),
                                     numericInput("axe2", "Axe Y", 2, min = 1, max = NA, step = 1),
                                     scatterD3Output("CA"),
                                     icon = icon("info"), width = "600px", size = "sm",
                                     right=FALSE
                                   )
                            ),
                            div(style="height:1260px;"),
                            selectInput("method", label = "Aggregation method:",
                                        choices = list(
                                          "ward.D2" = 2,
                                          "complete" = 4,
                                          "single" = 3,
                                          "average (= UPGMA)" = 5,
                                          "mcquitty (= WPGMA)" = 6 #,"median (= WPGMC)" = 7,"centroid (= UPGMC)" = 8
                                        ), selected = 2),




                            column(7,
                                   actionButton("subdivise",span(icon("scissors",lib="font-awesome"),"Subdivide"))
                            ),
                            column(5,
                                   actionButton("reset",span(icon("rotate-left",lib="font-awesome"),"Reset"))
                            ),div(style="height:40px;")


                            ),
               mainPanel(

                         ####################################
                         ##        SelectAlpha PLOT        ##
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
                                  tooltip = tooltipOptions(title = "Download")
                                )
                         ),
                         column(1,
                                # Help
                                dropdownButton(
                                  h3("Select Alpha"),
                                  hr(style="border-color: #222222;"),
                                  p("A criterion for the choice of alpha IN [0;1] must be determined by balancing the weights between the two sources of information in the final classification."),
                                  p("The CorCrit_alpha criterion represents the difference in absolute value between two cophenetic correlations (cophenetic correlation is defined as the correlation between two distance matrices. It is calculated by considering the half distance matrices as vectors. It measures how faithfully a dendrogram preserves the pairwise distances between the original unmodeled data points). The first correlation is related to the comparison between D1 and ultrametric distances from the clustering with alpha fixed, while the second compares D2 and ultrametric distances from the clustering with alpha fixed. As a compromise between the information provided by D1 and D2, we then decided to estimate alpha with hat(alpha) as the minimum."),
                                  circle = TRUE, status = "danger", icon = icon("question"), width = "500px",
                                  size = "sm",right=TRUE,
                                  tooltip = tooltipOptions(title = "Help", placement = "top")
                                )
                                # Help
                         ),br(),br(),

                         jqui_resizable(plotOutput("selectAlpha")),br(),


                         sliderInput("alpha", label = "Alpha", min = 0, max = 1, step = 0.01, value = 0),


                         br(),
                         uiOutput("formula"),
                         #HTML("<p><img style=\"width:200px;\" src=\"GS/formule.png\"></p>"),


                         ####################################
                         #            NB GROUPE             #
                         ####################################
                         br(),br(),
                         h3("Select the number of groups in the partition"),
                         hr(style="border-color: #222222;"),
                         column(width=6,
                                column(10,
                                       dropdownButton(
                                         downloadButton("silhouette.pdf","Save as pdf"),
                                         hr(style="border-color: #222222;"),
                                         downloadButton("silhouette.png","save as png"),
                                         circle = TRUE, status = "primary",
                                         icon = icon("floppy-disk", lib = "glyphicon"), width = "150px",
                                         size = "sm",
                                         tooltip = tooltipOptions(title = "Download")
                                       ),
                                ),
                                column(1,
                                       # Help
                                       dropdownButton(
                                         h3("Average silhouette widths"),
                                         hr(style="border-color: #222222;"),
                                         p("This graph shows the average silhouette width of each partition (Rousseeuw 1987). The silhouette width is a bounded index between -1 and 1, calculated for each observation. The closer the value is to 1, the better the observation is classified. We look for the average value for a partition that is closest to 1."),
                                         circle = TRUE, status = "danger", icon = icon("question"), width = "500px",
                                         size = "sm",right=TRUE,
                                         tooltip = tooltipOptions(title = "Help", placement = "top")
                                       )
                                       # Help
                                ),

                                plotOutput("silhouette")
                         ),

                         column(width=6,
                                column(10,
                                      dropdownButton(
                                        downloadButton("wss.pdf","Save as pdf"),
                                        hr(style="border-color: #222222;"),
                                        downloadButton("wss.png","save as png"),
                                        circle = TRUE, status = "primary",
                                        icon = icon("floppy-disk", lib = "glyphicon"), width = "150px",
                                        size = "sm",
                                        tooltip = tooltipOptions(title = "Download")
                                      )
                                ),
                                column(1,
                                       # Help
                                       dropdownButton(
                                         h3("WSS"),
                                         hr(style="border-color: #222222;"),
                                         p("This is the plot of the within group sum of squares against the number of clusters. The within group sum of squares decreases as the number of clusters increases. In this plot, the best partition is when adding one or more clusters doesn't decrease the WSS value. It's called the elbow method."),
                                         circle = TRUE, status = "danger", icon = icon("question"), width = "500px",
                                         size = "sm",right=TRUE,
                                         tooltip = tooltipOptions(title = "Help", placement = "top")
                                       )
                                       # Help
                                ),
                                plotOutput("WSS")
                         ),div(style="height:500px;"),

                         sliderInput("k",label ="Number of class",min = 2,max = 20,step = 1, value = 2),


                         ####################################
                         #         Dendrogram PLOT          #
                         ####################################
                         h3("Dendrogram"),
                         hr(style="border-color: #222222;"),
                         verbatimTextOutput("click_info"),
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
                                  tooltip = tooltipOptions(title = "Download")
                                )
                         ),
                         column(1
                         ),br(),br(),
                         shinyjqui::jqui_resizable(plotOutput("dendrogramme")),br(),br(),



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
                                  tooltip = tooltipOptions(title = "Download")
                                )
                         ),
                         column(1,
                                # Help
                                dropdownButton(
                                  h3("Silhouette plot"),
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
                                  tooltip = tooltipOptions(title = "Help", placement = "top")
                                )
                                # Help
                         ),br(),br(),
                         jqui_resizable(plotOutput("SilPlot"))
                         )
             ),
             tabPanel(id="tab13","Results visualization",
                      div(id="nodata",
                          HTML("<br><p><strong>Warnings:</strong> You must go to the clustering settings tab to set up the clustering first.</p>")
                          ),
                      div(id="toggle",
                          sidebarPanel(style = 'max-width:300px;',
                                       div(style="height:200px;"),
                                       prettyToggle(
                                         inputId = "permute",
                                         label_on = "Seriation",
                                         label_off = "Seriation",
                                         icon_on = icon("check"),
                                         icon_off = icon("trash"),
                                         value = TRUE
                                       ),
                                       selectInput("show", label = "Show", choices = list("both" = 1,
                                                                                          "frequencies" = 2,
                                                                                          "EPPM" = 3),
                                                   selected = 1),
                                       actionButton("sort",span(icon("sort",lib="font-awesome"),"Sort cluster"))
                          ),
                          mainPanel(style = 'max-width:1000px;',
                                    ####################################
                                    ##        Seriograph PLOT         ##
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
                                             tooltip = tooltipOptions(title = "Download")
                                           )
                                    ),
                                    column(1,
                                           # Help
                                           dropdownButton(
                                             h3("The seriograph (Desachy 2004)"),
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
                                               for other contingency data with hclustcompro object."),
                                             circle = TRUE, status = "danger", icon = icon("question"),
                                             width = "500px",
                                             size = "sm",right=TRUE,
                                             tooltip = tooltipOptions(title = "Help", placement = "top")
                                           )
                                           # Help
                                    ),br(),br(),

                                    jqui_resizable(plotOutput("EPPM")),
                                    downloadButton("table_serio.csv", "Download the contingency table"),
                                    div(style="height:20px;"),
                                    ####################################
                                    ##         timerange PLOT         ##
                                    ####################################
                                    conditionalPanel(
                                      condition = "1 == 0",
                                      checkboxInput("test","test",value = TRUE)
                                    ),
                                    conditionalPanel(
                                      condition = "input.datatype == 2 && input.test",
                                      br(),br(),
                                      h3("Timerange plot"),
                                      hr(style="border-color: #222222;"),
                                      column(11

                                      ),
                                      column(1,
                                             # Help
                                             dropdownButton(
                                               h3("Timerange clust"),
                                               hr(style="border-color: #222222;"),
                                               p("Each observation has a timerange estimation. We display the
                                                 timeranges grouping by cluster"),
                                               circle = TRUE, status = "danger", icon = icon("question"),
                                               width = "300px",
                                               size = "sm",right=TRUE,
                                               tooltip = tooltipOptions(title = "Help", placement = "top")
                                             )
                                             # Help
                                      ),br(),br(),
                                      jqui_resizable(plotlyOutput("timerange")),
                                      div(style="height:20px;"),

                                    )

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
    ),
    #-----------------------------------#
    #    Wiki
    #-----------------------------------#
    tabPanel(value = "4",span(icon("question-sign", lib = "glyphicon"),strong("Wiki")),style="max-width:1200px;",
             HTML("<style>body {text-align: justify}</style><h1>Get started with the hclustcompro application</h1>
                                      <h2>Table of content</h2>
                                      <div id=\"TOC\">
                                      <ul>
                                      <li><a href=\"#introduction\">Introduction</a>
                                        <ul>
                                          <li><a href=\"#hclustcompro\">hclustcompro</a></li>
                                          <li><a href=\"#CAdist\">Construction of D1</a></li>
                                          <li><a href=\"#adjacency\">Construction of D2 for spatial data</a></li>
                                          <li><a href=\"#overlap\">Construction of D2 for temporal data</a></li>
                                          <li><a href=\"#corcrit\">Correlation Criterion</a></li>
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
                                          <li><a href=\"#show1\">Vizualisation</a></li>
                                          <li><a href=\"#sort1\">Sort clusters</a></li>
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
                                          <li><a href=\"#timerangeclust\">Timerange</a></li>

                                          <li><a href=\"#CAp\">Correspondences analysis</a></li>

                                        </ul>
                                      </li>

                                      <li><a href=\"#references\">References</a></li>
                                      </ul>
                                      </div>

<hr>
<h2 id=\"introduction\">Introduction</h2>
<hr>
<p>Hierarchical clustering method by compromise. The method uses two sources of information.</p>

<h3 id=\"hclustcompro\">hclustcompro</h3>
<p>The merging of the two data sources is done by a parameter (alpha) that allows to weight each source.</p>
<p>$$ \\boldsymbol{D}_\\alpha=\\alpha \\boldsymbol{D}_1+(1-\\alpha) \\boldsymbol{D}_2 $$</p>

<p>The first is a contingency table. The rows must be individuals (site, ...) and the columns must be categories (type, ...). The second concerns the same individuals and can be stratigraphic connection type or time range data. The stratigraphic link object is a data frame with two columns. The first column contains the network elements (the same number as the number of lines in D1) and the second column contains a list of all the other elements connected to it. The list is a string of element names separated by commas.
</p><p><img src=\"GS/network.png\"></p>
<p>The data frame for time range data contains the same first column. The second column contains the lower time limit and the third column contains the upper time limit.</p>
<p><img src=\"GS/timerange.png\"></p>

<h3 id=\"CAdist\">Construction of D1</h3>
<p>Our first source of information is a contingency table. We need to construct a distance matrix from this table. We perform a correspondence analysis (CA) on the contingency table and then use the distances of the components (chi-square metric).</p>

<p>You must select the number of axes (CA components) to be used to construct the distance matrix.</p>

<p>By examining the eigenvalues, it is possible to determine the number of principal axes to be considered. The eigenvalues correspond to the amount of information retained by each axis.</p>

<p>A heuristic method is often used: one constructs the graph of eigenvalues sorted in descending order, retaining the eigenvalues (and associated principal components) that precede the first 'elbow'.</p>

<p>Another approach is based on maintaining an approximation quality of the original data, measured by the explained inertia rate (e.g. 85%). As many axes as necessary are chosen so that the sum of the corresponding eigenvalues exceeds the target inertia rate.</p>

<p><img src=\"GS/ncp.png\"></p>
<p><img src=\"GS/ncpslider.png\"></p>

<h3 id=\"adjacency\">Construction of D2 for stratigraphic connection</h3>
<p>From this data set, we construct an adjacency matrix. Based on this matrix, we generate a dissimilarity matrix. The matrix contains only 0 or 1, 1 if there is no relationship and 0 if there is a relationship.</p>


<h3 id=\"overlap\">Construction of D2 for time range data</h3>
<p>The overlap index is the ratio between the overlap or separation and the cumulative extend of both individuals. We define the cumulative extend of both as: the minimum of the lower limits of the pair of individuals and the maximum of the upper limits. We define the overlap or separation as: the maximum of the lower limits and the minimum of the upper limits.</p>
<p>From this ratio, we construct a dissimilarity by setting the value between 0 and 1, ensuring that the closer the value is to 1, the further apart the time ranges are.</p>
<p><img style=\"width:500px;\" src=\"GS/overlap.png\"></p>

<h3 id=\"corcrit\">Correlation Criterion</h3>
<p>A criterion for choosing alpha IN [0;1] must be determined by balancing the weights between the two
sources of information in the final classification. To obtain alpha, we define the following criterion:
</p><p>$$\\operatorname{CorCrit}_\\alpha=\\left|\\operatorname{Cor}\\left(\\mathbf{D}_\\alpha^{\\text {coph }}, \\mathbf{D}_1\\right)-\\operatorname{Cor}\\left(\\mathbf{D}_\\alpha^{\\text {coph }}, \\mathbf{D}_2\\right)\\right|
$$</p>
The CorCrit_alpha criterion in (1) represents the difference in absolute value between two cophenetic correlations (cophenetic correlation is defined as the correlation between two distance matrices. It is calculated by considering the half distance matrices as vectors. It measures how faithfully a dendrogram preserves the pairwise distances between the original unmodeled data points). The first correlation is related to the comparison between D1 and ultrametric distances from the clustering with alpha fixed, while the second compares D2 and ultrametric distances from the clustering with alpha fixed. Then, in order to compromise between the information provided by D1 and D2, we decided to estimate alpha with hat(alpha) such that:</p>
<p>$$
\\hat{\\alpha}=\\min _\\alpha \\text { CorCrit }{ }_\\alpha
$$</p>



<hr>
<h2 id=\"data1\">The inputs</h2>
<hr>
<h3 id=\"alpha\">Alpha</h3>
<p>The slider input allows you to select the value of the mixing parameter. You select the weight to be given to the first dissimilarity matrix (alpha) and the second (1-alpha). By default, the value is the one previously defined when estimating using the correlation criterion.</p>
<p><img src=\"GS/alpha.png\"></p>

<h3 id=\"k\">Number of class</h3>
<p>The slider input defines the number of clusters you want. To help you decide, you can view the WSS plot and silhouette information.</p>
<p><img src=\"GS/k.png\"></p>

<h3 id=\"aggregation\">Aggregation method</h3>
<p>The hierarchical clustering can use several aggregation methods when constructing the dendrogram. You can choose any of the methods, except the ones with the risk of inversion (if you know what you are doing, you can use the package version for this).</p>
<p><img src=\"GS/aggregation.png\"></p>

<h4>Ward.D2</h2>
<p>Ward's method is a criterion used in hierarchical cluster analysis. Ward's minimum variance method is a special case of the objective function approach originally proposed by Joe H. Ward, Jr. Ward proposed a general agglomerative hierarchical clustering procedure where the criterion for selecting the pair of clusters to merge at each step is based on the optimal value of an objective function.</p>

<h4>complete</h2>
<p>At each step, the two clusters separated by the shortest distance are combined. In complete-linkage clustering, the link between two clusters contains all element pairs, and the distance between clusters is equal to the distance between the two elements (one in each cluster) that are farthest apart. The shortest of these links remaining at each step causes the two clusters whose elements are involved to merge.</p>

<h4>single</h2>
<p>The distance between two clusters is determined by a single pair of elements, namely the two elements (one in each cluster) that are closest to each other. The shortest of these links remaining at any step causes the two clusters whose elements are involved to merge.</p>

<h4>average (UPGMA)</h2>
<p>UPGMA (unweighted pair group method with arithmetic mean) is a simple agglomerative (bottom-up) hierarchical clustering method. It is generally credited to Sokal and Michener. The UPGMA method is similar to its weighted variant, the WPGMA method. Note that the term 'unweighted' indicates that all distances contribute equally to each mean computed, and does not refer to the mathematics by which it is obtained.
<br><br>
At each step, the two closest clusters are combined to form a higher level cluster. The distance between any two clusters A and B, of size (i.e. cardinality) | A | and | B | respectively, is taken to be the average of all distances d ( x , y ) between pairs of objects x in A and y in B, i.e. the average distance between the elements of each cluster.</p>

<h4>mcquitty (WPGMA)</h2>
<p>WPGMA (Weighted Pair Group Method with Arithmetic Mean) is a simple agglomerative (bottom-up)
hierarchical clustering method, generally attributed to Sokal and Michener.
The WPGMA method is similar to its unweighted variant, the UPGMA method.
<br><br>
At each step, the closest two clusters, say A and B, are combined into a higher-level cluster A ∪ B. Then,
its distance to another cluster C, d ( A U B , C ), is simply the arithmetic mean of the distances d ( A , C )
and d ( B , C ).</p>

<h3 id=\"subdivide\">Subdivide</h3>
<p>You can subdivide a cluster by clicking on the button, selecting the cluster you want to subdivide and selecting the number of subclusters. You can reset the dendrogram using the reset button next to the dendrogram.</p>
<p><img src=\"GS/subdivide.png\"></p>

<h3 id=\"seriation\">Seriation</h3>
<p>This input allows you to enable or disable the seriation of the columns. The matrix permutation uses an algorithm called 'reciprocal averages'. Each row is assigned a rank from 1 to n, the number of rows. For each column, a barycentre is calculated by weighting according to the row rank. Finally, the columns are reorganised by sorting them according to their barycentre.</p>
<p><img src=\"GS/seriation.png\"></p>


<h3 id=\"show1\">Visualization</h3>
<p>This input allows you to select the element to be plotted. There are tree options: plot the positive deviation from the average percentage (EPPM in French), plot the frequency or plot both. The average percentage is calculated for each category (columns) on the total number of accounts (all classes combined). From the average percentage, we obtain for each category and for each row the difference between the percentage of the category in the class and the average percentage. The EPPM corresponds to the notion of independence deviation (between rows and columns, between categories and time classes) in a chi-square test approach. Although this approach is fundamental in statistical analysis, the independence deviations here are purely indicative and are not associated with a p_value that could determine the significance of the deviations.</p>
<p><img src=\"GS/show.png\"></p>

<h3 id=\"sort1\">Sort clusters</h3>
<p>The rows are initially in the order in which they appear on the dendrogram. It is possible to rearrange the classes in time. In the interface you can drag and drop the class (green) to change the order. You can add a gap between two clusters and also remove a cluster by dragging and dropping it in the remove section.</p>
<p><img src=\"GS/sort.png\"></p>

<h3 id=\"import\">Import your data</h3>
<p>You can import your data. There are two steps. You need to upload a first csv file for the first source of information.
The data must be a contingency table.</p>
<p><img src=\"GS/import1.png\"></p>
<p>The second step is to upload the second source. You can choose between two types of data (stratigraphic link or time range data). The stratigraphic link object is a data frame with two columns. The first column contains the elements (the same number as the number of rows in the contingency table) and the second column contains a list of all the other elements connected to it. The list is a string consisting of the names (those in the first column) of the elements separated by a comma.</p>
<p><img src=\"GS/network_ui.png\"></p><p>The data frame for time range data contains the same first column. The second column contains the lower time limit and the third column contains the upper time limit.</p>
<p><img src=\"GS/timerange_ui.png\"></p>

<p><img src=\"GS/import3.png\"></p>
<p>The settings allow you to import different data frame organisation (header, column separator, ...).</p>
<h4 id=\"header\">Header</h4>
<p>Yes or no option. Do you have headers on your columns?</p>
<h4 id=\"rownames\">Rownames</h4>
<p>Yes or no option. Do you have row names on your rows?</p>
<h4 id=\"separator\">Separator</h4>
<p>Select the character you want to use to separate the columns.</p>
<h4 id=\"quote\">Quote</h4>
<p>Select the quotation marks to use on strings.</p>
<h4 id=\"dec\">Decimal</h4>
<p>Select the character to use to indicate the decimal point.</p>

<h4 id=\"csv\">CSV Format and write.table</h4>
<p>It is a data.frame with columns separated by semicolons \";\".</p>
<p>The input format for importing data is the .csv format, but also supports the .txt format as a .csv file.</p>
<p>In R, you can export your data frame to a csv file using write.csv2 or write.table. In a csv you can choose a character to separate the columns. In the same way, you can define the character to indicate the decimal point.</p>
<code>
write.table(data,file=\"path/to/name_file.csv\",sep=\";\",dec=\".\",row.names=FALSE,quote=FALSE)
</code>
<p>In Excel you can save in csv format in order to import your data frame.</p>
<p>The import interface allows you to set these values using the 'header', 'decimal', 'separator' and 'quote' options.</p>


<hr>
<h2 id=\"output\">The outputs</h2>
<hr>

<h3 id=\"eval\">Evaluation Plot</h2>

<p>It is essential to be able to evaluate the different partitions of the hierarchy in order to identify the most relevant one(s). In our case, the number k of clusters of the partition to be retained is based on several indicators calculated for different values of K: the total within-class sum of squares (WSS) and the global average of the silhouette widths.</p>
<p>When running the hclustcompro method, you have to make a choice. You must cut the dendrogram. This operation selects the partition. To compare all the possibilities, you can see the evaluation plots (WSSPlot and AveSilPlot). These two plots evaluate the relative quality of the partition.</p>
<h3 id=\"wssplot\">Within Sum of Square Plot (WSSPlot)</h3>
<p>This is the plot of the within group sum of squares against the number of clusters. The within group sum of squares decreases as the number of clusters increases. In this plot, the best partition is when adding one or more clusters doesn't decrease the WSS value. It's called the elbow method.</p>
<h4>Example:</h4>
<p><img src=\"GS/WSS.png\"></p>
<p>On this graph, we start by looking at the value for the lowest number of classes: 2. If I add a third class, we see that the WSS value will decrease (from 2.4 to 1.8). If I add another class, this value decreases again (from 1.8 to 1.4). We are looking for the moment when it is not interesting to add a class. In this case, we can keep a partition with 8 or 9 classes.</p>
<h3 id=\"avesilplot\">Average silhouette Plot</h3>
<p>This graph shows the average silhouette width of each partition (ROUSSEEUW 1987). The silhouette width is a bounded index between -1 and 1, calculated for each observation. The closer the value is to 1, the better the observation is classified. We look for the average value for a partition that is closest to 1.</p>
<h4>Example:</h4>
<p><img src=\"GS/AveSil.png\"></p>
<p>On this graph we look for the maximum value. The best evaluation corresponds to the division into 9 classes. Looking at the second best partition, we identify the one with 10 classes.</p>

<h3 id=\"selectalpha\">Select alpha</h3>
<p>This plot shows the correlation criterion for each alpha.</p>
<p><img src=\"GS/selectalpha.png\"></p>

<h3 id=\"dendrogram\">Dendrogram</h3>
<p>The dendrogram is the result of the clustering process. You can see the clusters and sub-clusters.</p>
<p><img src=\"GS/dendrogram.png\"></p>

<h3 id=\"seriograph\">Seriograph</h3>
<p><p>We have chosen the serigraph (DESACHY 2004). This tool makes it possible to highlight the evolution of ceramics over time as well as to understand the commercial relations thanks to the imported ceramics. The percentages of each category of ceramics per set are displayed. The percentages are calculated independently for each set (row). The display of the percentages allows comparison of the different sets but does not provide information on the differences in numbers. To fill this gap, the proportion of the numbers in each class is displayed on the seriograph (weight column).</p>
<p>We have generalized this representation for other contingency data with hclustcompro object.</p>
<p><img src=\"GS/seriograph.png\"></p>
<p>In order to facilitate the exploitation of the data tables, we propose here a computerised graphic processing tool (EPPM serigraph - for Ecart Positif aux Pourcentages Moyens - positive deviation from the average percentage), which does not require specialised statistical skills and is adapted to the case of stratified sites, where the study of the evolution of artefacts can be based on the relative chronology provided by the excavation.</p>
<p>The treatment consists firstly of transforming this table of counts into a table of percentages, the total number in each set (each row) being reduced to 100; these are the proportions, or frequencies, of the types in the sets are thus compared.</p>
<p>The display of positive deviations from the mean percentages (EPPM) shows in black on a lighter background the percentage shares that are higher than the mean percentage of the variable, so as to highlight the most significant part of the values in the table.This display is simply adapted to the seriograph: when a percentage is greater than the average percentage of the type, the excess share (called here EPPM: positive deviation from the average percentage) is shown in black, centred around the axis of the type, on the grey background of the percentage bar.</p>
<p>The table is then transformed into a graphic matrix where these percentages are expressed, for each type, by horizontal bars centred on the column axis. When the rows are ordered chronologically, the silhouette formed by the superposition of these frequency bars bars makes it possible to visualise the evolution over time of the type concerned.</p>
<p>The display of the percentages allows comparison of the different sets but does not provide information on the differences in numbers. To fill this gap, the proportion of the numbers in each class is displayed on the seriograph (weight column).</p>
<p>The processing technique applies to sets whose chronological order is not known; the lines of the graph are to be reorganised so as to obtain boat-shaped silhouettes following the hypothesis of a chronological evolution corresponding to the seriation model.</p>


<h3 id=\"timerangeclust\">Timerange plot</h3>
<p>Show the different time ranges of observations for each cluster (you will need to import time range data).</p>
<p><img src=\"GS/timerangeclust.png\"></p>


<h3 id=\"CAp\">Correspondences analysis</h3>
<p>You can see the different planes of correspondence analysis used to construct the first dissimilarity matrix D1.</p>
<p><img src=\"GS/CA.png\"></p>


<hr>
<h2 id=\"references\">References</h2>
<hr>
<p><br>
  Rousseeuw, P.J. (1987). “Silhouettes: a graphical aid to the interpretation and validation of cluster analysis”.
  In : J. Comput. Appl. Math. 20, 53–65. <a href=\"https://doi.org/10.1016/0377-0427(87)90125-7\" target=\"_blank\">⟨doi:10.1016/0377-0427(87)90125-7⟩</a>.
  <br><br>Ward, J. H., Jr. (1963). \"Hierarchical Grouping to Optimize an Objective Function\".
  Journal of the American Statistical Association, 58, 236–244.
  <br><br>Sokal R. and Michener C. (1958). \"A statistical method for evaluating systematic
  relationships\". University of Kansas Science Bulletin. 38: 1409–1438.
  <br><br>Desachy B. (2004). \"Le sériographe EPPM : un outil informatisé de sériation graphique pour tableaux de comptages. In: Revue archéologique de Picardie, n°3-4, 2004. Céramiques domestiques et terres cuites architecturales. Actes des journées d'étude d'Amiens (2001-2002-2003) pp. 39-56 <a href=\"https://doi.org/10.3406/pica.2004.2396\" target=\"_blank\">⟨doi:10.3406/pica.2004.2396⟩</a>.</p>
<br><br><br><br>



")
    )
  )
)
