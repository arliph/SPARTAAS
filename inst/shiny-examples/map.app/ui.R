library(shinythemes)
library(shinyjs)
library(leaflet)
library(markdown)
library(sp)           #for grid ...
library(cluster)      #for silhouette
library(ggplot2)
library(ape)          #for Moran.i
library(ks)           #for kde Hscv
library("shinyWidgets")
library(mapview)      #mapshot
library(htmltools)

load(file = "./data/datarcheo.RData")
load(file = "./data/datacancer.RData")

#----------------------------------------------------------------------------------------------#
###########################                  - UI -                  ###########################
#----------------------------------------------------------------------------------------------#
ui <- fluidPage(useShinyjs(),style="padding-top: 150px;",theme = shinytheme("flatly"),
                tags$head(
                  tags$title("mapclust")
                ),
                tags$head(tags$style(
                  type = "text/css",
                  ".leaflet-left {
                    bottom:30px;
                  }
                  .btn-default {
                    background-color:#5da3a8;
                    border-color:#5da3a8;
                  }
                  .dropdown-menu li a{
                    background-color:#fff;
                    border-color:#dce4ec;
                  }"
                )),
                #header
                absolutePanel(class = "panel panel-default",
                              style="z-index: 2000;padding: 8px; background: #ecf0f1; opacity: 1;border-bottom: 1px solid #2c3e50;",
                              top = 0, left = 0, right = 0,
                              fixed = TRUE,
                              titlePanel(h2("Divisive Hierarchical Clustering using Spatials Patches")),
                              div(span(strong("SPARTAAS | MapClust ")),span("v1.0 ",style ="font-size:14px;"))
                ),
                #footer
                absolutePanel(style="z-index: 2000;padding: 0px; border-bottom: 0px solid #CCC; background: #fff;opacity: 1;",
                              bottom = 0, left = 0,
                              fixed = TRUE,
                              div(span(strong("SPARTAAS")),span("[Bellanger,Coulon,Husi]",style ="font-size:14px;"))
                ),
                #go top panel
                absolutePanel(class = "panel-hover",
                              style="opacity:0.8;z-index: 1000;padding: 8px; border-bottom: 1px solid #CCC; background: #FFFFEE;width:30px;",
                              bottom = 20, right = 0,
                              fixed = TRUE,
                              HTML("<a href=\"#top\"><span class=\"glyphicon glyphicon-chevron-up\"></span></a>")
                ),
                #-----------------------------------#
                #    NavBar
                #-----------------------------------#
                navbarPage("MapClust",id="navbar",
                           #-----------------------------------#
                           #    Home
                           #-----------------------------------#
                           tabPanel(value = "1",span(icon("home", lib = "glyphicon"),strong("Home")),style="max-width:1200px;",
                                    br(),
                                    sidebarLayout(
                                      #-----------------------------------#
                                      sidebarPanel(style="border-width:0px;border-color:#222222;",
                                                   h1("R Package:"),
                                                   br(),
                                                   HTML("<p>This method is part of the <a href='https://spartaas.frama.io/r-package/index.html'>SPARTAAS</a> package.</p>
                                                        <p>If you are interested you can install our R package avaible on the <a href='https://cran.r-project.org/package=SPARTAAS'>CRAN</a> and on <a href='https://github.com/arliph/SPARTAAS'>GitHub</a>.</p>

                                                        ")

                                      ),
                                      #-----------------------------------#
                                      mainPanel(
                                         HTML("
                                              <h3>Introduction</h3>
                                      <p>
                                              Hierarchical classification can also be considered by introducing constraints that allow, for example, geographical proximity to be taken into account.
                                              </p>
                                              <h4>Spatialpatches</h4>
                                              <p>
                                              The spatial distribution of data can be heterogeneous and present local aggregations or spatial patches. P. Petitgas proposed an algorithm to identify them in the context of fish population density (WOILLEZ et al., 2009). The parameters are the geographical coordinates (X,Y) and a variable of interest (var).
                                              </p>
                                              <p>
                                              The position of a patch is then determined by its center of gravity. The algorithm starts with the highest value of var and then considers each observation in decreasing order of value of var. The highest value initiates the first patch. Then, the observation considered is assigned to the nearest patch, provided that its distance from the center of gravity of the patch is smaller than the threshold distance dlim. Otherwise, the observation forms a new patch. The results on spatial patches are of course influenced by the choice of the dlim threshold and the location of the highest values of var.
                                              </p>
                                              <h4>MapClust</h4>
                                              <p>
                                              In order to better understand the spatial distribution of individuals, spatial patch construction is applied hierarchically top-down by varying the maximum acceptable distance (dlim) between observation points and the centre of gravity of a patch. This approach leads to a Hierarchical Top-down Classification (EVERITT et al., 2001) based at each step on the Spatialpatches algorithm: at each node (parent group) of the hierarchy we divide into two patches (child nodes). Spatialpatches has been developed to describe the spatial distribution patterns of a fish population based on density data, so the variable of interest is assumed to be positive. var can therefore be related to a positive count, frequency or real variable.
                                              </p>
                                              <p>
                                              In the case where the variable of interest var is real, means that it can take negative values, it was necessary to adapt our MapClust classification algorithm. We no longer work with the values of var but with the values of the probability density f associated with var. In the same way we were able to adapt the method to the case where var is multidimensional by working on a Kernel density estimate.
                                              </p>
                                              ")
                                      )
                                    )
                           ),
                           #-----------------------------------#
                           #    MapClust
                           #-----------------------------------#
                           tabPanel(value = "2",span(icon("map"),strong("MapClust")),style="max-width:1200px;",
                                    h1("MapClust"),
                                    sidebarLayout(
                                      sidebarPanel(width = 5,style="border-width:0px;border-color:#222222;",
                                                   wellPanel(width = 5, style="border-width:0px;",
                                                             selectInput("select", label = h4("Dataset:"),
                                                                         choices = list("datarcheo" = 1, "datacancer" = 2, "your data" = 3),
                                                                         selected = 1),
                                                             checkboxInput("Plabel","Print the list of names on the dendrogram.", value = FALSE),
                                                             checkboxInput("extenddend","Merge the leaves located at the same level on the dendrogram (activate before running the method: play button).", value = FALSE),
                                                             hr(style="border-color: #222222;"),
                                                             fluidRow(
                                                               column(width = 3,
                                                                      h4("Run:"),
                                                                      actionButton("Run",icon("play", lib = "glyphicon"))
                                                               ),
                                                               column(width = 6,
                                                                      h4("Import data:"),
                                                                      actionButton("UL",icon("import", lib = "glyphicon"))
                                                               ),
                                                               column(width = 1,
                                                                      h4("Export:"),
                                                                      dropdownButton(
                                                                        downloadButton("all.pdf","Download plots as pdf"),
                                                                        hr(style="border-color: #222222;"),
                                                                        downloadButton("map.png", "Download map as png"),
                                                                        icon = icon("floppy-disk", lib = "glyphicon"), width = "220px",
                                                                        size = "sm",
                                                                        tooltip = tooltipOptions(title = "Download")
                                                                      )
                                                               )
                                                             )
                                                   ),hr(style="border-color: #222222;border-width:2px;"),
                                                   wellPanel(width = 5, style="border-width:0px;background:#dddddd;",
                                                             h4("Dendrogram"),
                                                             plotOutput("distPlot",click = "plot_click")
                                                   ),hr(style="border-color: #222222;border-width:2px;"),
                                                   wellPanel(width = 5, style="border-width:0px;background:#dddddd;",
                                                             h4("Partition evaluation"),
                                                             plotOutput("avgsilPlot"),
                                                             hr(style="border-color: #777777;"),
                                                             plotOutput("WSSPlot")
                                                   )
                                      ),
                                      mainPanel(width = 7,
                                                h4("Console"),
                                                hr(style="border-color: #222222;"),
                                                verbatimTextOutput("click_info"),
                                                h4("Map"),
                                                hr(style="border-color: #222222;"),
                                                leafletOutput("mymap",height = 500),br(),
                                                h4("Silhouette"),
                                                hr(style="border-color: #222222;"),
                                                plotOutput("silPlot",height = "auto"),
                                                h4("Summary"),
                                                hr(style="border-color: #222222;"),
                                                verbatimTextOutput("click_info2")
                                      )
                                    )
                           ),
                           #-----------------------------------#
                           #    Contact
                           #-----------------------------------#
                           tabPanel(value = "5",span(icon("envelope", lib = "glyphicon"),strong("Contact")),style="max-width:1200px;",
                                    HTML("<style>
                                       .carte{
                                       border-left-width: 4px;
                                       border-left-style:solid;
                                       border-color: #2ac0a2;
                                       padding-left: 10px;
                                       }
                                       </style>
                                       <h3>Authors:</h3>
                                       <hr>

                                       <div class=\"carte\"><h4>L. Bellanger</h4><h5>mail: &lt;lise.bellanger@univ-nantes.fr&gt;</h5></div>
                                       <div class=\"carte\"><h4>P. Husi</h4><h5>mail: &lt;philippe.husi@univ-tours.fr&gt;</h5></div><div class=\"carte\"><h4>A. Coulon</h4></div>

                                       <h3>Maintainer:</h3>
                                       <hr>
                                       <div class=\"carte\"><h4>A. Coulon</h4><h5>mail: &lt;arthur-coulon@outlook.fr&gt;</h5></div>
                                       ")
                                                ),
                           #-----------------------------------#
                           #    Wiki
                           #-----------------------------------#
                           tabPanel(value = "6",span(icon("question-sign", lib = "glyphicon"),strong("Wiki")),style="max-width:1200px;",
                                    HTML(
                                      "<style>body {text-align: justify}</style><h1>Get started with the MapClust application</h1>
                                      <h2>Table of content</h2>
                                      <div id=\"TOC\">
                                      <ul>
                                      <li><a href=\"#introduction\">Introduction</a><ul>
                                      <li><a href=\"#spatialpatches\">Spatialpatches</a></li>
                                      <li><a href=\"#mapclust\">MapClust</a></li>
                                      </ul></li>
                                      <li><a href=\"#data\">The App</a><ul>
                                      <li><a href=\"#first\">Run (first time)</a></li>
                                      <li><a href=\"#plabel\">Print label option</a></li>
                                      <li><a href=\"#import\">Import your data</a><ul>
                                      <li><a href=\"#csv\">CSV Format and write.table</a></li>
                                      <li><a href=\"#label\">Label</a></li>
                                      <li><a href=\"#header\">Header</a></li>
                                      <li><a href=\"#separator\">Separator</a></li>
                                      <li><a href=\"#quote\">Quote</a></li>
                                      <li><a href=\"#dec\">Decimal</a></li>
                                      <li><a href=\"#uni\">Univariate data</a></li>
                                      <li><a href=\"#multi\">Multivariate data</a></li>
                                      </ul></li>
                                      </ul></li>
                                      <li><a href=\"#eval\">Evaluation Plot</a><ul>
                                      <li><a href=\"#wssplot\">Within Sum of Square Plot</a></li>
                                      <li><a href=\"#avesilplot\">Average Silhouette Plot</a></li>
                                      </ul></li>
                                      </li>
                                      <li><a href=\"#output\">Output: select a partition</a><ul>
                                      <li><a href=\"#dendrogram\">Dendrogram</a></li>
                                      <li><a href=\"#map\">Map</a></li>
                                      <li><a href=\"#silhouette\">Silhouette</a></li>
                                      <li><a href=\"#summary\">Additionnal informations</a></li>
                                      </ul></li>
                                      <li><a href=\"#references\">References</a></li>
                                      </ul>
                                      </div>

<hr>
<h2 id=\"introduction\">Introduction</h2>
<hr>
<p>Hierarchical classification can also be considered by introducing constraints that allow, for example, geographical proximity to be taken into account.</p>
<h3 id=\"spatialpatches\">Spatialpatches</h3>
<p>The spatial distribution of data can be heterogeneous and present local aggregations or spatial patches. P. Petitgas proposed an algorithm to identify them in the context of fish population density (WOILLEZ et al., 2009). The parameters are the geographical coordinates (X,Y) and a variable of interest (var).</p>
<p>The position of a patch is then determined by its center of gravity. The algorithm starts with the highest value of var and then considers each observation in decreasing order of value of var. The highest value initiates the first patch. Then, the observation considered is assigned to the nearest patch, provided that its distance from the center of gravity of the patch is smaller than the threshold distance dlim. Otherwise, the observation forms a new patch. The results on spatial patches are of course influenced by the choice of the dlim threshold and the location of the highest values of var.</p>
<h3 id=\"mapclust\">MapClust</h3>
<p>In order to better understand the spatial distribution of individuals, spatial patch construction is applied hierarchically top-down by varying the maximum acceptable distance (dlim) between observation points and the centre of gravity of a patch. This approach leads to a Hierarchical Top-down Classification (EVERITT et al., 2001) based at each step on the Spatialpatches algorithm: at each node (parent group) of the hierarchy we divide into two patches (child nodes). Spatialpatches has been developed to describe the spatial distribution patterns of a fish population based on density data, so the variable of interest is assumed to be positive. var can therefore be related to a positive count, frequency or real variable.</p>
<p>In the case where the variable of interest var is real, means that it can take negative values, it was necessary to adapt our MapClust classification algorithm. We no longer work with the values of var but with the values of the probability density f associated with var. In the same way we were able to adapt the method to the case where var is multidimensional by working on a Kernel density estimate.</p>


<hr>
<h2 id=\"data\">The App</h2>
<hr>
<h3 id=\"first\">Run (first time)</h3>
<p>The function therefore requires as data geographical coordinates, longitude and latitude (<code>coord</code>) and the variable(s) of interest (<code>var</code>) in the same table.</p>
<p>For your first time you can use our dataset. The first one it's call datarcheo. Just select it on the MapClust tab(it will be select by default) and run. You also can try datacancer the other dataset.</p>
<h3 id=\"plabel\">Print label option</h3>
<p>Placed just under the choice of the dataset, it is used to define if the labels should be displayed on the dendrogram.</p>
<p><img src=\"GS/Plabel.png\">Select dataset and print label option</p>
<h3 id=\"import\">Import your data</h3>
<p>You can import your data. Since they may contain only one or more variables of interest, you must adjust the parameters.</p>
<p>Let's take a look of the interface. There are three main parts. The import space with only one import button. The second part is more complex with more options. This is the parameter part where you must configure the import tool according to your data. Finally, there is the preview area where you can see how your data should look and how your data is. There is still something left, the concordance indicator that allows you to quickly check if your data matches what it should look like.</p>
<p><img src=\"GS/import.png\"></p>
<h4 id=\"csv\">CSV Format and write.table</h4>
<p><img src=\"GS/csv.png\"><br>A csv file. It is a data.frame with 3 colunms with headers and separate by semicolon \";\".</p>
<p>The input format for importing data is the \".csv\" format but also supports the\".txt\" format as a csv file.</p>
<p>In R you can export your data frame into a csv file using write.csv2 or write.table. In a csv you can choose a character to separate the columns. In the same way, you can define the character to indicate the decimal point.</p>
<code>
write.table(data,file=\"path/to/name_file.csv\",sep=\";\",dec=\".\",row.names=FALSE,quote=FALSE)
</code>
<p>In Excel you can save as CSV format in order to import your data frame.</p>
<p>The import interface allows you to setup this values with the \"header\", \"decimal\", \"separator\" and \"quote\" option.</p>
<p>In order to compute exactly as you want we need to know some things. The first one is the presence or not of label in your data. The label have to be in the last colunm. We also need to know how many of variable of interest you have in addition of the two coordinates variables.</p>
<p>If the setup don't match with your real data frame you will not be able to import anything. You have to change the data or the parameters. When it was good a green check icon will appear if not you will see a red cross.</p>
<p><img src=\"GS/importsetting.png\"></p>
<h4 id=\"label\">Label</h4>
<p>Yes or not option. Do you have labels in your data ? If you have labels you have to put them on the last colunms.</p>
<h4 id=\"header\">Header</h4>
<p>Yes or not option. Do you have headers on your colunms ?</p>
<h4 id=\"separator\">Separator</h4>
<p>Choose the character use to separate the colunms.</p>
<h4 id=\"quote\">Quote</h4>
<p>Choose the quote use to strings.</p>
<h4 id=\"dec\">Decimal</h4>
<p>Choose the character use to indicate decimal.</p>
<h4 id=\"uni\">Univariate data</h4>
<p>Put the slider to \"Univariate\" then import your data. If you have labels let the option check, if not uncheck label.</p>
<p>Configure all settings for the csv format. Which symbol is used to separate the columns, which decimal symbol, which quotation marks to use for the character string and whether your columns contain a header.</p>
<p><img src=\"GS/importUni.png\"></p>
<h4 id=\"multi\">Multivariate data</h4>
<p>The only différence is that you have to put the slider on \"Multivariate\". When you do this, you hace acces to another slider. You have to indicate the number of variables of interest. Warning: You must not count the two coordinates variables in.</p>
<p><img src=\"GS/importMulti.png\"></p>


<hr>
<h2 id=\"eval\">Evaluation Plot</h2>
<hr>
<p>It is essential to be able to evaluate the different partitions of the hierarchy in order to identify the one(s) that is (are) most relevant. The number k of classes of the partition to be retained is based in our case on several indicators calculated for different values of K : the total within-class sum of square (<code>WSS</code>) and the global average of the silhouette widths.</p>
<p>During the execution of the MapClust method you have to make a choice. You have to cut the dendrogram. This operation select the partition. In order to compare all the possibilities you can see the evaluation plot (WSSPlot and AveSilPlot). This two plot evaluate the relative good quality of the partition.</p>
<h3 id=\"wssplot\">Within Sum of Square Plot (WSSPlot)</h3>
<p>This is the plot of within-groups sum of squares against number of clusters. The Within Sum of Square decrease when the number of cluster increase. In this plot the best partition is when add one or more clusters don’t decrease the WSS value. It’s call the Elbow method.</p>
<h4>Example:</h4>
<p><img style=\"width:600px;\" src=\"GS/WSS2.png\"></p>
<p>On this graph we start by looking at the value for the lowest number of groups: 2. if I add a third group we see that the WSS value will decrease (from 0.09 to 0.03). If I add another group I will decrease this value again (from 0.03 to 0.01). After that, adding a group no longer or only slightly affects the value. Adding a group is therefore not interesting, we keep a partition with 4 groups.</p>
<h3 id=\"avesilplot\">Average silhouette Plot</h3>
<p>This graph shows the average silhouette width of each partition (ROUSSEEUW 1987). The silhouette width is a limited index between -1 and 1, which is calculated for each observation. The closer the value is to 1, the better the observation is classified. We look for the average value for a partition closest to 1.</p>
<h4>Example:</h4>
<p><img style=\"width:600px;\" src=\"GS/AveSil2.png\"></p>
<p>On this graph we look for the maximum value. The best evaluation corresponds to the division into 7 groups. Looking at the second best partition we identify the one with 4 groups. Although the one with 7 groups is higher for this silhouette index we will select the partition in 4. we make this choice because of the WSSPlot which advised us the partition in 4.</p>


<hr>
<h2 id=\"output\">Output</h2>
<hr>
<h4>The selection of a partition is done by clicking on the dendrogram at the desired height.</h4>
<p>You can change it at any moment.</p>
<div id=\"dendrogram\" style=\"display:inline-block\">
  <img src=\"GS/dendro.png\">
  <img src=\"GS/dedroCut.png\">
</div><br>
<h3 id=\"map\">Map</h3>
<p><img src=\"GS/map.png\"></p><br>
<h3 id=\"silhouette\">Silhouette</h3>
<p><img src=\"GS/silhouette.png\"></p>
<p>You can see on this plot the silhouette index for each observation of the selected partition. The observation are sorted in decreasing order and by cluster.</p><br>
<h3 id=\"summary\">Additionnal informations</h3>
<p>above the map you can see a little description of the cluster. You can find the average value of each variable.</p>
<p><img src=\"GS/consol.png\"></p>
<p>Below the silhouette plot you have a summary of the diffrents values for each partition. You can find the dlim, the number of cluster, the WSS value the average sil_widht and if you are in univariate case the Moran index (MORAN, 1950) for each partition. The Moran index is a measure of spatial autocorrelation. It is calculated for each cluster of the partition and then the average is calculated to characterize the partition.</p>
<p><img src=\"GS/summary.png\"></p>


<hr>
<h2 id=\"references\">References</h2>
<hr>
<p>EVERITT, B.S., S. LANDAU et M. LEESE(2001). “Cluster Analysis”. In : London UK: Arnold.</p>
<p>MORAN, P. A. P. (1950). \"Notes on Continuous Stochastic Phenomena\". Biometrika. 37 (1): 17–23. <a href=\"https://doi.org/10.2307/2332142\" class=\"uri\">doi.org/10.2307/2332142</a>. JSTOR 2332142.</p>
<p>ROUSSEEUW, P.J. (1987). “Silhouettes: a graphical aid to the interpretation and validation of cluster analysis”. In : J. Comput. Appl. Math. 20, 53–65. <a href=\"DOI:https://doi.org/10.1016/0377-0427(87)90125-7\" class=\"uri\">DOI:https://doi.org/10.1016/0377-0427(87)90125-7</a>.</p>
<p>WOILLEZ, M., J.RIVOIRARD et P.PETITGAS(2009). “Notes on survey-based spatial indicators for monitoring fish populations”. In : Aquat. Living Resour.22, p. 155–164.</p>

                                      "
                                    )
                                    )
                           )
                           )
