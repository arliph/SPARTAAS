library(SPARTAAS)
suppressMessages(library(shinyjs))
library(shinyjqui)
suppressMessages(library(shinydashboard))
library(colorspace)
suppressMessages(library(dplyr))
library(tidyr)
library(ggplot2)
library(shinyWidgets)
library(shinythemes)

ui <- fluidPage(useShinyjs(),style="padding-top: 150px;",theme = shinytheme("flatly"),
                tags$style(HTML('#sort,#data,#subdivise,#tri,#sub1,#sub2,
         #import2,#import3,#import3,#Finishimport,#back1,#back2{color:white;background-color:#337ab7}
         #back{color:white;background-color:#f0ad4e}
         #tabsetperso{background-color:#e5e5e5;}
         #help{text-align;right;}
         #tabsetperso{margin-bottom:20px;}
         .dropdown-menu > li > a {
                    background-color:#fff;
                    border-color:#dce4ec;
         }
                  .btn-default {
                  border-color:#337ab7;
                  }
    ')),

                tags$head(
                  tags$link(rel = "icon", type = "image/gif", href = "https://spartaas.gitpages.huma-num.fr/r-package/img/lambda.png"),
                  tags$title("seriograph")
                ),
                #header panel
                absolutePanel(class = "panel panel-default",
                              style="z-index: 2000;padding: 8px; background: #ecf0f1; opacity: 1;border-bottom: 1px solid #2c3e50;",
                              top = 0, left = 0, right = 0,
                              fixed = TRUE,
                              h2("seriograph: a graphic seriation tool for count tables"),
                              div(span(strong("SPARTAAS | seriograph")))
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
                navbarPage("seriograph",id="navbar",
                           tabPanel(value = "1",span(icon("home", lib = "glyphicon"),strong("Home")),style="max-width:1200px;",
                                    br(),
                                    sidebarLayout(
                                      #-----------------------------------#
                                        sidebarPanel(style="",
                                                     h1("R Package:"),
                                                     br(),
                                                     HTML("<p>This method is part of the <a href='https://spartaas.gitpages.huma-num.fr/r-package/index.html' target='_blank'>SPARTAAS</a> package.</p>
                                                        <p>If you are interested you can install our R package available on the <a href='https://cran.r-project.org/package=SPARTAAS' target='_blank'>CRAN</a> and on <a href='https://github.com/arliph/SPARTAAS' target='_blank'>GitHub</a>.</p>
                                                        <p>There is also a macro version in LibreOffice Calc. You can find it <a href='https://abp.hypotheses.org/le-programme-bassin-parisien/les-projets/les-projets-associes-au-programme/outils-danalyse-graphique-des-donnees' target='_blank'>here</a>.</p>
                                                        <hr style=\"border-color: #222222;\">
                                                        <p><b>Desachy B.</b> (2004) Le sériographe EPPM : un outil informatisé de sériation graphique pour tableaux de comptages. In: Revue archéologique de Picardie, n°3-4, 2004. Céramiques domestiques et terres cuites architecturales. Actes des journées d'étude d'Amiens (2001-2002-2003) pp. 39-56 <a href=\"https://doi.org/10.3406/pica.2004.2396\" target=\"_blank\">⟨doi:10.3406/pica.2004.2396⟩</a>.</p>")

                                        ),
                                      #-----------------------------------#
                                      mainPanel(
                                        HTML("
                                              <h3>Introduction</h3>
                                              <p>In order to facilitate the exploitation of the data tables, we propose here a computerised graphic processing tool (EPPM serigraph - for Ecart Positif aux Pourcentages Moyens - positive deviation from the average percentage), which does not require specialised statistical skills and is adapted to the case of stratified sites, where the study of the evolution of artefacts can be based on the relative chronology provided by the excavation.</p>
                                              <p>The treatment consists firstly of transforming this table of counts into a table of percentages, the total number in each set (each row) being reduced to 100; these are the proportions, or frequencies, of the types in the sets are thus compared.</p>
                                              <p>The display of positive deviations from the mean percentages (EPPM) shows in black on a lighter background the percentage shares that are higher than the mean percentage of the variable, so as to highlight the most significant part of the values in the table.This display is simply adapted to the seriograph: when a percentage is greater than the average percentage of the type, the excess share (called here EPPM: positive deviation from the average percentage) is shown in black, centred around the axis of the type, on the grey background of the percentage bar.</p>
                                              <p>The table is then transformed into a graphic matrix where these percentages are expressed, for each type, by horizontal bars centred on the column axis. When the rows are ordered chronologically, the silhouette formed by the superposition of these frequency bars bars makes it possible to visualise the evolution over time of the type concerned.</p>
                                              <p>The display of the percentages allows comparison of the different sets but does not provide information on the differences in numbers. To fill this gap, the proportion of the numbers in each class is displayed on the seriograph (weight column).</p>
                                              <p>The processing technique applies to sets whose chronological order is not known; the lines of the graph are to be reorganised so as to obtain boat-shaped silhouettes following the hypothesis of a chronological evolution corresponding to the seriation model.</p>

                                              ")
                                      )
                                    )
                           ),
                           tabPanel(value = "2",span(icon("chart-area", lib = "font-awesome"),strong("seriograph")),style="max-width:1200px;",
                                    sidebarPanel(style = 'max-width:300px;',
                                                 h4("Import your data:"),
                                                 column(10,
                                                        actionButton("data",
                                                                                      icon("import", lib = "glyphicon"))
                                                 ),
                                                 column(1,
                                                        # Help
                                                        dropdownButton(
                                                          h3("data format"),
                                                          hr(style="border-color: #222222;"),
                                                          HTML("<p>The data is in the form of a count table with the contexts in rows and the categories (GT) in columns.</p>
                                                       <p>Example:</p>
                                                       <img style='width:280px;' src='GS/exemple.png' alt='data'>"),
                                                          circle = TRUE, status = "danger", icon = icon("question"), width = "300px",
                                                          size = "sm",right=FALSE,
                                                          tooltip = tooltipOptions(title = "Help", placement = "top")
                                                        )
                                                        # Help
                                                 ),
                                                 div(style="height:200px;"),
                                                 prettyToggle(
                                                   inputId = "permute",
                                                   label_on = "Seriation",
                                                   label_off = "Seriation",
                                                   icon_on = icon("check"),
                                                   icon_off = icon("times"),
                                                   value = TRUE
                                                 ),
                                                 selectInput("show", label = "Show", choices = list("both" = 1,
                                                                                                    "frequencies" = 2,
                                                                                                    "EPPM" = 3),
                                                             selected = 1),
                                                 actionButton("sort",span(icon("sort",lib="font-awesome"),"Sort rows"))
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
                                               numericInput("width", "width", 900, min = 0, max = NA,
                                                            step = 1, width = NULL),
                                               numericInput("height", "height", 500, min = 0, max = NA,
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
                                               h3("The seriograph (B. DESACHY)"),
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
                                               for other contingency data with hclustcompro object.
                                               The seriograph can be used in a strictly deductive way, by blocking the order of the rows (archaeological contexts) when their succession is known, and by reordering only the columns, in order to examine the chronological behaviour of the variables, and thus to see in particular what does not obey a serial evolution: cyclical phenomena, or the problems of intrusion and residuality."),
                                               circle = TRUE, status = "danger", icon = icon("question"),
                                               width = "500px",
                                               size = "sm",right=TRUE,
                                               tooltip = tooltipOptions(title = "Help", placement = "top")
                                             )
                                             # Help
                                      ),br(),br(),

                                      jqui_resizable(plotOutput("EPPM")),
                                      div(style="height:20px;"),
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

                                      <li><a href=\"#output\">The outputs</a>
                                        <ul>
                                          <li><a href=\"#seriograph\">Seriograph</a></li>
                                        </ul>
                                      </li>

                                      <li><a href=\"#data1\">The inputs</a>
                                        <ul>
                                          <li><a href=\"#seriation\">Seriation</a></li>
                                          <li><a href=\"#weight\">Weight color indicator</a></li>
                                          <li><a href=\"#show1\">Vizualisation</a></li>
                                          <li><a href=\"#sort1\">Sort periods</a></li>
                                          <li><a href=\"#import\">Import your data</a></li>
                                          <li><a href=\"#csv\">CSV Format and write.table</a></li>
                                        </ul>
                                      </li>


                                      <li><a href=\"#references\">References</a></li>


                                      </ul>
                                      </div>

<hr>
<h2 id=\"output\">The outputs</h2>
<hr>
<h3 id=\"seriograph\">Seriograph</h3>
<p>We have chosen the serigraph (Desachy 2004). This tool makes it possible to highlight the evolution of ceramics over time as well as to understand the commercial relations thanks to the imported ceramics. The percentages of each category of ceramics per set are displayed. The percentages are calculated independently for each set (row). The display of the percentages allows comparison of the different sets but does not provide information on the differences in numbers. To fill this gap, the proportion of the numbers in each class is displayed on the seriograph (weight column).</p>
<p>The seriograph can be used in a strictly deductive way, by blocking the order of the rows (archaeological contexts) when their succession is known, and by reordering only the columns, in order to examine the chronological behaviour of the variables, and thus to see in particular what does not obey a serial evolution: cyclical phenomena, or the problems of intrusion and residuality.</p>

<p><img src=\"GS/seriograph.png\"></p>

<p>In order to facilitate the exploitation of the data tables, we propose here a computerised graphic processing tool (EPPM serigraph - for Ecart Positif aux Pourcentages Moyens - positive deviation from the average percentage), which does not require specialised statistical skills and is adapted to the case of stratified sites, where the study of the evolution of artefacts can be based on the relative chronology provided by the excavation.</p>
<p>The treatment consists firstly of transforming this table of counts into a table of percentages, the total number in each set (each row) being reduced to 100; these are the proportions, or frequencies, of the types in the sets are thus compared.</p>
<p>The display of positive deviations from the mean percentages (EPPM) shows in black on a lighter background the percentage shares that are higher than the mean percentage of the variable, so as to highlight the most significant part of the values in the table.This display is simply adapted to the seriograph: when a percentage is greater than the average percentage of the type, the excess share (called here EPPM: positive deviation from the average percentage) is shown in black, centred around the axis of the type, on the grey background of the percentage bar.</p>
<p>The table is then transformed into a graphic matrix where these percentages are expressed, for each type, by horizontal bars centred on the column axis. When the rows are ordered chronologically, the silhouette formed by the superposition of these frequency bars bars makes it possible to visualise the evolution over time of the type concerned.</p>
<p>The display of the percentages allows comparison of the different sets but does not provide information on the differences in numbers. To fill this gap, the proportion of the numbers in each class is displayed on the seriograph (weight column).</p>
<p>The processing technique applies to sets whose chronological order is not known; the lines of the graph are to be reorganised so as to obtain boat-shaped silhouettes following the hypothesis of a chronological evolution corresponding to the seriation model.</p>



<hr>
<h2 id=\"data1\">The inputs</h2>
<hr>
<h3 id=\"seriation\">Seriation</h3>
<p>This input allows you to enable or disable the seriation of the columns. The matrix permutation uses an algorithm called 'reciprocal averages'. Each row is assigned a rank from 1 to n, the number of rows. For each column, a barycentre is calculated by weighting according to the row rank. Finally, the columns are reorganised by sorting them according to their barycentre.</p>
<p><img src=\"GS/seriation.png\"></p>


<h3 id=\"show1\">Visualization</h3>
<p>This input allows you to select the element to be plotted. There are tree options: plot the positive deviation from the average percentage (EPPM in French), plot the frequency or plot both. The average percentage is calculated for each category (columns) on the total number of accounts (all classes combined). From the average percentage, we obtain for each category and for each row the difference between the percentage of the category in the class and the average percentage. The EPPM corresponds to the notion of independence deviation (between rows and columns, between categories and time classes) in a chi-square test approach. Although this approach is fundamental in statistical analysis, the independence deviations here are purely indicative and are not associated with a p_value that could determine the significance of the deviations.</p>
<p><img src=\"GS/show.png\"></p>

<h3 id=\"sort1\">Sort periods</h3>
<p>The rows are initially in the order of the data table. It is possible to reorder the rows in a temporal way. In the interface you can drag and drop the row (green) to change the order.</p>
<p><img src=\"GS/sort.png\"></p>

<h3 id=\"import\">Import your data</h3>
<p>You can import your data. You will need to upload a csv for the contingency table.</p>
<p>The data is in the form of a count table with the contexts in rows and the categories (GT) in columns.</p>
<img src='GS/exemple.png'>
<p><img src=\"GS/import1.png\"></p>

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
<p>It is a data.frame with colunms separated by semicolons \";\".</p>
<p>The input format for importing data is the .csv format, but also supports the .txt format as a .csv file.</p>
<p>In R, you can export your data frame to a csv file using write.csv2 or write.table. In a csv you can choose a character to separate the columns. In the same way, you can define the character to indicate the decimal point.</p>
<code>
write.table(data,file=\"path/to/name_file.csv\",sep=\";\",dec=\".\",row.names=FALSE,quote=FALSE)
</code>
<p>In Excel you can save in csv format in order to import your data frame.</p>
<p>The import interface allows you to set these values using the 'header', 'decimal', 'separator' and 'quote' options.</p>




<hr>
<h2 id=\"references\">References</h2>
<hr>
<p><br>
  <br><br>Desachy B. (2004). Le sériographe EPPM : un outil informatisé de sériation graphique pour tableaux de comptages. In: Revue archéologique de Picardie, n°3-4, 2004. Céramiques domestiques et terres cuites architecturales. Actes des journées d'étude d'Amiens (2001-2002-2003) pp. 39-56 <a href=\"https://doi.org/10.3406/pica.2004.2396\" target=\"_blank\">⟨doi:10.3406/pica.2004.2396⟩</a>.
</p>
<br><br><br><br>


")
                           )
                )
)
