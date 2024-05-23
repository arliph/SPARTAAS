estimateDensityMulti <- function(data){
  cat("\nEstimation of multidimensionnal density ..")
  #Smooth cross-validation bandwidth selector
  H.scv <- Hscv(data)
  fhat <- kde(data, H.scv, eval.points=data)
  res <- fhat$estimate*100
  cat(". Done")
  return(res)
}

IW <- function(x)
{
  #===============================================================================
  # SOMMES DES ECARTS AUX CARRE
  #===============================================================================
  sum(scale(x, scale = FALSE)^2)
}

spatialpatches2 <- function (x, y, z, w, Lim.D = 100, A.li = 0, modproj = NA, mlong = NA, mlat = NA)
{
  #===============================================================================
  # IDENTIFICATION OF SPATIAL PATCHES
  # AND COUNT OF CENTERS OF GRAVITY OF SPATIAL PATCHES
  #
  # Routine from EU program Fisboat, DG-Fish, STREP n? 502572
  # Authors : P.Petitgas (Ifremer), M.Woillez and J.Rivoirard (Mines-ParisTech)
  # Last update : 01 march 2008
  #
  # Arguments:
  # x	      The x-coordinate of the first population
  #	        Can be a vector, a matrix or an xyZ list (-,|,[]).
  # y	      The y-coordinates of the first population
  # z	      The regionalised variable of the first population.
  #	        If missing, the results of 'cgi' will concern the samples only.
  # w	      Optional. A weight or an area of influence of the first population.
  #	        Set to 1 if missing
  # Lim.D   Select minimum distance from sample to patch centre: to
  #         identify patches (units are those of coordinates)
  # A.li    Visualisation of gravity centres for those patches with
  #         abundance > A.li (in %)
  # modproj	Optional. Indicates the type of projection to perform.
  # mlong   mean longitude in DEGREES of the data set to be transformed
  # mlat    mean latitude in DEGREES of the data set to be transformed
  #	        See 'dg2nm' for precisions
  #
  #===============================================================================

  patch.id <- function(xx,yy,zz,ww,dli,ali,modproj,mlong,mlat)
  {
    ############################################################################
    # inputs : xx , yy , ww & zz ranked by decreasing order of zz
    #          dli is a distance limit from the patch gravity centre to define
    #          the patch border
    #
    # the function starts from the richest zz value and considers each sample
    # in decreasing order of zz. It tests whether the current value is spatially
    # close enough to the gravity centre of previously formed patches. if not a
    # new patch is considered. and so on until the last value. Patches of nul
    # values are returned with centres as NA and code 0 and their areas are
    # summed.
    #
    # outputs : the patch number for each sample,
    #           the gravity center for each patch,
    #           the percent abundance and area for each patch
    ############################################################################

    g <- c(1,rep(0,(length(xx)-1)))
    for (j in 1:(length(xx)-1)) {
      xg <- tapply(ww[1:j]*zz[1:j]*xx[1:j],g[1:j],sum,na.rm=T) /
        tapply(ww[1:j]*zz[1:j],g[1:j],sum,na.rm=T)
      yg <- tapply(ww[1:j]*zz[1:j]*yy[1:j],g[1:j],sum,na.rm=T) /
        tapply(ww[1:j]*zz[1:j],g[1:j],sum,na.rm=T)
      d <- sqrt(((xg-x1[j+1]))^2 + (yg-y1[j+1])^2)
      o <- order(d)
      if (d[o[1]] < dli){
        g[j+1] <- o[1]
      }else{
        g[j+1] <- max(g[1:j], na.rm=T)+1
      }
    }
    #second passage des obs pour affiner centres de gravite
    for(i in 1:length(xx)) {
      xg <- tapply(ww*zz*xx,g,sum,na.rm=TRUE)/tapply(ww*zz,g,sum,na.rm=TRUE)
      yg <- tapply(ww*zz*yy,g,sum,na.rm=TRUE)/tapply(ww*zz,g,sum,na.rm=TRUE)
      d <- sqrt(((xg-x1[i]))^2+(yg-y1[i])^2)
      o <- order(d)
      g[i] <- o[1]
    }

    xg <- tapply(ww*zz*xx,g,sum,na.rm=TRUE)/tapply(ww*zz,g,sum,na.rm=TRUE)
    yg <- tapply(ww*zz*yy,g,sum,na.rm=TRUE)/tapply(ww*zz,g,sum,na.rm=TRUE)
    pb <- tapply(ww*zz,g,mean,na.rm=TRUE)
    inertiaTot <- sum(scale(zz, scale = FALSE)^2)
    InertiaWithin <- tapply(zz,g,IW)
    sel <- (tapply(ww*zz,g,sum,na.rm=TRUE)/sum(ww*zz,na.rm=TRUE))*100>=ali
    res <- max(unique(g)[sel],na.rm=TRUE)
    n <- sort(unique(g))
    selna <- is.na(xg)
    if (sum(selna) > 1) {
      n <- c(n[!selna],0)
      pb <- c(pb[!selna],0)
      pa <- c(pa[!selna],sum(pa[selna],na.rm=T))
      xg <- c(xg[!selna],NA)
      yg <- c(yg[!selna],NA)
      g[!is.element(g,n)] <- 0
    }
    return(list(n=g,mat=cbind(n=n,xg=round(xg,2),yg=round(yg,2),zg=round(pb*100,2)),
                InertiaWithin=InertiaWithin,
                InertiaWithintot=sum(InertiaWithin),
                InertiaBetween=inertiaTot - sum(InertiaWithin),
                inertiaTot=inertiaTot,
                nsp=res))
  }

  # prepare data
  Xsta <- x		#Xsta <- bid$x
  Ysta <- y		#Ysta <- bid$y
  Zsta <- z
  Wsta <- w/sum(w,na.rm=TRUE)

  # exclude NA values
  SEL <- is.na(Xsta)+is.na(Ysta)+is.na(Zsta)+is.na(Wsta)
  Xsta <- Xsta[!SEL]
  Ysta <- Ysta[!SEL]
  Zsta <- Zsta[!SEL]
  Wsta <- Wsta[!SEL]

  # order data by decreasing order of value z
  zz1 <- sort(Zsta,decreasing=TRUE,index.return=TRUE)
  w1 <- Wsta[zz1$ix]
  x1 <- Xsta[zz1$ix]
  y1 <- Ysta[zz1$ix]
  z1 <- zz1$x

  # identify patches around high values
  SP <- patch.id(x1,y1,z1,w1,Lim.D,A.li,modproj=modproj,mlong=mlong,mlat=mlat)

  # graphical visualisation : patch identification
  #shows data points, circles for data values, number of patch &
  #crosses for patch gravity centres with abundance > A.li of total
  o <- sort(zz1$ix,index.return=TRUE)
  SP$n <- SP$n[o$ix]
  SP
}

WSS <- function(hist,data)
{
  #====================================================================================
  # WSS        Return the within cluster sum of square for each partition.
  #
  # hist       The historique of the partitions
  # data       The geographical data and third variable Z
  #====================================================================================

  out <- list()
  n <- length(hist[1,])
  Itot <- calctot(data)
  for (i in 1:n){
    WSS <- c()
    #recup compo of one partition i
    compo <- hist[,i]

    #calculer Within-cluster Sum of Square
    X <- data [,1]
    Y <- data [,2]
    Z <- data [,3]
    XG <- tapply(X * Z,compo,sum)/tapply(Z,compo,sum)
    YG <- tapply(Y * Z,compo,sum)/tapply(Z,compo,sum)
    nbgrp <- length(XG)
    for(j in 1:nbgrp){
      subX <- X[compo==j]
      subY <- Y[compo==j]
      subZ <- Z[compo==j]
      DX <- subX - XG[j]
      DY <- subY - YG[j]
      D <- sqrt(DX^2 + DY^2)
      INERT <- sum(subZ * (D^2))/sum(Z)
      WSS <- c(WSS,INERT)
    }
    out <- c(out,list(WSS))
  }
  out
}

calctot <- function(data)
{
  X <- data[,1]
  Y <- data[,2]
  Z <- data[,3]
  XG <- sum(X * Z)/sum(Z)
  YG <- sum(Y * Z)/sum(Z)
  DX <- X - XG
  DY <- Y - YG
  D <- sqrt(DX^2 + DY^2)
  INERT <- sum(Z * (D^2))/sum(Z)
  INERT
}

dist_m <- function(data)
{
  #===================================================================================
  # DISTANCE MATRIX   Weighted Eucidiean distance
  #
  # data      data.frame or matrix with 3 colunms (X, Y, Z)
  #===================================================================================
  square <-function(x){return(x*x)}
  Dmatrix <- matrix(ncol=length(data[,1]), nrow=length(data[,1]))
  for(i in 1:length(data[,1])){
    for (j in 1:length(data[,1])){
      #Weighted Eucidiean distance
      distance_pond <- sqrt( ( square(data[i,1]-data[j,1]) + square(data[i,2]-data[j,2]) ) *
                               square(data[i,3]-data[j,3])
      )
      Dmatrix[i,j] <- distance_pond
      Dmatrix[j,i] <- distance_pond
    }
  }
  return(Dmatrix)
}

transfo_z <- function(z)
{
  #====================================================================================
  # Transfo Z  Confidence interval on Empirical Cumulative Distribution Function (ECDF)
  #            We only can use z data positive like density, count or frequencies.
  #
  # z         vector of values
  #====================================================================================
  cat("\nEstimation of univariate density ..")
  dx <- density(z)
  xnew <- z # find density value on a new point
  dx.approx <- approx(dx$x,dx$y,xout=xnew)
  pi <- dx.approx$y*100 # densite estimee * 100 en chaque point du jeu de donnees
  cat(". Done")
  return(pi)
}

addLegendCustom <- function(map, colors, labels, sizes, opacity = 0.5, title)
{
  #===================================================================================
  # Custom legend for diffrent sizes of point in final map
  #
  # Arguments:
  # map       The leaflet map object
  # colors    The colors of the points
  # labels    The labels of the points
  # sizes     The sizes signification of the points in the map
  # opacity   The opacity of the background
  #
  # You can change the title and the position of the legend directly in the code. (last line of the function)
  #===================================================================================

  colorAdditions <- paste0(colors, ";border-radius: 50%;margin-top:10%; width:", 2*sizes, "px; height:", 2*sizes, "px")
  labelAdditions <- paste0("<div style='display: inline-block;height: ",
                           2*sizes, "px;margin-top: 4px;line-height: 0px;'>", labels, "</div>")

  return(addLegend(map, colors = colorAdditions, labels = labelAdditions, opacity = opacity,
                   title = title, position = "topright"))
}

"ecartmatrice" <-
  function (historique = "matrix composition of groups for each dlim")
  {
    #===================================================================================
    # ECART MATRICE
    #
    # Author : A. Coulon (Universite de Nantes)
    # Last update : 22 march 2018
    #
    # Arguments:
    # historique   The composition of groups for each dlim and each initial elements
    #
    #===================================================================================

    historique <- as.matrix(historique)
    #variable securite
    OK <- 1

    if(length(historique[,1]) <= 1){
      cat("No Cluster")
    }
    ecart <- matrix(ncol = length(historique[,1]), nrow = length(historique[,1]))
    #On parcours le tableau pour tous les j on les compare a tous les i (compare tous les elt entre eux)
    for (j in 1:length(historique[,1])) {
      for (i in 1:length(historique[,1])) {
        if(i > j){
          v <- historique[i,]==historique[j,]
          somme <- sum(v %in% FALSE)
          ecart[i,j] <- somme
        }

      }
    }
    ecart <- replace(ecart,is.na(ecart),Inf)
    ecart
  }

IndiceMoran <- function(hist,data) {
  out <- list()
  X <- data[,1]
  Y <- data[,2]
  Z <- data[,3]

  n <- length(hist[1,])
  for (i in 1:n){
    #recup compo of one partition i
    compo <- hist[,i]

    nb <- length(unique(compo))

    for(j in 1:nb){
      subX <- X[compo==j]
      subY <- Y[compo==j]
      subZ <- Z[compo==j]

      if(length(subX) > 3){

        if(length(unique(subZ)) > 2 && length(unique(subX)) > 2 && length(unique(subY)) > 2){
          data.dists <- as.matrix(dist(cbind(subX, subY)))
          Moran <- Moran.I(subZ, data.dists)
        }else{
          Moran <- list(observed= NA, p.value = NA)
        }

      }else{
        Moran <- list(observed= NA, p.value = NA)
      }

      if(!exists("M.I")){
        M.I <- data.frame(grp=c(Moran$observed,Moran$p.value),row.names = c("Moran.I","p.value"))
      }else{
        M.I <- cbind(M.I,grp=c(Moran$observed,Moran$p.value))
      }
    }
    out <- c(out,list(M.I))
    rm(M.I)
  }
  out
}

load(file = "./data/datarcheo.RData")
load(file = "./data/datacancer.RData")

#----------------------------------------------------------------------------------------------#
###########################                - SERVER -                ###########################
#----------------------------------------------------------------------------------------------#
server <- function(input, output, session) {

  dlmap <- reactiveValues(map = 0)
  values <- reactiveValues()
  plot.dend <- reactiveValues(main=NULL, add=NULL, add2=NULL)
  values[["path"]] <- "wait.png"
  values[["multi"]] <- "wait"
  values[["lonlat"]] <- TRUE
  values[["silhouetteplot"]] <- NULL
  values[["labelall"]] <- NULL
  values[["isLabel"]] <- TRUE

  #------------------------------------------------------#
  # Download
  #------------------------------------------------------#
  # output$map.png <- downloadHandler(
  #   filename = "map.png",
  #   content = function(file) {
  #     showModal(waitModal())
  #     mapshot(dlmap$map, file = file)
  #     removeModal()
  #   }
  # )

  #------------------------------------------------------#
  #  modalDialog zdata ask
  #------------------------------------------------------#
  # Return the UI for a modal dialog with data selection input. If 'failed' is
  # TRUE, then display a message that the previous value was invalid.
  dataModal <- function(failed = FALSE, data = 1,multi = FALSE) {
    if(data == 3 && multi=="Univariate"){
      modalDialog(title="datatypemessage",size = c("s"),

                    h4("Data:",style="margin-left:10px;"),hr(style="border-color: #222222;"),
                    span("This classification requires geographical coordinates and a third variable: var."),
                    br(),
                    span("We use var to weight the distance and var must be a frequency or density and always positive."),
                    br(),br(),
                    radioButtons("z_select", "Choose data type of your var",
                                 choices = c(positive = "un",
                                             other = "de"),
                                 selected = "de"),

                  if (failed)
                    div(tags$b("You have to choose", style = "color: red;")),
                  footer = tagList(
                    modalButton("Cancel"),
                    actionButton("dismiss_modal","OK")
                  )
      )
    }else{
      modalDialog(title="datatypemessage",size = c("s"),
                  h4("Are you sure ? "),hr(style="border-color: #222222;"),
                  footer = tagList(
                    modalButton("Cancel"),
                    actionButton("dismiss_modal","OK")
                  )
      )
    }
  }
  #------------------------------------------------------#
  #  modalDialog cuttree message
  #------------------------------------------------------#
  # Return the UI for a modal dialog with data selection input. If 'failed' is
  # TRUE, then display a message that the previous value was invalid.
  cuttreeModal <- function(failed = FALSE) {
    modalDialog(title="cuttreemessage",size = c("l"),
                h4("Cut the tree:",style="margin-left:10px;"),hr(style="border-color: #222222;"),
                span("You can clic on the dendrogram "),img(src="tree.png",height = 23, width = 26),span(" in order to select the dlim where to cut the tree."),
                span("The map and the silhouette plot will be automatically updated."),

                if (failed)
                  div(tags$b("You have to choose", style = "color: red;")),
                footer = tagList(
                  modalButton("OK")
                )
    )
  }
  #------------------------------------------------------#
  #  modalDialog wait message
  #------------------------------------------------------#
  waitModal <- function(failed = FALSE) {
    modalDialog(title="waitmessage",size = c("s"),
                h4("Wait . . . ",style="margin-left:100px;"),hr(style="border-color: #222222;"),
                footer = NULL
    )
  }

  #------------------------------------------------------#
  #  modalDialog upload
  #------------------------------------------------------#
  uploadModal <- function(failed = FALSE) {
    modalDialog(title = "upload",size=c("l"),
                h4("Upload file"),
                # Input: Select a file ----
                fileInput("file1", span(icon("import", lib = "glyphicon"),"Choose CSV File"),
                          multiple = FALSE,
                          accept = c("text/csv",
                                     "text/comma-separated-values,text/plain",
                                     ".csv")
                ),

                # Horizontal line ----
                hr(style="border-color: #222222;"),
                actionButton(inputId = "toggle", span(icon("cogs",lib = "font-awesome"),"Settings")),
                fluidRow(id="myBox",
                  column(width = 5,
                         span(strong("Non-geographic variables.")),
                         sliderTextInput("UniMulti","Your case is:",choices=c("Univariate","Multivariate"),width = "100px"),
                         hr(style="border-color: #C6C0AA;"),
                         sliderInput("nb_var","How many (non-geographic) variables do you have?",min=2,max=6,value=2)
                  ),
                  column(width = 3,
                         # Input: Select Label ----
                         span(strong("Your data have:")),
                         checkboxInput("label", "Label", TRUE),
                         checkboxInput("lonlat", "Longitude\nLatitude", TRUE),
                         hr(style="border-color: #C6C0AA;"),
                         # Input: Select separator ----
                         radioButtons("sep", "Separator",
                                      choices = c(Comma = ",",
                                                  Semicolon = ";",
                                                  Tab = "\t",
                                                  Space = " "),
                                      selected = ";")
                  ),
                  column(width = 3,
                         # Input: Checkbox if file has header ----
                         span(strong("Header")),
                         checkboxInput("header", "Header", TRUE),
                         hr(style="border-color: #C6C0AA;"),
                         # Input: Select quotes ----input$dec
                         radioButtons("quote", "Quote",
                                      choices = c(None = "",
                                                  "Double Quote" = '"',
                                                  "Single Quote" = "'"),
                                      selected = '"'),
                         radioButtons("dec", "Decimal",
                                      choices = c("Dot" = '.',
                                                  "Comma" = ","),
                                      selected = '.')
                  )
                ),
                # Horizontal line ----
                hr(style="border-color: #C6C0AA;"),
                fluidRow(
                  column(width = 8
                  ),
                  column(width = 2,
                         h4("Match:")
                  ),
                  column(width = 1,
                         imageOutput("image2",width = "40px", height = "40px")
                  )
                ),
                h4("How your data must be organized"),
                hr(style="border-color: #222222;"),
                tableOutput("table"),
                # Output: Data file ----
                h4("Your data"),hr(style="border-color: #222222;"),
                tableOutput("contents"),
                if (failed)
                  div(tags$b("Your data must be match with the example above ! Change parameter or file.", style = "color: red;")),
                footer = tagList(
                  modalButton("Cancel"),
                  actionButton("okUL", "OK")
                )
    )
  }
  #------------------------------------------------------#
  #  Uploading File
  #------------------------------------------------------#
  output$image2 <- renderImage({

    if (values[["path"]] == "yes.jpg") {
      return(list(
        src = "www/yes.jpg",
        width="40px",
        contentType = "image/jpg",
        alt = "yes"
      ))
    }else{
      return(list(
        src = "www/no.jpg",
        width="40px",
        contentType = "image/jpg",
        alt = "no"
      ))
    }

  }, deleteFile = FALSE)

  observe({
    if(!is.null(input$UniMulti)){
      if(input$UniMulti == "Univariate"){
        shinyjs::hide("nb_var")
      }else{
        shinyjs::show("nb_var")
      }
    }
  })

  observeEvent(input$UL,{
    showModal(uploadModal())
    shinyjs::hide("myBox")
    shinyjs::hide("nb_var")
  })

  observe({
    if(!is.null(values[["test"]])){
      if(values[["test"]] == length(values[["DF"]][1,])){
        values[["path"]] <- "yes.jpg"
      }else{
        values[["path"]] <- "no.jpg"
      }
    }
  })

  observeEvent(input$okUL,{
    if(values[["test"]] == length(values[["DF"]][1,])){
      values[["nb_var"]] <- input$nb_var
      values[["multi"]] <- input$UniMulti
      values[["lonlat"]] <- input$lonlat

      data <- values[["DF"]]
      multi <- values[["multi"]]
      lonlat <- values[["lonlat"]]

      if(!is.numeric(data[,1])){

        showNotification("coordinates must be numeric! (Check the decimal character)",type = "warning", duration=15)
        showModal(uploadModal(failed = TRUE))
        shinyjs::hide("nb_var")
      }else{
        if(!is.numeric(data[,2])){

          showNotification("coordinates must be numeric! (Check the decimal character)",type = "warning", duration=15)
          showModal(uploadModal(failed = TRUE))
          shinyjs::hide("nb_var")
        }else{
          if(!is.numeric(data[,3])){

            showNotification("variable must be numeric! (Check the decimal character)",type = "warning", duration=15)
            showModal(uploadModal(failed = TRUE))
            shinyjs::hide("nb_var")
          }else{
            #change data select
            updateSelectInput(session, "select", selected = 3)
            removeModal()
          }

        }


      }



    }else{
      showModal(uploadModal(failed = TRUE))
      shinyjs::hide("nb_var")
    }
  })

  output$contents <- renderTable({

    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    req(input$file1)
    inFile <- input$file1

    if(is.null(inFile))
      return(NULL)

    tryCatch({
        df <- read.csv(input$file1$datapath,
                       header = input$header,
                       sep = input$sep,
                       dec = input$dec,
                       quote = input$quote)
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
    })

    if (!is.null(inFile)){
      DF <- df
      values[["DF"]] <- DF
    }

    return(head(df))

  })

  output$table <- renderTable(
    if(input$UniMulti == "Univariate"){
      if(input$label){
        values[["test"]] <- 4
        values[["isLabel"]] <- TRUE
        data.frame(coord.1="X",coord.2="Y",var="var",label_optional="label")
      }else{
        values[["test"]] <- 3
        values[["isLabel"]] <- FALSE
        data.frame(coord.1="X",coord.2="Y",var="var")
      }

    }else{
      if(input$nb_var == 2){
        if(input$label){
          values[["test"]] <- 5
          values[["isLabel"]] <- TRUE
          data.frame(coord.1="X",coord.2="Y",var.1="var1",var.2="var2",label_optional="label")
        }else{
          values[["test"]] <- 4
          values[["isLabel"]] <- FALSE
          data.frame(coord.1="X",coord.2="Y",var.1="var1",var.2="var2")
        }
      }else{
        if(input$nb_var == 3){
          if(input$label){
            values[["test"]] <- 6
            values[["isLabel"]] <- TRUE
            data.frame(coord.1="X",coord.2="Y",var.1="var1",var.2="var2",var.3="var3",label_optional="label")
          }else{
            values[["test"]] <- 5
            values[["isLabel"]] <- FALSE
            data.frame(coord.1="X",coord.2="Y",var.1="var1",var.2="var2",var.3="var3")
          }
        }else{
          if(input$nb_var == 4){
            if(input$label){
              values[["test"]] <- 7
              values[["isLabel"]] <- TRUE
              data.frame(coord.1="X",coord.2="Y",var.1="var1",var.2="var2",var.3="var3",var.4="var4",label_optional="label")
            }else{
              values[["test"]] <- 6
              values[["isLabel"]] <- FALSE
              data.frame(coord.1="X",coord.2="Y",var.1="var1",var.2="var2",var.3="var3",var.4="var4")
            }
          }else{
            if(input$nb_var == 5){
              if(input$label){
                values[["test"]] <- 8
                values[["isLabel"]] <- TRUE
                data.frame(coord.1="X",coord.2="Y",var.1="var1",var.2="var2",var.3="var3",var.4="var4",var.5="var5",label_optional="label")
              }else{
                values[["test"]] <- 7
                values[["isLabel"]] <- FALSE
                data.frame(coord.1="X",coord.2="Y",var.1="var1",var.2="var2",var.3="var3",var.4="var4",var.5="var5")
              }
            }else{
              if(input$label){
                values[["test"]] <- 9
                values[["isLabel"]] <- TRUE
                data.frame(coord.1="X",coord.2="Y",var.1="var1",var.2="var2",var.3="var3",var.4="var4",var.5="var5",var.6="var6",label="label")
              }else{
                values[["test"]] <- 8
                values[["isLabel"]] <- FALSE
                data.frame(coord.1="X",coord.2="Y",var.1="var1",var.2="var2",var.3="var3",var.4="var4",var.5="var5",var.6="var6")
              }

            }

          }

        }

      }

    }
  )

  observe({
    if(!values[["isLabel"]]){
      shinyjs::hide(id = "Plabel")
    }else{
      shinyjs::show(id = "Plabel")
    }
  })

  observeEvent(input$toggle, {
    if(input$toggle %% 2 == 0){
      shinyjs::hide(id = "myBox")
    }else{
      shinyjs::show(id = "myBox")
    }
  })

  #------------------------------------------------------#
  #  MapClust
  #------------------------------------------------------#
  observe({
    if(input$select==3){
      if (is.null(values[["DF"]])){
        showNotification("Oops ! you forgot a little thing ! To plot your data, you need to follow a few things:
                       You must import your own data! Re Run after import!",type = "warning", duration = 100)
      }
      if(input$UniMulti == "Univariate"){
        if(input$label){
          values[["isLabel"]] <- TRUE
        }else{
          values[["isLabel"]] <- FALSE
        }

      }else{
        if(input$nb_var == 2){
          if(input$label){
            values[["isLabel"]] <- TRUE
          }else{
            values[["isLabel"]] <- FALSE
          }
        }else{
          if(input$nb_var == 3){
            if(input$label){
              values[["isLabel"]] <- TRUE
            }else{
              values[["isLabel"]] <- FALSE
            }
          }else{
            if(input$nb_var == 4){
              if(input$label){
                values[["isLabel"]] <- TRUE
              }else{
                values[["isLabel"]] <- FALSE
              }
            }else{
              if(input$nb_var == 5){
                if(input$label){
                  values[["isLabel"]] <- TRUE
                }else{
                  values[["isLabel"]] <- FALSE
                }
              }else{
                if(input$label){
                  values[["isLabel"]] <- TRUE
                }else{
                  values[["isLabel"]] <- FALSE
                }

              }

            }

          }

        }

      }
    }else{
      values[["isLabel"]] <- TRUE
    }
  })

  ##select data and run calcul
  observeEvent(input$Run,{
    if (input$select==3 && is.null(values[["DF"]])){
      showNotification("Oops ! you forgot a little thing ! To plot your data, you need to follow a few things:
                       You must import your own data and after press again the run button!",type = "error", duration = 100)
      showModal(uploadModal())
    }else{
      if (input$select==3){showModal(dataModal(data=3,multi=input$UniMulti))}
      if (input$select==2){showModal(dataModal(data=2))}
      if (input$select==1){showModal(dataModal(data=1))}
    }
  })

  observeEvent(input$dismiss_modal,{
    removeModal()
    showModal(waitModal())

    if(input$select==1){
      data <- datarcheo
      ctrl <- "de"
      multi <- "Univariate"
      lonlat <- TRUE
    }
    if(input$select==2){
      data <- datacancer
      ctrl <- "un"
      multi <- "Univariate"
      lonlat <- TRUE
    }
    OK <- TRUE
    if(input$select==3){
      if (is.null(values[["DF"]])){
        showNotification("There is no upload data. We select datarcheo by default",type = "warning", duration = 100)
        data <- datarcheo
      }else{
        data <- values[["DF"]]
        multi <- values[["multi"]]
        lonlat <- values[["lonlat"]]

        if(!is.numeric(data[,1])){
          OK <- FALSE
          showNotification("coordinates must be numeric!",type = "warning", duration=15)
        }
        if(!is.numeric(data[,2])){
          OK <- FALSE
          showNotification("coordinates must be numeric!",type = "warning", duration=15)
        }

        if(multi == "Multivariate"){
          for(i in 3:(2+values[["nb_var"]])){
            if(!is.numeric(data[,i])){
              OK <- FALSE
              showNotification("variable must be numeric!",type = "warning", duration=15)
            }
          }
        }else{
          if(!is.numeric(data[,3])){
            OK <- FALSE
            showNotification("variable must be numeric!",type = "warning", duration=15)
          }
        }
      }
    }
    if(OK == TRUE){

      values[["multi2"]] <- multi
      if(input$Plabel==T){Plabel = T}
      if(input$Plabel==F){Plabel = F}
      if(input$extenddend==T){extenddend = T}
      if(input$extenddend==F){extenddend = F}
      if(input$select==3 && multi == "Univariate"){
        if(input$z_select=="un"){ctrl="un"}else{ctrl="de"}
      }
      values[["data"]] <- data

      values[["index"]] <- NULL
      binf = 0
      bsup = max(dist(data[,1:2]))
      pas = (max(dist(data[,1:2]))/20)
      plot.dend$add <- NULL
      plot.dend$add2 <- NULL

      progress1 <- Progress$new(session, min = 0, max = 2)
      on.exit(progress1$close())
      progress1$set(message = 'Initialisation and calcul of the distance matrix ..')
      progress1$set(value = 1)

      progress2 <- Progress$new(session,style = "notification", min = 0, max = round((bsup - binf) / pas,0))
      on.exit(progress2$close())
      progress2$set(message = 'Calculation in progress',
                    detail = 'This may take a while...')


      #Si Z negatif : reduction
      if(multi == "Multivariate"){
        Z <- try(Z <- estimateDensityMulti(data[,3:(values[["nb_var"]]+2)]))
        values[["mdata"]] <- data[,3:(values[["nb_var"]]+2)]
        if(!is.numeric(Z)){
          Z <- transfo_z(data[,3])
          values[["multi2"]] <- "Univariate"
          multi <- "Univariate"
          showNotification("Error in chol.default(S) :
    the leading minor of order 1 is not positive definite.",type="error")
          showNotification("There was a error during the kernel density estimation. We consider only the first variable for the rest of the process.",type="warning",duration=40)
          saveZ <- data[,3]
          data[,3] <- Z
        }else{
          saveZ <- rep(7,length(data[,3]))
        }

      }else{
        Z <- data[,3]
        if(!ctrl=="un"){
          Z <- transfo_z(Z)
        }else{
          if(min(Z)<0){
            Z <- transfo_z(Z)
            showNotification("The var data have been changed due to negative value in third column.",type="error")
          }
        }
        saveZ <- data[,3]
        data[,3] <- Z
      }


      #on enregistre dans l'element renvoyer a l'utilisateur les parametres utilises
      X <- data[,1]
      Y <- data[,2]
      L <- data[,length(data[1,])]
      data <- data.frame(X,Y,Z)
      xg <- c()
      groupes=1
      dlim=1
      xend=1
      yend=1

      res <- NULL
      values[["Size"]] <- saveZ
      values[["C"]] <- ifelse(saveZ<0,"2, 5","1")
      W = rep(1 , length(Y))


      if (bsup <= binf) {
        showNotification("The bsup value must be higher than binf",type = "error", duration = 10)
        stop("The bsup value must be higher than binf")
      }

      label <- as.character(L)

      id <- c(1:length(data[,1]))

      if(values[["isLabel"]]){
        values[["labelall"]] <- label
      }else{
        values[["labelall"]] <- id
        label <- as.character(id)
      }


      data=cbind(data,id)
      data=as.data.frame(data)

      if( bsup > max(dist(data[,1:2])) ) {
        cat("The maximal value of dlim is: ",round(max(dist(data[,1:2])),2))
        cat("bsup have been changed in order to not run high dlim (refers to data)")
        bsup <- max(dist(data[,1:2]))
      }

      #remise a zero des objects
      names <- c()
      #initialisation matrice historique des groupes
      historique <- matrix(ncol=0, nrow=length(X))

      cat("\nInitialisation and calcul of the distance matrix ..")
      DiMatrix <- dist_m(data)
      progress1$set(value = 2)
      cat(". Done\n")
      cat("Running Progress:\n")
      values[["DiMatrix"]] <- DiMatrix
      #creation barre de progression(pb)
      pb <- txtProgressBar( min = 0, max = round((bsup - binf) / pas,0) , style = 3 )
      progress1$close()

      pre_k=0
      k = 1
      pass =0
      grp=c(1)
      data= cbind(data,grp)
      dlim <- bsup
      #running spatialpatches (WOILLEZ et al., 2009)
      Pch <- spatialpatches2(X , Y , Z , w = rep(1 , length(Y)) , Lim.D = dlim)
      grp1 <- Pch$n
      names <- c(round(dlim,2))
      historique <- cbind(historique,grp1)

      data$grp=grp

      if (!(length(unique(data$grp)) <= 1)) {
        showNotification("The maximal value of dlim is not enough to create one patch at the top of the tree.
                         The groups may be merged too early. Try higher bsup.",type = "message", duration = 10)
      }

      while ( dlim >= binf) {

        #update progress bar
        Sys.sleep(0.1)
        setTxtProgressBar(pb, k)
        progress2$set(value = k)

        #if k have not been increment
        if(k == pre_k){
          dlim <- nextdlim
        }else{
          dlim = bsup - k * pas
          #save value of k for this loop. need for chack next loop time
          pre_k <- k
        }

        if (round(dlim,2) > 0) {

          nb_grp=length(unique(data$grp))
          if(nb_grp>1){
            #More than one patch so we separate the data
            grp_tmp <- c()

            for (a in 1:nb_grp){
              # a take the different value corresponding to  the numero of cluster

              # eval permit to run the text into (in order to run, in the loop, data1 and data2 and not dataa)
              # parse permit to pass text to expression in eval(...)
              # Separtion of data into sub dataset (data1, data2, ...) order by the cluster vector Pch$n
              eval(parse(text=paste0("data",a,"=data[data$grp==",a,",]")))

              # Run spatialpatche for each of the new dataset (data1, data2, ...)
              # if number of element i is sup to 2
              if(eval(parse(text=paste0("length(data",a,"[,1])>3")))){
                eval(parse(text=paste0("Pch",a," <- spatialpatches2(data",a,"[,1],data",a,"[,2],data",a,"[,3],w = rep(1 , length(data",a,"[,2])),Lim.D = dlim)")))
              }else{
                eval(parse(text=paste0("Pch",a,"<-NULL")))
                eval(parse(text=paste0("Pch",a,"$n<-rep(1,length(data",a,"[,1]))")))
              }

              #recuperation of the cluster vector for each sub dataset (data1, data2, ...)
              eval(parse(text=paste0("data",a,"$grp=Pch",a,"$n")))
              #On ajoute a en dizaine au nombre pour ne pas generer de doublons entre les differents a
              eval(parse(text=paste0("data",a,"$grp=as.numeric(paste(",a,",data",a,"$grp,sep=''))")))

              #Si 3 groupes.
              if( eval(parse(text=paste0("length(unique(data",a,"$grp)) == 3"))) ){
                eval(parse(text=paste0("XX <- data",a,"[,1]")))
                eval(parse(text=paste0("YY <- data",a,"[,2]")))
                eval(parse(text=paste0("ZZ <- data",a,"[,3]")))

                eval(parse(text=paste0("WW <- rep(1 , length(data",a,"[,2]))")))
                eval(parse(text=paste0("xg <- tapply(WW*ZZ*XX,data",a,"$grp,sum,na.rm=T)/tapply(WW*ZZ,data",a,"$grp,sum,na.rm=T)")))
                eval(parse(text=paste0("yg <- tapply(WW*ZZ*YY,data",a,"$grp,sum,na.rm=T)/tapply(WW*ZZ,data",a,"$grp,sum,na.rm=T)")))
                diff1_2 <- abs(xg[1]-xg[2])
                diff1_3 <- abs(xg[1]-xg[3])
                diff2_3 <- abs(xg[2]-xg[3])
                if(diff2_3 > diff1_2 && diff2_3 > diff1_3){
                  eval(parse(text=paste0("data",a,"$grp <- ifelse(data",a,"$grp==",a,"2,",a,"2,",a,"1)")))
                }else{
                  eval(parse(text=paste0("data",a,"$grp <- ifelse(data",a,"$grp==",a,"1,",a,"1,",a,"2)")))
                }
                nextdlim <- dlim - (pas/3)
                #pass take value 1 in order to re rsun one time the loop witjhout change k and replace dlim by nextdlim.
                #This nextdlim permit to separate in two times (one dlim and one nextdlim near to dlim) the 3 patchs.
                pass <- 1
              }
              #Si plus de 3 groupes.
              if( eval(parse(text=paste0("length(unique(data",a,"$grp)) >= 4"))) ){
                eval(parse(text=paste0("XX <- data",a,"[,1]")))
                eval(parse(text=paste0("YY <- data",a,"[,2]")))
                eval(parse(text=paste0("ZZ <- data",a,"[,3]")))

                eval(parse(text=paste0("WW <- rep(1 , length(data",a,"[,2]))")))
                eval(parse(text=paste0("xg <- tapply(WW*ZZ*XX,data",a,"$grp,sum,na.rm=T)/tapply(WW*ZZ,data",a,"$grp,sum,na.rm=T)")))
                eval(parse(text=paste0("yg <- tapply(WW*ZZ*YY,data",a,"$grp,sum,na.rm=T)/tapply(WW*ZZ,data",a,"$grp,sum,na.rm=T)")))
                diff1_2 <- abs(xg[1]-xg[2])
                diff1_3 <- abs(xg[1]-xg[3])
                diff2_3 <- abs(xg[2]-xg[3])
                if(diff2_3 > diff1_2 && diff2_3 > diff1_3){
                  eval(parse(text=paste0("data",a,"$grp <- ifelse(data",a,"$grp==",a,"2,",a,"2,",a,"1)")))
                }else{
                  eval(parse(text=paste0("data",a,"$grp <- ifelse(data",a,"$grp==",a,"1,",a,"1,",a,"2)")))
                }
              }
              #fusion of the cluster vectors
              eval(parse(text=paste0( "grp_tmp <- rbind(grp_tmp,data",a,")" )))
            }

            grp_tmp <- grp_tmp[order(grp_tmp$id,decreasing=F), ]

            #change the cluster names (11, 21, ...) into new (1, 2, ...) more simple
            namegrp<-c()
            namegrp<-as.factor(grp_tmp$grp)
            levels(namegrp) <- c(1:length(levels(namegrp)))
            data$grp<-as.numeric(namegrp)

          }else{
            # only one patch in all so we continu on all the data
            Pch <- spatialpatches2(X , Y , Z , w = rep(1 , length(Y)) , Lim.D = dlim)
            grp <- Pch$n
            data$grp=grp
          }

          if(round(dlim,1) > 0){
            j=length(historique[1,])
            #compare last and new compo
            v <- data$grp==historique[,j]
            somme <- sum(v %in% FALSE)
            if(somme != 0){
              #if data$grp is different than the last compo hist
              #add the dlim to name and the hist
              names <- c(names , round(dlim,2))
              historique <- cbind(historique , data$grp)
            }else{
              if(round(dlim,2) < names[length(names)]){
                #IF no difference : replace the last dlim name by the new
                names[length(names)] <- round(dlim,2)
              }
            }
          }
        }
        if(pass == 1){
          #if pass egal to 1 don't increment K and re run one time the while loop. pass egal 0 for next time of loop
          pass <- 0
        }else{
          k = k + 1
        }
      }#end while

      dlim <- as.numeric(colnames(historique)[length(colnames(historique))])

      historique <- cbind(grp1,historique)
      names <- c(round(bsup,2),names)
      historique <- as.matrix(historique)
      colnames(historique) <- names
      res$dlim <- names


      eval_class <- "test"
      output$log <- renderPrint({
        eval_class
      })

      #silhouette
      Average_Sil_Width <- c()
      for(i in 1:length(historique[1,])){
        cluster <- historique[,i]
        if(length(unique(cluster))>1){
          sil <- silhouette(cluster,DiMatrix)
          avg <- mean(sil[(length(cluster)*2):(length(cluster)*3)])
          Average_Sil_Width <- c(Average_Sil_Width,avg)
        }
      }
      number <- apply(historique,2,unique)
      number <- lapply(number,length)

      Number_of_groups <- as.vector(unlist(number))[number>1]

      Average_Sil_Width <- rbind(Number_of_groups,Average_Sil_Width)
      values[["sil_width"]] <- Average_Sil_Width

      historique <- cbind(id,historique)
      hist<-c()

      cpy_historique <- historique
      cpy_historique<-cbind(label,cpy_historique)

      #On trie l'historique par la colonne du dlim minimal
      historique<-historique[order(historique[,length(historique[1,])],decreasing=F), ]
      cpy_historique<-cpy_historique[order(as.numeric(cpy_historique[,length(cpy_historique[1,])]),decreasing=F), ]


      label_ini_grp <- c()
      label_tmp <- c()
      label_avant <- c()
      label_grp <- c()

      if(extenddend){

        for(b in 2:(length(cpy_historique[,1]))){
          c=b-1
          if(!(FALSE %in% (historique[b,2:length(historique[1,])]==historique[c,2:length(historique[1,])]))){
            #Doublon
            label_tmp <- as.character(cpy_historique[b,1])
            label_avant <- as.character(cpy_historique[c,1])
            if(is.null(label_grp)){
              label_grp <- c(label_avant,label_tmp)
            }else{
              label_grp <- c(label_grp,label_avant,label_tmp)
            }
          }else{
            #Pas doublon
            hist<-rbind(hist,historique[c,])
            label_grp <- c(label_grp,as.character(cpy_historique[c,1]))
            label_grp <- paste(unique(label_grp), collapse = ".")
            label_ini_grp<-c(label_ini_grp,label_grp)
            label_grp <- c()
          }
        }

        #add the last elt
        hist<-rbind(hist,historique[length(historique[,1]),])
        label_grp <- paste(unique(label_grp), collapse = ".")
        label_ini_grp<-c(label_ini_grp,label_grp)

      }else{
        #----------------
        for(b in 2:(length(cpy_historique[,1]))){
          c=b-1
          label_grp <- c(label_grp,as.character(cpy_historique[c,1]))
          label_ini_grp<-c(label_ini_grp,label_grp)
          label_grp <- c()
        }
        hist <- historique
        #add last elt
        label_grp <- unique(as.character(cpy_historique[length(historique[,1]),1]))
        label_ini_grp <- c(label_ini_grp,label_grp)
        #----------------
      }


      #We cut the label colunms (Not use here because we work on the entire cluster the label and id of
      # each members of one cluster are still save in data)
      hist<-as.matrix(hist)
      hist<-hist[,2:length(historique[1,])]

      ecart=ecartmatrice(hist)
      values[["hist"]] <- hist
      #if un seul groupe pas de div gerer ce cas

      ##Les segments
      #creation colonne x
      x <- c()
      for (b in 1:length(hist[,1])) {
        x <- c(x , b)
      }
      hist<-as.matrix(hist)
      #creation colonne y
      y <- rep(binf , length(hist[,1]))

      #creation d'une matrice ou on stockera les coordonnees pour les segments du drendrogramme
      coord <- cbind(x , y)

      #initialisation des segments
      segment <- c()

      #boucle se repetant (nombre de grp initiaux -1) fois 3 pour ecrire les segments
      #Pour 6 grps il faut (6-1)*3=15 segments. n grp il faut (n-1)*3 segments.
      for (a in 1:(length(hist[,1])-1)) {

        #on recherche min (premi?re fusion)
        mini <- min(ecart)
        #On recherche les indices
        index_vect <- which.min(ecart)
        index <- arrayInd(index_vect,c(length(ecart[,1]),length(ecart[,1])))
        grp1 <- index[1]
        grp2 <- index[2]

        #on recupere dlim dans l'objet res$dlim qui contient tous les dlim utilises. length(..)-mini pour
        # partir de la fin et remonter. mini est le nombre de dlim a passer vers la gauche.

        if(length(res$dlim)-mini!=0){
          dlim_tmp <- res$dlim[length(res$dlim)-mini]

          #on construit 3 segments (deux qui monte et un qui fusionne branche horizontal)

          # segments verticaux
          # x = l'abscisse du tab coord pour grp correspondant  y = coord[grp_,2]    xend = x    yend = dlim_tmp
          ######################################################################################################
          segment <- rbind(segment , c(coord[grp1,1] ,     coord[grp1,2] ,  coord[grp1,1] ,         dlim_tmp))
          segment <- rbind(segment , c(coord[grp2,1] ,     coord[grp2,2] ,  coord[grp2,1] ,         dlim_tmp))

          # segment horizontal
          ##################     x = coord[grp1,1]     y = dlim_tmp     xend = coord[grp2,1]    yend = dlim_tmp
          ######################################################################################################
          segment <- rbind(segment,c(coord[grp1,1] ,       dlim_tmp ,       coord[grp2,1] ,         dlim_tmp))


          ## on sauvegarde la hauteur (y / dlim) a laquelle les grps sont rendus
          #MAJ valeur de y (dlim) grp1
          coord[grp1,2] <- dlim_tmp
          #MAJ ligne grp1, valeur en X du nouveau grp (fusion) moyenne des deux
          coord[grp1,1] <- (coord[grp1,1]+coord[grp2,1])/2
          #on marque la ligne du grp2 en mettant x=0, pour ne plus la prendre en compte dans prochaine
          # fusion puisque ce grp a fusionne avec grp1
          coord[grp2,1] <- 0
          #on elimine la ligne/colonne du grp2 dans les valeurs (historique grp)
          ecart[,grp2] <- Inf
          ecart[grp2,] <- Inf
        }
      }
      groupes=1
      dlim=1
      xend=1
      yend=1
      #on rajoute un dernier segment pour la racine de l'arbre
      segment <- rbind(segment , c(max(coord[,1]) , max(coord[,2]) , max(coord[,1]) , (max(res$dlim)+pas)))
      segment <- data.frame(segment)
      #on prepare tableau segment
      names(segment) <- c("groupes" , "dlim" , "xend" , "yend")
      row.names(segment) <- NULL


      #on creer un tableau pour les noms (implant au bout des branches y=0)
      label <- matrix(ncol = 3 , nrow = length(hist[,1]))
      for (a in 1:length(hist[,1])) {
        label[a,1] <- a
        label[a,2] <- 0
        label[a,3] <- label_ini_grp[a]
      }
      label <- data.frame(label)

      values[["label"]] <- label
      values[["dlim"]]<-res$dlim

      #dendrogramme affiche avec indication pour le clic de cut.
      plot.dend$main <- ggplot(segment , col = 'white') +
        geom_segment(data = segment,aes(x = groupes , y = dlim , xend = xend , yend = yend )) +
        ylim(0 , max(res$dlim) + pas)
      #Title subtitle caption
      plot.dend$main <- plot.dend$main + labs(title = "Hierarchical Clustering" , subtitle = "Spatialpatches")
      #axes x
      plot.dend$main <- plot.dend$main + scale_x_continuous(breaks = c())

      plot.dend$main <- plot.dend$main + theme(
        panel.grid = element_line(colour = "grey")
      )

      historique<-historique[order(historique[,1],decreasing=F), ]
      #WSS
      Sum_of_Square <- WSS(historique[,-1],data)
      Within_Sum_of_Square <- c()
      for (i in Sum_of_Square) {
        if(length(i)>1){
          Within_Sum_of_Square <- c(Within_Sum_of_Square,sum(i))
        }
      }

      WSS2 <- c()
      for (i in Sum_of_Square) {
        WSS2 <- c(WSS2,sum(i))
      }
      values[["WSS2"]] <- WSS2

      Within_Sum_of_Square <- rbind(Number_of_groups,Within_Sum_of_Square)
      values[["WSS"]] <- Within_Sum_of_Square

      if(multi == "Univariate"){
        Moran <- IndiceMoran(historique[,-1],data)
        values[["Moran"]] <- Moran
      }


      values[["historique"]] <- historique
      values[["data"]] <- data

      output$distPlot <- renderPlot({
        values[["dendplot"]] <- plot.dend$main + plot.dend$add + plot.dend$add2
        plot.dend$main + plot.dend$add + plot.dend$add2
      })

      output$avgsilPlot <- renderPlot({
        #plot eval sil
        dataplot <- t(values[["sil_width"]])
        dataplot <- as.data.frame(dataplot)
        colnames(dataplot) <- c("number_of_groups","Average_Sil_Width")
        minisil <- min(dataplot$Average_Sil_Width)
        if(minisil > 0 ){
          minisil <- 0
        }else{
          minisil <- -1
        }
        plot <- ggplot(data=dataplot, aes(x=number_of_groups, y=Average_Sil_Width))
        plot <- plot + geom_line(linetype = 3)
        plot <- plot +   geom_point(aes(color=Average_Sil_Width), size=3) + scale_colour_gradientn(colours=c("brown1","chartreuse3"))
        plot <- plot +   scale_x_continuous(breaks = dataplot$number_of_groups,
                                            labels = dataplot$number_of_groups )
        plot <- plot +   labs(title = "Average Silhouette Width" , subtitle = "Partition evaluation")
        plot <- plot + ylim(minisil,1)
        plot <- plot + theme(
          legend.position="bottom",
          plot.background = element_rect(
            colour = "#222222",
            linewidth = 1
          )
        )
        values[["avesilplot"]] <- plot
        plot
      })

      output$WSSPlot <- renderPlot({
        #plot WSS
        dataplot <- t(values[["WSS"]])
        dataplot <- as.data.frame(dataplot)
        colnames(dataplot) <- c("number_of_groups","Within_Sum_of_Square")
        maxwws <- max(dataplot$Within_Sum_of_Square)
        plot <- ggplot(data=dataplot, aes(x=number_of_groups, y=Within_Sum_of_Square))
        plot <- plot + geom_line(linetype = 3)
        plot <- plot + geom_point(aes(color=Within_Sum_of_Square), size=3) + scale_colour_gradientn(colours=c("chartreuse3","brown1"))
        plot <- plot + scale_x_continuous(breaks = dataplot$number_of_groups ,
                                          labels = dataplot$number_of_groups )
        plot <- plot + labs(title = "Within Sum of Square" , subtitle = "Partition evaluation")
        plot <- plot + theme(
          legend.position="bottom",
          plot.background = element_rect(
            colour = "#222222",
            linewidth = 1
          )
        )
        values[["WSSplot"]] <- plot
        plot
      })
      removeModal()
      showModal(cuttreeModal())
      #TODO : add download log graph map and all data
    }else{
      removeModal()
    }

  })

  ##dendrogramme click
  output$click_info <- renderPrint({

    if(is.null(values[["index"]])){
      cat("Hierarchical Divise Clustering")
    }else{
      multi <- values[["multi2"]]
      cat(multi, " ")
      nb_var <- values[["nb_var"]]
      cat(nb_var)

      size <- c()
      data <- values[["data"]]
      if(multi=="Univariate"){
        Z <- values[["Size"]]
      }else{
        Z <-values[["mdata"]]
      }
      cluster <- values[["historique"]][,values[["index"]]+1]
      nb <- length(unique(cluster))
      tot <- length(cluster)
      X_mean <- tapply(data[,1],cluster,mean)
      Y_mean <- tapply(data[,2],cluster,mean)
      cat("The selected partition contains",nb,"clusters among",tot,"observations.\ncluster:\t")
      for (k in 1:nb){
        size <- c(size,length(cluster[cluster == k]))
        cat(k,"\t")
      }
      cat("\nSize:\t\t")
      for(y in 1:nb){
        cat(size[y],"\t")
      }
      cat("\nX_mean:\t\t")
      for(y in 1:nb){
        cat(round(X_mean[y],2),"\t")
      }
      cat("\nY_mean:\t\t")
      for(y in 1:nb){
        cat(round(Y_mean[y],2),"\t")
      }
      if(multi == "Univariate"){
      cat("\nvar_mean:\t")
        var_mean <- tapply(Z,cluster,mean)
        for(y in 1:nb){
          cat(round(var_mean[y],2),"\t")
        }
      }else{
        for(i in 1:nb_var){
          var_mean <- tapply(Z[,i],cluster,mean)
          cat("\nvar",i,"_mean:\t",sep="")
          for(y in 1:nb){
            cat(round(var_mean[y],2),"\t")
          }
        }
      }
    }
  })

  output$click_info2 <- renderPrint({

    multi <- values[["multi2"]]

    if(is.null(values[["index"]])){
      cat("Summary")
    }else{
      cat("> summary()")
      cat("\n  -----------------------------------\n  dlim\t  n\tWSS\tsil\tMoran\n  -----------------------------------\n")
      j <- 1
      for (i in 1:length(values[["dlim"]])){
        dlim <- values[["dlim"]][i]
        cluster <- values[["historique"]][,i+1]
        number_grp <- length(unique(cluster))
        if(number_grp == 1){
          sil <- NA
        }else{
          sil <- values[["sil_width"]][2,j]
          j <- j+1
        }
        WSS <- values[["WSS2"]][i]
        if(multi == "Univariate"){
          Moran <- mean(na.omit(unlist(values[["Moran"]][[i]][1,])))
        }
        cat(" ",dlim,"\t")
        cat(" ",number_grp,"\t")
        cat(round(WSS,2),"\t")
        cat(round(sil,2),"\t")
        if(multi == "Univariate"){
          cat(round(Moran,2),"\t")
        }
        cat("\n")
      }
      cat("  -----------------------------------\n")
    }
  })



  observeEvent(input$Plabel,{
    if(input$Plabel){
      label <- values[["label"]]
      if(!is.null(values[["index"]])){
        type <- values[["hist"]][,values[["index"]]]
      }else {
        type <- rep(1,length(values[["hist"]][,1]))
      }
      couleur <- c("slateblue2", "blue", "chartreuse1","darkolivegreen2",
                   "#222222", "darkmagenta", "mediumpurple4",
                   "magenta","cyan3", "plum4","darkorange1", "pink1", "red",
                   "lightgoldenrod", "violetred3", "brown4", "tan1","yellow", "gray50",
                   "#967272","deepskyblue4","salmon", "cyan","palegreen", "green3")
      pal <- colorFactor(couleur,as.factor(type))
      col <- pal(as.factor(type))
      if(!is.null(plot.dend$add)){
        plot.dend$add2 <- geom_text(data = label, aes(x = as.numeric(paste(label[,1])) ,
                                                      y = as.numeric(paste(label[,2])) ,
                                                      label = as.character(label[,3])) ,
                                    fontface = "bold" , nudge_x = -0.12 , colour = col ,
                                    hjust = 0 , vjust = 0.25 , size = 4 , angle = 90)
      }else{
        plot.dend$add2 <- geom_text(data = label, aes(x = as.numeric(paste(label[,1])) ,
                                                      y = as.numeric(paste(label[,2])) ,
                                                      label = as.character(label[,3])) ,
                                    fontface = "bold" , nudge_x = -0.12 , colour = "royalBlue4" ,
                                    hjust = 0 , vjust = 0.25 , size = 4 , angle = 90)
      }
    }else{
      plot.dend$add2 <- NULL
    }

  })


  observeEvent(input$plot_click,{
    diff <- values[["dlim"]]-input$plot_click$y
    #which.max (<=0) return index of the first value <= 0
    if((TRUE %in% (diff <= 0))){
      index <- which.max(diff <= 0)
    }else{
      index <- which.min(diff)
    }
    values[["cutdlim"]] <- values[["dlim"]][index]
    values[["index"]] <- index
    showNotification(paste("cutdlim : ",values[["cutdlim"]]),type = "message",duration = 10)

    plot.dend$add
    plot.dend$add2
    plot.dend$add <- geom_hline(yintercept = values[["cutdlim"]] + 0.02 , col = "dimgrey" , lty = 2 ,lwd = 0.5)
    if (input$Plabel) {
      if(!is.null(values[["index"]])){
        type <- values[["hist"]][,values[["index"]]]
      }else {
        type <- rep(1,length(values[["hist"]][,1]))
      }
      couleur <- c("slateblue2", "blue", "chartreuse1","darkolivegreen2",
                   "#222222", "darkmagenta", "mediumpurple4",
                   "magenta","cyan3", "plum4","darkorange1", "pink1", "red",
                   "lightgoldenrod", "violetred3", "brown4", "tan1","yellow", "gray50",
                   "#967272","deepskyblue4","salmon", "cyan","palegreen", "green3")
      pal <- colorFactor(couleur,as.factor(type))
      col <- pal(as.factor(type))
      label <- values[["label"]]
      plot.dend$add2 <- geom_text(data = label, aes(x = as.numeric(paste(label[,1])) ,
                                                    y = as.numeric(paste(label[,2])) ,
                                                    label = as.character(label[,3])) ,
                                                   fontface = "bold" , nudge_x = -0.12 , colour = col ,
                                                   hjust = 0 , vjust = 0.25 , size = 4 , angle = 90)
    }else{
      plot.dend$add2 <- NULL
    }
  })

  ##leaflet map
  output$mymap <- renderLeaflet({
    #leaflet
    if(input$select == 3){
      lonlat <- values[["lonlat"]]
    }else{
      lonlat <- TRUE
    }
    if(lonlat){
      esri <- structure(c("Esri.WorldStreetMap", "Esri.DeLorme", "Esri.WorldTopoMap",
                          "Esri.WorldImagery", "Esri.WorldTerrain", "Esri.WorldShadedRelief",
                          "Esri.WorldPhysical", "Esri.OceanBasemap", "Esri.NatGeoWorldMap",
                          "Esri.WorldGrayCanvas"),
                        .Names = c("Esri.WorldStreetMap", "Esri.DeLorme", "Esri.WorldTopoMap",
                                   "Esri.WorldImagery", "Esri.WorldTerrain", "Esri.WorldShadedRelief",
                                   "Esri.WorldPhysical", "Esri.OceanBasemap", "Esri.NatGeoWorldMap",
                                   "Esri.WorldGrayCanvas"))
      m <- leaflet()
      for (provider in esri) {
        m <- m %>% addProviderTiles(provider, group = provider)
      }
      m <- m %>%
        addLayersControl(baseGroups = names(esri),
                         options = layersControlOptions(collapsed = TRUE),
                         position = "topright"
        ) %>%
        addMiniMap(tiles = esri[[1]], toggleDisplay = TRUE,minimized = TRUE,
                   position = "bottomright"
        ) %>%
        htmlwidgets::onRender("
                              function(el, x) {
                                var myMap = this;
                                myMap.on('baselayerchange',
                                function (e) {
                                  myMap.minimap.changeLayer(L.tileLayer.provider(e.name));
                                })
                              }"
        )
    }else{
      m <- leaflet()
    }
    dlmap$map <- m
  })

  observe({
    if(input$select == 3){
      lonlat <- values[["lonlat"]]
    }else{
      lonlat <- TRUE
    }
    multi <- values[["multi2"]]
    couleur <- c("slateblue2", "blue", "chartreuse1","darkolivegreen2",
                 "#222222", "darkmagenta", "mediumpurple4",
                 "magenta","cyan3", "plum4","darkorange1", "pink1", "red",
                 "lightgoldenrod", "violetred3", "brown4", "tan1","yellow", "gray50",
                 "#967272","deepskyblue4","salmon", "cyan","palegreen", "green3")
    if(is.null(values[["data"]])){
      data <- datarcheo
    }else{
      data <- values[["data"]]
    }
    if(!is.null(values[["index"]])){
      type <- values[["historique"]][,values[["index"]]+1]
    }else {
      type <- rep(1,240)
    }
    Size <- values[["Size"]]
    if(!is.null(Size)){
      sign = FALSE
      if(min(Size)<0){sign = TRUE}
      C <- values[["C"]]
      data[,3] <- as.numeric(Size)

      if(multi == "Univariate"){
        Size <- 10*(0.1+(abs(data[,3])/max(abs(data[,3]))))+5
      }

      pal <- colorFactor(couleur,as.factor(type))
      col <- pal(as.factor(type))
      col2 <- data.frame(col,type)
      col2 <- col2[order(col2[,2],decreasing=F), ]
      values[["col"]] <- as.character(col2[,1])
      leafletProxy("mymap", data = data) %>%
        clearMarkers() %>% clearControls
        if(input$Plabel){
          leafletProxy("mymap", data = data) %>% addCircleMarkers(lng = data[,1], lat =data[,2], popup = ~htmlEscape(values[["labelall"]]), radius = Size, fillColor =  col,
                                      fill = TRUE, fillOpacity = 0.8, opacity=1, stroke = TRUE,weight = 1.5, color="black", dashArray=C
          )
        }else{
          leafletProxy("mymap", data = data) %>% addCircleMarkers(lng = data[,1], lat = data[,2], radius = Size, fillColor = col,
                                      fill = TRUE, fillOpacity = 0.8, opacity=1, stroke = TRUE,weight = 1.5, color="black", dashArray=C
          )
        }
        leafletProxy("mymap", data = data) %>%
        fitBounds(lat1 = max(data[,2]), lng1 = max(data[,1]), lat2 = min(data[,2]), lng2 = min(data[,1])
        ) %>%
        addLegend("bottomleft", pal = pal, values = as.factor(type),
                  title = "Cluster",
                  opacity = 1
        )

      if(multi == "Univariate"){
        leafletProxy("mymap", data = data) %>%
          addLegendCustom(colors = c("black", "black", "black"),
                          labels = c(round(min(abs(data[,3])),2),
                                     round(min(abs(data[,3])) + (max(abs(data[,3]))-min(abs(data[,3])))/2  ,2),
                                     round(max(abs(data[,3])),2)),
                          sizes = c(10*(0.1+( min(abs(data[,3])) /max(abs(data[,3]))))+5,
                                    10*(0.1+((min(abs(data[,3])) + (max(abs(data[,3]))-min(abs(data[,3])))/2) /max(abs(data[,3]))))+5,
                                    10*(0.1+( max(abs(data[,3])) /max(abs(data[,3]))))+5),
                          title = "Size"
          )
      }
    }
    if(!is.null(Size)){
      if(sign){
        leafletProxy("mymap", data = data) %>%
          addLegendCustom(colors = c("white", "white"), labels = c("line: positive","dotted-line: negative"),
                          sizes = 1, title="Sign"
          )
      }
    }
    if(!is.null(Size)){
      if(lonlat){
        esri <- structure(c("Esri.WorldStreetMap", "Esri.DeLorme", "Esri.WorldTopoMap",
                            "Esri.WorldImagery", "Esri.WorldTerrain", "Esri.WorldShadedRelief",
                            "Esri.WorldPhysical", "Esri.OceanBasemap", "Esri.NatGeoWorldMap",
                            "Esri.WorldGrayCanvas"),
                          .Names = c("Esri.WorldStreetMap", "Esri.DeLorme", "Esri.WorldTopoMap",
                                     "Esri.WorldImagery", "Esri.WorldTerrain", "Esri.WorldShadedRelief",
                                     "Esri.WorldPhysical", "Esri.OceanBasemap", "Esri.NatGeoWorldMap",
                                     "Esri.WorldGrayCanvas"))
        m = leaflet()
        for (provider in esri) {
          m <- m %>% addProviderTiles(provider, group = provider)
        }
        m = m %>%
          addLayersControl(baseGroups = names(esri),
                           options = layersControlOptions(collapsed = TRUE),
                           position = "topright"
          ) %>%
          addMiniMap(tiles = esri[[1]], toggleDisplay = TRUE,minimized = TRUE,
                     position = "bottomright"
          ) %>%
          htmlwidgets::onRender("
                                function(el, x) {
                                var myMap = this;
                                myMap.on('baselayerchange',
                                function (e) {
                                myMap.minimap.changeLayer(L.tileLayer.provider(e.name));
                                })
                                }")
      }else{
        m <- leaflet()
      }

      m <- m %>%
      addCircleMarkers(lng = data[,1], lat = data[,2], radius = Size, fillColor =  col,
                       fill = T, fillOpacity = 0.8, opacity=1, stroke = T,weight = 1.5, color="black", dashArray=C
      ) %>%
      fitBounds(lat1 = max(data[,2]), lng1 = max(data[,1]), lat2 = min(data[,2]), lng2 = min(data[,1])
      ) %>%
      addLegend("bottomleft", pal = pal, values = as.factor(type),
                title = "Cluster",
                opacity = 1
      )
      if(multi == "Univariate"){
        m <- m %>%
          addLegendCustom(colors = c("black", "black", "black"),
                          labels = c(round(min(abs(data[,3])),2),
                                     round(min(abs(data[,3])) + (max(abs(data[,3]))-min(abs(data[,3])))/2  ,2),
                                     round(max(abs(data[,3])),2)),
                          sizes = c(10*(0.1+( min(abs(data[,3])) /max(abs(data[,3]))))+5,
                                    10*(0.1+((min(abs(data[,3])) + (max(abs(data[,3]))-min(abs(data[,3])))/2) /max(abs(data[,3]))))+5,
                                    10*(0.1+( max(abs(data[,3])) /max(abs(data[,3]))))+5),
                          title = "Size"
          )
      }
      if(sign){
        m <- m %>%
          addLegendCustom(colors = c("white", "white"), labels = c("line: positive","dotted-line: negative"),
                          sizes = 1, title="Sign"
          )
      }
      dlmap$map <- m
    }

  })

  ##silhouette plot
  observe({
    if(!is.null(values[["index"]])){
      type <- values[["historique"]][,values[["index"]]+1]
      DiMatrix = values[["DiMatrix"]]
      if(length(as.numeric(type)) != dim(DiMatrix)[1]){
        output$silPlot <- renderPlot({},height = 1)
      }else{
        if(length(unique(type)) == 1){
          showNotification("You have to create more than one patch. The \"Silhouette plot\" can't be calculate for one cluster.", type="warning", duration=10)
          output$silPlot <- renderPlot({},height = 1)
        }else{
          output$silPlot <- renderPlot({
            cluster = type
            DiMatrix = values[["DiMatrix"]]
            sil <- silhouette(cluster,DiMatrix)
            plot(sil,col=values[["col"]])
            values[["silhouetteplot"]] <- sil
          },height = function(){2.5*length(type)})
        }
      }
    }else{
      output$silPlot <- renderPlot({},height = 1)
    }
  })



  #EXPORT BUTTON
  output$all.pdf <- downloadHandler(
    filename = "mapclust.pdf",
    content = function(file) {
      cairo_pdf(file,onefile=T)
      print(values[["dendplot"]])
      print(values[["avesilplot"]])
      print(values[["WSSplot"]])
      if(!is.null(values[["silhouetteplot"]])){
        plot(values[["silhouetteplot"]],col=values[["col"]])
      }
      dev.off()
    }
  )

  }
