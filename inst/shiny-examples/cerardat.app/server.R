load(file = "./data/datacerardat.RData")

press <- function(linear.model) {
  # calculate the predictive residuals
  pr <- residuals(linear.model)/(1-lm.influence(linear.model)$hat)
  # calculate the PRESS
  PRESS <- mean(pr^2)
  return(PRESS)
}

log10Tck <- function(side, type){
  lim <- switch(side,
                x = par('usr')[1:2],
                y = par('usr')[3:4],
                stop("side argument must be 'x' or 'y'"))
  at <- floor(lim[1]) : ceiling(lim[2])
  return(switch(type,
                minor = outer(1:9, 10^(min(at):max(at))),
                major = 10^at,
                stop("type argument must be 'major' or 'minor'")
  ))
}

CA_var_data <- function(res, xax = 1, yax = 2, var_sup = TRUE, var_hide = "None",
                        var_lab_min_contrib = 0) {
  tmp_x <- res$vars %>%
    filter(Axis == xax) %>%
    select("Level", "Position", "Type", "Class", "Coord", "Contrib", "Cos2", "Count")
  tmp_y <- res$vars %>%
    filter(Axis == yax) %>%
    select("Level", "Position", "Type", "Class", "Coord", "Contrib", "Cos2", "Count")
  if (!var_sup) {
    tmp_x <- tmp_x %>% filter(Type == 'Active')
    tmp_y <- tmp_y %>% filter(Type == 'Active')
  }
  if (var_hide != "None") {
    tmp_x <- tmp_x %>% filter(Position != var_hide)
    tmp_y <- tmp_y %>% filter(Position != var_hide)
  }
  tmp <- tmp_x %>%
    left_join(tmp_y, by = c("Level", "Position", "Type", "Class", "Count")) %>%
    mutate(Contrib = Contrib.x + Contrib.y,
           Cos2 = Cos2.x + Cos2.y,
           tooltip = paste(paste0("<strong>", Level, "</strong><br />"),
                           paste0("<strong>",
                                  gettext("Position", domain = "R-explor"),
                                  ":</strong> ", Position, "<br />"),
                           paste0("<strong>Axis ",xax," :</strong> ", Coord.x, "<br />"),
                           paste0("<strong>Axis ", yax," :</strong> ", Coord.y, "<br />"),
                           ifelse(is.na(Cos2), "",
                                  paste0("<strong>",
                                         gettext("Squared cosinus", domain = "R-explor"),
                                         ":</strong> ", Cos2, "<br />")),
                           ifelse(is.na(Contrib), "",
                                  paste0("<strong>",
                                         gettext("Contribution:", domain = "R-explor"),
                                         "</strong> ", Contrib, "<br />")),
                           ifelse(is.na(Count), "",
                                  paste0("<strong>",
                                         gettext("Count:", domain = "R-explor"),
                                         "</strong> ", Count))),
           Lab = ifelse(Contrib >= as.numeric(var_lab_min_contrib) |
                          (is.na(Contrib) & as.numeric(var_lab_min_contrib) == 0), Level, ""))
  data.frame(tmp)
}

server <- function(input, output, session) {
  #------------------------------------------------------#
  #  modalDialog wait message
  #------------------------------------------------------#
  waitModal <- function() {
    modalDialog(title="Initialisation", size = c("s"),
                div(style="height:25px;"),
                h4("Calculation in progress"),
                span('This may take a while...'),
                hr(style="border-color: #222222;"),
                footer = NULL
    )
  }
  #------------------------------------------------------#
  #  modalDialog import
  #------------------------------------------------------#
  import1Modal <- function(failed = FALSE) {
    modalDialog(title="Reference dataset:",size=c("l"),
                div(style="height:15px;"),
                h3("Reference dataset",style="text-align:center;"),hr(style="border-color: #222222;"),
                p("reference data: count table plus a first column giving dates if known."),
                h4("Upload file"),
                # Input: Select a file ----
                fileInput("fileD1", span(icon("import", lib = "glyphicon"),"Choose CSV File"),
                          multiple = FALSE,
                          accept = c("text/csv",
                                     "text/comma-separated-values,text/plain",
                                     ".csv")
                ),
                actionButton(inputId = "toggle", span(icon("gears",lib = "font-awesome"),"Settings")),
                fluidRow(shinydashboard::box(id="myBox",
                                             fluidRow(
                                               column(width = 6,
                                                      # Input: Checkbox if file has header ----
                                                      span(strong("Header")),
                                                      checkboxInput("headerD1", "Header (colnames)", TRUE),
                                                      checkboxInput("rownamesD1", "rownames", TRUE)
                                               ),
                                               column(width = 6,
                                                      # Input: Select quotes ----input$dec
                                                      radioButtons("quoteD1", "Quote",
                                                                   choices = c(None = "",
                                                                               "Double Quote" = '"',
                                                                               "Single Quote" = "'"),
                                                                   selected = '"')
                                               )
                                             ),
                                             fluidRow(
                                               column(width = 6,
                                                      # Input: Select separator ----
                                                      radioButtons("sepD1", "Separator",
                                                                   choices = c(Comma = ",",
                                                                               Semicolon = ";",
                                                                               Tab = "\t",
                                                                               Space = " "),
                                                                   selected = ";")
                                               ),
                                               column(width = 6,
                                                      radioButtons("decD1", "Decimal",
                                                                   choices = c("Dot" = '.', "Comma" = ","),
                                                                   selected = '.')
                                               )
                                             )
                )),
                tableOutput("contentsD1"),
                if (failed)
                  div(tags$b("Your data must be match with the example above ! Change parameter or file.",
                             style = "color: red;")),
                footer = tagList(
                  modalButton("Cancel"),
                  actionButton("import2", "Next")
                )
    )
  }

  import2Modal <- function(failed = FALSE) {
    modalDialog(title="Supplementary dataset:",size=c("l"),
                div(style="height:15px;"),
                h3("Supplementary dataset",style="text-align:center;"),hr(style="border-color: #222222;"),
                p("Supplementary data: count table (same number of columns) plus a column for dates if known."),
                h4("Upload file"),
                # Input: Select a file ----
                fileInput("fileD2", span(icon("import", lib = "glyphicon"),"Choose CSV File"),
                          multiple = FALSE,
                          accept = c("text/csv",
                                     "text/comma-separated-values,text/plain",
                                     ".csv")
                ),
                actionButton(inputId = "toggle", span(icon("gears",lib = "font-awesome"),"Settings")),
                fluidRow(shinydashboard::box(id="myBox",
                                             fluidRow(
                                               column(width = 6,
                                                      # Input: Checkbox if file has header ----
                                                      span(strong("Header")),
                                                      checkboxInput("headerD2", "Header (colnames)", TRUE),
                                                      checkboxInput("rownamesD2", "rownames", TRUE)
                                               ),
                                               column(width = 6,
                                                      # Input: Select quotes ----input$dec
                                                      radioButtons("quoteD2", "Quote",
                                                                   choices = c(None = "",
                                                                               "Double Quote" = '"',
                                                                               "Single Quote" = "'"),
                                                                   selected = '"')
                                               )
                                             ),
                                             fluidRow(
                                               column(width = 6,
                                                      # Input: Select separator ----
                                                      radioButtons("sepD2", "Separator",
                                                                   choices = c(Comma = ",",
                                                                               Semicolon = ";",
                                                                               Tab = "\t",
                                                                               Space = " "),
                                                                   selected = ";")
                                               ),
                                               column(width = 6,
                                                      radioButtons("decD2", "Decimal",
                                                                   choices = c("Dot" = '.', "Comma" = ","),
                                                                   selected = '.')
                                               )
                                             )
                )),
                tableOutput("contentsD2"),
                if (failed)
                  div(tags$b("Your data must be match with the example above ! Change parameter or file.",
                             style = "color: red;")),
                footer = tagList(
                  actionButton("import1", "Back"),
                  actionButton("Finishimport", "Finish")
                )
    )
  }

  values <- reactiveValues()

  #raw data
  values[["date"]] = datacerardat$date
  values[["row.sup"]] = datacerardat$row.sup
  values[["df"]] = datacerardat$df
  values[["row.ref"]] = NULL
  values[["cerar_res"]] = NULL
  values[["D1"]] = NULL
  values[["D2"]] = NULL


  shinyjs::disable(selector='a[data-value="Results"]')
  shinyjs::disable(selector='a[data-value="Visualization"]')

  observeEvent(input$data,{
    showModal(import1Modal())
    shinyjs::hide(id = "myBox")
  })

  observeEvent(input$import1,{
    showModal(import1Modal())
    shinyjs::hide(id = "myBox")
  })

  observeEvent(input$import2,{
    showModal(import2Modal())
    shinyjs::hide(id = "myBox")
  })

  observeEvent(input$toggle, {
    if(input$toggle %% 2 == 0){
      shinyjs::hide(id = "myBox")
    }else{
      shinyjs::show(id = "myBox")
    }
  })

  observeEvent(input$Finishimport,{
    values[["date"]] = c(values[["D1"]][,1],values[["D2"]][,1])
    values[["row.sup"]] = (length(values[["D1"]][,1])+1):(length(values[["D1"]][,1])+length(values[["D2"]][,1]))
    values[["df"]] = rbind(values[["D1"]][,-c(1)],values[["D2"]][,-c(1)])
    values[["row.ref"]] = NULL
    values[["cerar_res"]] = NULL
    shinyjs::disable(selector='a[data-value="Results"]')
    shinyjs::disable(selector='a[data-value="Visualization"]')

    removeModal()
  })

  output$contentsD1 <- renderTable({

    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.

    req(input$fileD1)
    inFile <- input$fileD1

    if(is.null(inFile)){
      return(NULL)
    }

    tryCatch({
      if(input$rownamesD1){
        df <- read.table(input$fileD1$datapath,
                         header = input$headerD1,
                         sep = input$sepD1,
                         dec = input$decD1,
                         quote = input$quoteD1,
                         row.names = 1,
                         check.names = TRUE,
                         stringsAsFactors=FALSE)
      }else{
        df <- read.table(input$fileD1$datapath,
                         header = input$headerD1,
                         sep = input$sepD1,
                         dec = input$decD1,
                         quote = input$quoteD1,
                         check.names = TRUE,
                         stringsAsFactors=FALSE)
      }

    },
    error = function(e) {
      # return a safeError if a parsing error occurs
      stop(safeError(e))
    })
    values[["D1"]] <- df
    return(head(df))
  },rownames=TRUE)

  output$contentsD2 <- renderTable({

    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.

    req(input$fileD2)
    inFile <- input$fileD2

    if(is.null(inFile))
      return(NULL)

    tryCatch({
      if(input$rownamesD2){
        df <- read.table(input$fileD2$datapath,
                         header = input$headerD2,
                         sep = input$sepD2,
                         dec = input$decD2,
                         quote = input$quoteD2,
                         row.names = 1,
                         check.names = TRUE,
                         stringsAsFactors=FALSE)
      }else{
        df <- read.table(input$fileD2$datapath,
                         header = input$headerD2,
                         sep = input$sepD2,
                         dec = input$decD2,
                         quote = input$quoteD2,
                         check.names = TRUE,
                         stringsAsFactors=FALSE)
      }
    },
    error = function(e) {
      # return a safeError if a parsing error occurs
      stop(safeError(e))
    })
    values[["D2"]] <- df
    return(head(df))
  },rownames=TRUE)

  output$eigenvalue <- renderPlotly({
    date = values[["date"]]
    df = values[["df"]]
    row.sup = values[["row.sup"]]
    row.ref = which(!(1:length(df[,1]) %in% row.sup))

    #nettoyage des GT ie GT<5
    GT_rm_ref = which(colSums(df[row.ref,])<5)
    #?????????????????
    #GT_rm_sup = which(colSums(df[-c(GT_rm_ref),row.sup])<5)

    #data reference
    if(sum(GT_rm_ref)==0){
      data_ref = df[row.ref,]
    }else{
      data_ref = df[row.ref,-c(GT_rm_ref)]
    }

    DOF = min(ncol(data_ref)-1,nrow(data_ref)-1)
    ca.NMI.init <- ade4::dudi.coa(data_ref,scannf=FALSE,nf=DOF)

    inertia = ade4::inertia.dudi(ca.NMI.init, row.inertia=T, col.inertia=T)
    data_df = data.frame(
      Dim = 1:DOF,
      Inertia = inertia$tot.inertia[,1] / sum(inertia$tot.inertia[,1]),
      Cumul = inertia$tot.inertia[,3]
    )
    p = ggplot(data_df, aes(x=Dim, y=Inertia,text = paste("Percentage:",arrondi(Inertia,4)*100,"%<br>Cumul:",arrondi(Cumul,2),"%"))) +
      geom_bar(stat="identity",fill="steelblue")+
      theme_minimal()
    p
  })

  output$RMSE <- renderPlot({
    df = values[["df"]]
    date = values[["date"]]
    row.sup = values[["row.sup"]]
    row.ref = which(!(1:length(df[,1]) %in% row.sup))

    #vector for ref data (every not in row.sup)

    #nettoyage des GT ie GT<5
    GT_rm_ref = which(colSums(df[row.ref,])<5)
    #?????????????????
    #GT_rm_sup = which(colSums(df[-c(GT_rm_ref),row.sup])<5)

    #data reference
    if(sum(GT_rm_ref)==0){
      data_ref = df[row.ref,]
    }else{
      data_ref = df[row.ref,-c(GT_rm_ref)]
    }

    #DOF Degrees of freedom
    DOF = min(ncol(data_ref)-1,nrow(data_ref)-1)

    #Correspondance analysis
    ca.NMI <- ade4::dudi.coa(data_ref,scannf=FALSE,nf=DOF)

    #df with coord ens (li) de reference + date monnaie
    DATA_REF = cbind(ca.NMI$li, date=date[row.ref])
    max = min(length(DATA_REF[!is.na(DATA_REF$date),1])-2,DOF)


    #specify the cross-validation method
    formula <- "date ~"
    #ctrl <- caret::trainControl(method = "LOOCV")
    MSE <- c()
    PRESS <- c()
    R_sq <- c()

    for(i in 1:max){
      if(i == 1){
        formula <- paste0(formula,' Axis',i)
      }else{
        formula <- paste0(formula,' + Axis',i)
      }

      lm = lm(as.formula(formula),data = DATA_REF)
      MSE <- c(MSE,mean(lm$residuals^2))
      PRESS <- c(PRESS,press(lm))

    }
    PRESS_max=max(PRESS[is.finite(PRESS)])

    if(input$quid == "both"){
      plot(x=1:max,y=PRESS,log="y",yaxt='n',type="o",pch=1,ylab="PRediction Error Sum Of Squares (PRESS)",
           xlab="Number of component in lm()",ylim = c(min(MSE),PRESS_max))
      points(x = which(PRESS == min(PRESS)), y=min(PRESS),pch=16,cex=1.1)
      axis(2, at=log10Tck('y','major'), tck= 0.02, labels=format(log10Tck('y','major'), scientific = T))
      axis(2, at=log10Tck('y','minor'), tck= 0.01, labels=NA)
      abline(v=which(PRESS == min(PRESS)),lty=3,col="black")
      points(x=1:max,y=MSE,type="o",pch=1,ylab="",xlab="",col="#cc525b")
      mtext("Mean Squared Error (MSE)",side=4,col="#cc525b",line=.5)
    }
    if(input$quid == "PRESS"){
      plot(x=1:max,y=PRESS,log="y",yaxt='n',type="o",pch=1,ylab="PRediction Error Sum Of Squares (PRESS)",
           xlab="Number of component in lm()")
      points(x = which(PRESS == min(PRESS)), y=min(PRESS),pch=16,cex=1.1)
      axis(2, at=log10Tck('y','major'), tck= 0.02, labels=format(log10Tck('y','major'), scientific = T))
      axis(2, at=log10Tck('y','minor'), tck= 0.01, labels=NA)
      abline(v=which(PRESS == min(PRESS)),lty=3,col="black")
    }
    if(input$quid == "MSE"){
      plot(x=1:max,y=MSE,type="o",pch=1,ylab="Mean Squared Error (MSE)",col="#cc525b",
           xlab="Number of component in lm()",ylim = c(0,max(MSE)))
    }



    shinyjs::enable(selector='a[data-value="Results"]')
    shinyjs::enable(selector='a[data-value="Visualization"]')

  })

  output$R_sq <- renderPlot({
    df = values[["df"]]
    date = values[["date"]]
    row.sup = values[["row.sup"]]
    row.ref = which(!(1:length(df[,1]) %in% row.sup))

    #vector for ref data (every not in row.sup)

    #nettoyage des GT ie GT<5
    GT_rm_ref = which(colSums(df[row.ref,])<5)
    #?????????????????
    #GT_rm_sup = which(colSums(df[-c(GT_rm_ref),row.sup])<5)

    #data reference
    if(sum(GT_rm_ref)==0){
      data_ref = df[row.ref,]
    }else{
      data_ref = df[row.ref,-c(GT_rm_ref)]
    }

    #DOF Degrees of freedom
    DOF = min(ncol(data_ref)-1,nrow(data_ref)-1)

    #Correspondance analysis
    ca.NMI <- ade4::dudi.coa(data_ref,scannf=FALSE,nf=DOF)

    #df with coord ens (col) de reference + date monnaie
    DATA_REF = cbind(ca.NMI$li, date=date[row.ref])
    max = min(length(DATA_REF[!is.na(DATA_REF$date),1])-2,DOF)


    #specify the cross-validation method
    formula <- "date ~"
    #ctrl <- caret::trainControl(method = "LOOCV")
    MSE <- c()
    MSE_sd <- c()
    PRESS <- c()
    R_sq <- c()

    for(i in 1:max){
      if(i == 1)
        formula <- paste0(formula,' Axis',i)
      else
        formula <- paste0(formula,' + Axis',i)

      #model <- caret::train(as.formula(formula), data = DATA_REF, method = "lm", trControl = ctrl, na.action=na.omit)
      #RMSE <- c(RMSE,model$results$RMSE)
      #R_sq <- c(R_sq,model$results$Rsquared)

      lm = lm(as.formula(formula),data = DATA_REF)
      R_sq <- c(R_sq,summary(lm)$adj.r.squared)

    }
    min_R = min(R_sq)
    if(min(R_sq) > 0)min_R= 0
    if(min(R_sq) > 0.25)min_R= 0.25
    if(min(R_sq) > 0.5)min_R= 0.5
    if(min(R_sq) > 0.75)min_R= 0.75

    plot(x=1:max,y=R_sq,type="o",pch=1,main="Adjusted R²",ylab="Adj.R²",
           xlab="Number of component in lm()",ylim = c(min_R,1))

  })


  #output$checkbox <- renderText({
  #  row.ref = which(!(1:length(values[["df"]][1,]) %in% values[["row.sup"]]))
  #  date = values[["date"]]
  #
  #  updateCheckboxGroupInput(session, "row.refbox",
  #                           choiceNames = as.list(colnames(values[["df"]])),
  #                           choiceValues = as.list(1:ncol(values[["df"]])),
  #                           selected = row.ref
  #  )
  #  ""
  #})

  output$cerar <- renderText({
    date = values[["date"]]
    df = values[["df"]]
    row.sup = values[["row.sup"]]
    row.ref = which(!(1:length(df[,1]) %in% row.sup))

    #todo: while convergence

    #check rowsum ens
    if(sum(rowSums(df)<5)!=0){
      warning(paste0("The sums of rows ",capture.output(cat(row.names(df)[which(rowSums(df)<5)]))," are less than 5. They were suppressed from the analysis."))
      #retire index de la col dans row.sup
      tmp.row.status = data.frame(
        index = c(row.ref,row.sup),
        status = c(rep("ref",length(row.ref)),rep("sup",length(row.sup)))
      )

      tmp.row.status = tmp.row.status[order(tmp.row.status$index),]
      #new row.sup
      tmp.row.status = tmp.row.status[!rowSums(df)<5,]
      #retire la col de date
      date = date[!rowSums(df)<5]
      #retire la col de df
      df = df[!rowSums(df)<5,]

      row.ref = which(tmp.row.status$status == "ref")
      row.sup = which(tmp.row.status$status == "sup")

      values[["df"]] = df
      values[["date"]] = date
      values[["row.sup"]] = row.sup
    }

    values[["cerar_res"]] = cerardat(df, row.sup, date, nf = input$k, graph=F, confidence = input$confidence)
    ""
  })

  output$lm <- renderPrint({
    summary(values[["cerar_res"]]$lm)
  })

  output$SW <- renderPrint({
    print(values[["cerar_res"]]$Shapiro_Wilks)
  })

  output$DW <- renderPrint({
    print(values[["cerar_res"]]$Durbin_Watson)
  })

  output$BP <- renderPrint({
    print(values[["cerar_res"]]$Breusch_Pagan)
  })

  output$sw <- renderText({
    if(values[["cerar_res"]]$Shapiro_Wilks$p.value < 0.05){
      "Warning: The Shapiro-Wilks test indicates a problem with the normality of the residuals."
    }
  })

  output$dw <- renderText({
    if(values[["cerar_res"]]$Durbin_Watson$p.value < 0.05){
      "Warning: The Durbin-Watson test indicates a first order autocorrelation problem."
    }
  })

  output$bp <- renderText({
    if(values[["cerar_res"]]$Breusch_Pagan$p.value < 0.05){
      "Warning: The Breusch-Pagan test indicates a heteroskedasticity problem."
    }
  })

  output$lmplot <- renderUI({
    if(is.null(values[["cerar_res"]])){
      cerar_res = cerardat(df, row.sup, date, nf = input$k, graph=F, confidence = input$confidence)
    }else{
      cerar_res = values[["cerar_res"]]
    }

    plot_output_list <- lapply(c(1:4), function(i) {
      name <- paste("lm_plot", i, sep="")
      output[[name]] <- renderPlot({
        plot(cerar_res$lm,which = i,ask=FALSE)
      },width = 450,height = 450)
    })
    shinycssloaders::withSpinner(do.call(tagList,plot_output_list),color='#18bc9c')

  })


  output$evaldateref <- renderPlotly({
    date = values[["date"]]
    df = values[["df"]]

    row.sup = values[["row.sup"]]
    row.ref = which(!(1:length(df[,1]) %in% row.sup))

    if(is.null(values[["cerar_res"]])){
      cerar_res = cerardat(df, row.sup, date, nf = input$k, graph=F, confidence = input$confidence)
    }else{
      cerar_res = values[["cerar_res"]]
    }


    tmp = cerar_res$prediction[row.ref,]
    tmp = tmp[!is.na(tmp$date),]


    df = data.frame(
      names = rownames(tmp),
      date = tmp$date,
      lwr = tmp$lower_Ev,
      upr = tmp$upper_Ev,
      fit = tmp$Fit_dateEv
    )

    #generation du graph
    graph = suppressWarnings(ggplot(df) +
      geom_point(aes(x = names, y = date, text=paste("lower:",lwr,"<br>upper:",upr)),colour="red",shape=3, size=2) +
      geom_errorbar( aes(x=names, ymin=lwr, ymax=upr), width=.15, colour="#3f0b18", alpha=0.9, linewidth=1)+
      ylab("Date (year)") + xlab("Site") + ggtitle("") +
      theme_classic() +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
      theme(
        panel.grid.major.x = element_line(color = rgb(0.5,0.5,0.5,0.3),
                                          linewidth = .1,
                                          linetype = 2),
        panel.grid.major.y = element_line(color = rgb(0.5,0.5,0.5,0.3),
                                          linewidth = .1,
                                          linetype = 2),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
      ))
    graph

  })

  output$evaldatesup <- renderPlotly({
    date = values[["date"]]
    df = values[["df"]]
    row.sup = values[["row.sup"]]

    if(is.null(values[["cerar_res"]])){
      cerar_res = cerardat(df, row.sup, date, nf = input$k, graph=F, confidence = input$confidence)
    }else{
      cerar_res = values[["cerar_res"]]
    }


    tmp = cerar_res$prediction[row.sup,]
    tmp = tmp[!is.na(tmp$date),]

    if(dim(tmp)[1] == 0){
      graph = ggplot() +
        geom_blank() +
        theme_classic()
      graph
    }else{
      df = data.frame(
        names = rownames(tmp),
        date = tmp$date,
        lwr = tmp$lower_Ev,
        upr = tmp$upper_Ev,
        fit = tmp$Fit_dateEv
      )

      #generation du graph
      graph = suppressWarnings(ggplot(df) +
                                 geom_point(aes(x = names, y = date, text=paste("lower:",lwr,"<br>upper:",upr)),colour="red",shape=3, size=2) +
                                 geom_errorbar( aes(x=names, ymin=lwr, ymax=upr), width=.15, colour="#3f0b18", alpha=0.9, linewidth=1)+
                                 ylab("Date (year)") + xlab("Site") + ggtitle("") +
                                 theme_classic() +
                                 scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
                                 theme(
                                   panel.grid.major.x = element_line(color = rgb(0.5,0.5,0.5,0.3),
                                                                     linewidth = .1,
                                                                     linetype = 2),
                                   panel.grid.major.y = element_line(color = rgb(0.5,0.5,0.5,0.3),
                                                                     linewidth = .1,
                                                                     linetype = 2),
                                   axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
                                 ))
      graph
    }



  })

  observeEvent(input$axe1,{
    values[["axe1"]] <- input$axe1
  })

  observeEvent(input$axe2,{
    values[["axe2"]] <- input$axe2
  })

  observeEvent(input$confidence,{
    df = values[["df"]]
    date = values[["date"]]
    row.sup = values[["row.sup"]]

    values[["cerar_res"]] = cerardat(df, row.sup, date, nf = input$k, graph=F, confidence = input$confidence)
  })



  output$CA <- renderScatterD3({
    df = values[["df"]]
    row.sup = values[["row.sup"]]
    row.ref = which(!(1:length(df[,1]) %in% row.sup))

    #nettoyage des GT ie GT<5
    GT_rm_ref = which(colSums(df[row.ref,])<5)

    #?????????????????
    #GT_rm_sup = which(colSums(df[-c(GT_rm_ref),row.sup])<5)

    #data reference
    if(sum(GT_rm_ref)==0){
      data_ref = df[row.ref,]
    }else{
      data_ref = df[row.ref,-c(GT_rm_ref)]
    }
    DOF = min(ncol(data_ref)-1,nrow(data_ref)-1)

    updateNumericInput(session, "axe1", max = DOF)
    updateNumericInput(session, "axe2", max = DOF)

    res <- explor::prepare_results(dudi.coa(data_ref,nf=DOF,scannf = FALSE))

    xax = values[["axe1"]]
    yax = values[["axe2"]]
    var_sup = TRUE
    var_hide = "None"
    var_lab_min_contrib = 0
    point_size = 64
    col_var = "Position"
    symbol_var = NULL
    size_var = NULL
    size_range = c(10,300)
    zoom_callback = NULL
    in_explor = FALSE

    ## Settings changed if not run in explor
    html_id <- if(in_explor) "explor_var" else  NULL
    dom_id_svg_export <- if(in_explor) "explor-var-svg-export" else NULL
    dom_id_lasso_toggle <- if(in_explor) "explor-var-lasso-toggle" else NULL
    lasso <- if(in_explor) TRUE else FALSE
    lasso_callback <- if(in_explor) explor_multi_lasso_callback() else NULL
    zoom_callback <- if(in_explor) explor_multi_zoom_callback(type = "var") else NULL

    var_data <- CA_var_data(res, xax, yax)

    scatterD3::scatterD3(
      x = var_data[, "Coord.x"],
      y = var_data[, "Coord.y"],
      xlab = names(res$axes)[res$axes == xax],
      ylab = names(res$axes)[res$axes == yax],
      lab = var_data[, "Lab"],
      point_size = point_size,
      point_opacity = 1,
      col_var = if (is.null(col_var)) NULL else var_data[,col_var],
      col_lab = col_var,
      symbol_var = if (is.null(symbol_var)) NULL else var_data[,symbol_var],
      symbol_lab = symbol_var,
      size_var = if (is.null(size_var)) NULL else var_data[,size_var],
      size_lab = size_var,
      size_range = if (is.null(size_var)) c(10,300) else c(30,400) * point_size / 32,
      tooltip_text = var_data[, "tooltip"],
      type_var = ifelse(var_data[, "Class"] == "Quantitative", "arrow", "point"),
      unit_circle = var_sup && "Quantitative" %in% var_data[,"Class"],
      key_var = paste(var_data[,"Position"], var_data[,"Level"], sep = "-"),
      fixed = TRUE,
      html_id = html_id,
      dom_id_svg_export = dom_id_svg_export,
      dom_id_lasso_toggle = dom_id_lasso_toggle,
      lasso = lasso,
      lasso_callback = lasso_callback,
      zoom_callback = zoom_callback
    )
  })


  output$cerardat_plot_ref <- renderUI({
    df = values[["df"]]

    row.sup = values[["row.sup"]]
    row.ref = which(!(1:length(df[,1]) %in% row.sup))


    progress <- Progress$new(session, min=1, max=length(row.ref))
    on.exit(progress$close())

    progress$set(message = 'Calculation in progress',
                 detail = 'This may take a while...')

    if(is.null(values[["cerar_res"]])){
      cerar_res = cerardat(df, row.sup, date, nf = input$k, graph=F, confidence = input$confidence)
    }else{
      cerar_res = values[["cerar_res"]]
    }

    GT_date_sd = arrondi(cbind(
      Pred=data.frame(cerar_res$predict_obj_col$fit)[,1],
      data.frame(E_T=cerar_res$predict_obj_col$se.fit)
    ),1)

    ecart = max(cerar_res$date_gt)-min(cerar_res$date_gt)
    xlim = c(min(cerar_res$date_gt-ecart*0.07),max(cerar_res$date_gt+ecart*0.07))

    tmp_ = c()
    for(i in 1:ncol(cerar_res$cont_gt) )
    {
      date_accumulation <- nor1mix::norMix(mu = GT_date_sd[,1], w = unlist(as.vector(cerar_res$cont_gt[i,])), sigma= GT_date_sd[,2])
      date_accumulation_density = nor1mix::dnorMixL(date_accumulation, xlim=xlim)

      tmp_ = c(tmp_,max(date_accumulation_density$y))
      progress$set(value = i)
    }
    ylim=c(0,max(tmp_))

    plot_output_list <- lapply(row.ref, function(i) {
      name <- paste("cera_plot", i, sep="")
      output[[name]] <- renderPlot({
        plot(cerar_res,which=i,xlim=xlim,ylim=ylim)
      },width = 450,height = 450)
    })
    shinycssloaders::withSpinner(do.call(tagList,plot_output_list),color='#18bc9c')

  })

  output$cerardat_plot_sup <- renderUI({
    df = values[["df"]]

    row.sup = values[["row.sup"]]
    row.ref = which(!(1:length(df[,1]) %in% row.sup))

    progress <- Progress$new(session, min=1, max=length(row.sup))
    on.exit(progress$close())

    progress$set(message = 'Calculation in progress',
                 detail = 'This may take a while...')

    if(is.null(values[["cerar_res"]])){
      cerar_res = cerardat(df, row.sup, date, nf = input$k, graph=F, confidence = input$confidence)
    }else{
      cerar_res = values[["cerar_res"]]
    }

    GT_date_sd = arrondi(cbind(
      Pred=data.frame(cerar_res$predict_obj_col$fit)[,1],
      data.frame(E_T=cerar_res$predict_obj_col$se.fit)
    ),1)

    ecart = max(cerar_res$date_gt)-min(cerar_res$date_gt)
    xlim = c(min(cerar_res$date_gt-ecart*0.07),max(cerar_res$date_gt+ecart*0.07))

    tmp_ = c()
    for(i in 1:ncol(cerar_res$cont_gt) )
    {
      date_accumulation <- nor1mix::norMix(mu = GT_date_sd[,1], w = unlist(as.vector(cerar_res$cont_gt[i,])), sigma= GT_date_sd[,2])
      date_accumulation_density = nor1mix::dnorMixL(date_accumulation, xlim=xlim)

      tmp_ = c(tmp_,max(date_accumulation_density$y))
      progress$set(value = i)
    }
    ylim=c(0,max(tmp_))

    plot_output_list <- lapply(row.sup, function(i) {
      name <- paste("cera_plot", i, sep="")
      output[[name]] <- renderPlot({
        plot(cerar_res,which=i,xlim=xlim,ylim=ylim)
      },width = 450,height = 450)
    })
    shinycssloaders::withSpinner(do.call(tagList,plot_output_list),color='#18bc9c')

  })


  #EXPORT BUTTON
  output$downloadData <- downloadHandler(
    filename = "output.pdf",
    content = function(file) {
      cairo_pdf(file,onefile=T)
      for (i in 1:ncol(values[["df"]])) {
        plot(values[["cerar_res"]],which=i)
      }
      dev.off()
    }
  )

  output$tablesite <- renderTable({
    if(is.null(values[["cerar_res"]])){
      df=data.frame()
    }else{
      df=values[["cerar_res"]]$prediction
    }
  },digits = 0, rownames = T, bordered = T, striped = T, hover = T)

  output$tableGT <- renderTable({
    if(is.null(values[["cerar_res"]])){
      df=data.frame()
    }else{
      df=values[["cerar_res"]]$date_gt
    }
  },digits = 0, rownames = T, bordered = T, striped = T, hover = T)


  output$formula_dateEv <- renderUI({
    withMathJax(paste0("Linear model: $$ \\text { dateEv }_i=\\beta_0 + \\sum_{k=1}^K \\beta_k\\left(F^k\\right)_i+\\varepsilon_i \\quad \\forall i=1, \\ldots, I $$ for K the number of component kept in the model and I the number of observations with a date."))
  })

  observeEvent(input$data,{
    shinyjs::disable("import2")
  })

  observeEvent(input$import2,{
    shinyjs::disable("Finishimport")
  })

  observeEvent(input$fileD1,{
    shinyjs::enable("import2")
  })

  observeEvent(input$fileD2,{
    shinyjs::enable("Finishimport")
  })




}
