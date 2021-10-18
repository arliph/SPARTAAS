server <- function(input, output, session) {
  #------------------------------------------------------#
  #  modalDialog import
  #------------------------------------------------------#
  import1Modal <- function(failed = FALSE) {
    modalDialog(title="Step 1:",size=c("l"),
                div(style="height:15px;"),
                h3("Step 1",style="text-align:center;"),hr(style="border-color: #222222;"),
                h4("Upload file contingency table"),
                # Input: Select a file ----
                fileInput("fileD1", span(icon("import", lib = "glyphicon"),"Choose CSV File"),
                          multiple = FALSE,
                          accept = c("text/csv",
                                     "text/comma-separated-values,text/plain",
                                     ".csv")
                ),
                actionButton(inputId = "toggle", span(icon("cogs",lib = "font-awesome"),"Settings")),
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
                  actionButton("Finishimport", "Finish")
                )
    )
  }
  #------------------------------------------------------#
  #  modalDialog sort periods
  #------------------------------------------------------#
  sortPeriod <- function(failed = FALSE) {
    modalDialog(title="Sort Periods:",size=c("m"),
                div(style="height:15px;"),
                h3("Sort periods",style="text-align:center;"),hr(style="border-color: #222222;"),
                h4("Sort each cluster from the oldest to the newest!"),
                helpText("Note: Drag and drop clusters to sort them.
                         Drag unused clusters into the trash to remove them from the seriograph.
                         You can insert 'Hiatus' by picking it up and dropping it off where you want."),
                orderInput('seq2',span(icon("sort-amount-up",lib = "font-awesome"),'Sort the cluster'),
                           items = labels(values[["data1"]])[[1]],
                           placeholder = 'There must be at least one item',
                           item_class = 'success'),
                hr(style="border-color: #b2b2b2;"),
                if (failed)
                  div(tags$b("You have to select the sequence!", style = "color: red;")),
                footer = tagList(
                  modalButton("Cancel"),
                  actionButton("tri", "OK")
                )
    )
  }

  values <- reactiveValues()
  values[["seriograph"]] <- NULL
  shinyjs::hide("toggle")
  #contingency table
  values[["data1"]] <- data.frame(
    Cat10 = c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
    Cat20 = c(4,8,0,0,0,0,0,0,0,6,0,0,0,0,0,0,0),
    Cat30 = c(18,24,986,254,55,181,43,140,154,177,66,1,24,15,0,31,37),
    Cat40.50 = c(17,121,874,248,88,413,91,212,272,507,187,40,332,174,17,288,224),
    Cat60 = c(0,0,1,0,0,4,4,3,0,3,0,0,0,0,0,0,0),
    Cat70 = c(3,1,69,54,10,72,7,33,74,36,16,4,40,5,0,17,13),
    Cat80 = c(4,0,10,0,12,38,2,11,38,26,25,1,18,4,0,25,7),
    Cat100.101 = c(23,4,26,51,31,111,36,47,123,231,106,21,128,77,10,151,114),
    Cat102 = c(0,1,2,2,4,4,13,14,6,6,0,0,12,5,1,17,64),
    Cat110.111.113 = c(0,0,22,1,17,21,12,20,30,82,15,22,94,78,18,108,8),
    Cat120.121 = c(0,0,0,0,0,0,0,0,0,0,66,0,58,9,0,116,184),
    Cat122 = c(0,0,0,0,0,0,0,0,0,0,14,0,34,5,0,134,281),
    row.names = c("AI01","AI02","AI03","AI04","AO03","AI05","AO01","AI07","AI08",
                  "AO02","AI06","AO04","APQR01","APQR02","AO05","APQR03","AI09")
  )

  observeEvent(input$width,{
    values[["width"]] <- input$width
  })
  observeEvent(input$height,{
    values[["height"]] <- input$height
  })

  output$EPPM <- renderPlot({
    cont <- values[["data1"]]
    show <- c("both","frequency","EPPM")[as.numeric(input$show)]
    seriograph(cont, permute=input$permute, col_weight = input$col_weight, show = show)
    values[["seriograph"]] <- recordPlot()
  })



  observeEvent(input$data,{
    showModal(import1Modal())
    shinyjs::hide(id = "myBox")
  })

  observeEvent(input$Finishimport,{
    values[["data1"]] <- values[["D1"]]
    removeModal()
  })

  observeEvent(input$sort,{
    showModal(sortPeriod())
  })

  output$contentsD1 <- renderTable({

    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.

    req(input$fileD1)
    inFile <- input$fileD1

    if(is.null(inFile))
      return(NULL)

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


  #EXPORT BUTTON
  output$seriograph.pdf <- downloadHandler(
    filename = "seriograph.pdf",
    content = function(file) {
      cairo_pdf(file,width=12)
      print(values[["seriograph"]])
      dev.off()
    }
  )

  output$seriograph.png <- downloadHandler(
    filename = "seriograph.png",
    content = function(file) {
      png(file,
          width = values[["width"]],
          height = values[["height"]]
      )
      print(values[["seriograph"]])
      dev.off()
    }
  )


  observeEvent(input$toggle, {
    if(input$toggle %% 2 == 0){
      shinyjs::hide(id = "myBox")
    }else{
      shinyjs::show(id = "myBox")
    }
  })

  observeEvent(input$tri,{
    order <- input$seq2_order

    new_order <- rep(NA,length(order[,2]))

    for(i in 1:length(order[,2])){
      new_order[i] <- which(labels(values[["data1"]])[[1]] == order[i,2])
    }
    values[["data1"]] <- values[["data1"]][new_order,]

    removeModal()
  })

}
