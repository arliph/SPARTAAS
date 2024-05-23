EPM <- function(x) {
  x <- as.matrix(x)

  x[is.na(x)] <- 0
  ni. <- rep(1,nrow(x))
  nj. <- margin.table(x, margin = 2)
  n. <- sum(nj.)

  TabInde <- (nj.%*%t(ni.)/n.)
  TabInde <- t(TabInde)
  x <- prop.table(x, margin = 1)
  TabEcart <- x - TabInde
  TabEcart <- as.data.frame(TabEcart)
  return(TabEcart)
}

plot_EPPM <- function(x,show,permute) {
  if(is.table(x)){
    x <- as.data.frame(x)
    x <- matrix(x[,3],
                ncol=sqrt(length(x[,2])),
                dimnames=list(unique(x[,1]),unique(x[,1])))
    x <- as.data.frame(x)
  }
  if(is.matrix(x)){x <- as.data.frame(x)}
  if(!is.data.frame(x)){stop("x is not a data frame.")}
  for (i in 1:length(x[1,])){
    if(!(sum(x[,i]) == round(sum(x[,i])))){stop("The data frame must contains integers.")}
  }
  .data <- c()
  Cluster <- c()
  labels <- factor(rownames(x),levels = rev(unique(rownames(x))))
  #permute type in x
  if(permute){

    if(show == "EPPM"){
      ecart <- as.data.frame(EPM(x))
      ecart[ecart < 0] <- 0
      i <- seq(nrow(ecart),1)
      barycenter <- colSums(ecart * i) / colSums(ecart)
      x <- x[,order(barycenter)]
    }else{
      i <- seq(nrow(x),1)
      barycenter <- colSums(x * i) / colSums(x)
      x <- x[,order(barycenter)]
    }

  }
  rowsum <- rowSums(x)
  #add Weight to x
  x <- cbind(x, Weight = rowSums(x) / sum(x) * rowSums(x))
  #data
  data <- as.data.frame(x)
  data <- mutate(data,ens = labels,rowsum = rowsum)
  data <- gather(data,key = "type", value = "frequency",-ens,-rowsum, factor_key = TRUE)
  data <- mutate(data,frequency = frequency / rowsum)
  #remove col rowsum
  data <- data[-2]
  data$frequency[is.na(data$frequency)] <- 0

  #EPPM
  ecart <- as.data.frame(EPM(x[-length(x[1,])]))
  save_ecart <- ecart
  ecart <- cbind(ecart,Weight=c(0))
  ecart[ecart < 0] <- 0
  ecart <- mutate(ecart,ens = labels)
  ecart <- gather(ecart,key = "type", value = "EPPM", -ens, factor_key = TRUE)
  ecart$EPPM[ecart$type == "Weight"] <- 0

  # Join data and EPPM
  data <- inner_join(data,ecart,by = c("ens", "type"))
  data <- mutate(data,frequency = frequency - EPPM)
  data <- gather(data,key = "legend", value = "frequency", -ens,
                 -type, factor_key = TRUE)
  #if insert hiatus in df replace NA by 0 (NA are generating by divising by 0 (rowsum of a hiatus is 0))
  data[is.na(data)] <- 0
  #center (div /2 et copy for one part - and the other +)
  data <- rbind(data,data)
  data <- mutate(data,frequency = frequency * c(rep(.5, nrow(data)/2), rep(-.5, nrow(data)/2)))

  #color for Weight
  tmp <- data$frequency[data$type == "Weight"]
  levels(data$legend) <- c(levels(data$legend),"1","2","3","4","5","6","7","8","9","10","Weight")
  data$legend[data$type == "Weight"] <- "Weight"

  breaks <- c("frequency","EPPM","Weight")

  #define the color
  freq <- c("frequency","frequencies","freq","fq","fr\u00E9quence","fr\u00E9quences","fr\u00E9q")
  EPPM <- c("EPPM","ecart","\u00E9cart","ecarts","\u00E9carts")
  if(show %in% freq){
    color <- c("frequency"="#bebebe","EPPM"="#bebebe","Weight"="#7f7f7f")
    #remove EPPM from legend
    breaks <- c("frequency","Weight")
    default_col <- "#bebebe"
  }else{
    if(show %in% EPPM){
      color <- c("frequency"=NA,"EPPM"="black","Weight"="#7f7f7f")
      #remove frequency from legend
      breaks <- c("EPPM","Weight")
      default_col <- NA
    }else{
      #default color grey: frequency, black: EPPM, my_palette: weight
      color <- c("frequency"="#bebebe","EPPM"="black","Weight"="#7f7f7f")
      default_col <- NA
    }
  }
  #ggplot
  p <- ggplot(data = data) +
    facet_grid( ~type, scales = "free", space = "free_x") +
    geom_col(aes(x = ens, y = frequency, fill = legend), width = 1) +
    scale_x_discrete(labels = rev(rownames(x))) +
    coord_flip() +
    scale_fill_manual(values = color,
                      aesthetics = "fill",
                      limits = breaks,
                      na.value = default_col
    ) +
    scale_y_continuous(breaks = c(-.1,.1),labels = c("","20% ")) +
    labs(
      title = "Seriograph",
      subtitle = "EPPM (Ecart Positif au Pourcentage Moyen) Positive Deviation at Average Percentage",
      caption = "Warnings: The frequencies are calculated independently for each element (row).
You can see the relative number of observations in the Weight column"
    ) +
    theme_grey(base_size = 15) +
    annotate("segment",x = 0, y = -.04, xend = 0, yend = .04,color = "white") +
    theme(legend.position = "bottom",
          panel.spacing = unit(.2, "lines"),
          strip.text.x = element_text(angle = 90),
          axis.text.x = element_text(hjust=1,margin = margin(t = -13)),
          axis.ticks.length = unit(.5, "cm"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.caption = element_text(face = "italic", hjust = 0)
    )
  plot(p)
}

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
    modalDialog(title="Sort clusters:",size=c("m"),
                div(style="height:15px;"),
                h3("Sort rows",style="text-align:center;"),hr(style="border-color: #222222;"),
                h4("Sort each row from the oldest to the newest!"),
                helpText("Note: Drag and drop rows to sort them.
                         Drag unused row into the trash to remove them from the seriograph.
                         You can insert 'Hiatus' by picking it up and dropping it off where you want."),
                orderInput('seq2',span(icon("sort-amount-up",lib = "font-awesome"),'Sort the row'),
                           items = labels(values[["data1"]])[[1]] [ values[["new_order"]] ] ,
                           placeholder = 'There must be at least one item',
                           connect = "rm",
                           item_class = 'success'),
                hr(style="border-color: #b2b2b2;"),
                column(4,
                       orderInput('hiatus', span(icon("plus",lib = "font-awesome"),'Add'), items = c("Hiatus"),
                                  connect = 'seq2', as_source = TRUE, item_class = 'warning', width = '70px')
                ),
                column(8,
                       orderInput('rm', span(icon("trash",lib = "font-awesome"),'Remove from seriograph'),
                                  items = values[["rm"]],
                                  connect = "seq2",
                                  placeholder = 'Drag the items here to delete them',
                                  item_class = 'danger')
                ),
                div(style="height:50px;"),
                if (failed)
                  div(tags$b("You have to select the sequence!", style = "color: red;")),
                footer = tagList(
                  modalButton("Cancel"),
                  actionButton("tri", "OK")
                )
    )
  }

  values <- reactiveValues()
  values[["hiatus"]] <- NULL
  values[["rm"]] <- NULL
  values[["seriograph"]] <- NULL
  shinyjs::hide("toggle")
  #contingency table
  df = data.frame(
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
  values[["data1"]] <- df
  values[["new_order"]] <- 1:17

  #url param
  o=observe({
    query <- parseQueryString(session$clientData$url_search)

    for (i in 1:(length(reactiveValuesToList(input)))) {
      nameval = names(reactiveValuesToList(input)[i])
      valuetoupdate = query[[nameval]]

      if (!is.null(query[[nameval]])) {
        if (is.na(as.numeric(valuetoupdate))) {
          #updateTextInput(session, nameval, value = valuetoupdate)
          if(nameval == 'data'){
            data1 = read.csv(valuetoupdate,header =TRUE, row.names=1, sep=",")
            if(ncol(data1)==0){
              data1 = read.csv(valuetoupdate,header =TRUE, row.names=1, sep=";")
            }
            values[["data1"]] = data1
            values[["new_order"]] <- 1:length(labels(data1)[[1]])
          }
        }
        else {
          #updateTextInput(session, nameval, value = as.numeric(valuetoupdate))
        }
      }

    }
    o$destroy()

  })

  observeEvent(input$width,{
    values[["width"]] <- input$width
  })
  observeEvent(input$height,{
    values[["height"]] <- input$height
  })



  output$EPPM <- renderPlot({
    cont <- values[["data1"]]
    show <- c("both","frequency","EPPM")[as.numeric(input$show)]

    hiatus <- values[["hiatus"]]
    new_order <- values[["new_order"]]

    #reorder and rm
    if(length(new_order) > 0){
      cont <- cont[new_order,]
    }

    #hiatus
    if(length(hiatus) > 0 && sum(is.na(hiatus)) == 0){
      #insert hiatus
      label <- "Hiatus"
      insert <- values[["hiatus"]]
      insert <- sort(insert)
      for (i in length(insert):1){
        if(insert[i] < length(cont[,1]) && insert[i] >= 1 && insert[i] == trunc(insert[i])){
          cont <- rbind(cont[1:insert[i],], label = rep(0,length(cont[1,])),
                         cont[(insert[i]+1):length(cont[,1]),])
          if(length(insert)== 1){
            row.names(cont)[insert[i]+1] <- label
          }else{
            row.names(cont)[insert[i]+1] <- paste0(label,".",(length(insert)+1 - i))
          }
        }else{
          warning("Some positions are outside the range. Ex: ",insert[i])
        }
      }
    }

    plot_EPPM(cont, permute=input$permute, show = show)
    values[["seriograph"]] <- recordPlot()
  })

  observeEvent(input$tri,{
    order <- input$seq2

    values[["rm"]] <- input$rm

    if(is.na(order)[1] || ( length(order) == 1 && order[1] == "Hiatus" )){
      removeModal()
    }else{
      hiatus <- which(order == "Hiatus")
      for (i in 1:length(hiatus)){
        hiatus[i] <- hiatus[i] - 1 * i
      }

      values[["hiatus"]] <- hiatus

      order <- order[order != "Hiatus"]

      new_order <- rep(NA,length(order))
      for(i in 1:length(order)){
        new_order[i] <- which(labels(values[["data1"]])[[1]] == order[i])
      }
      values[["new_order"]] <- new_order
    }
    removeModal()
  })



  observeEvent(input$data,{
    showModal(import1Modal())
    shinyjs::hide(id = "myBox")
  })

  observeEvent(input$Finishimport,{
    values[["data1"]] <- values[["D1"]]
    values[["new_order"]] <- labels(values[["D1"]])[[1]]
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



}
