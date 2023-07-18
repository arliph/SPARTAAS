validation <- function(D1,D2,valid){

  if(class(D2)[1] != "dist"){
    #matrix ?
    if(!is.matrix(D2)){
      if(!is.data.frame(D2)){return("stratigraphic connection data must be a data.frame.")}
      #number of col  ?
      if(valid == 1){
        if(length(D2[1,])!=2){
          return(
            "Stratigraphic connection must have 2 colunms. Change data type parameter or try to inverse contingency table and stratigraphic connection."
            )}
      }else{
        if(length(D2[1,])!=3){
          return("Time range data must have 3 colunms. Change data type parameter or try to inverse contingency table and time range data.")}
      }
      if(valid == 1){
        D2 <- adjacency(D2)
      }else{
        if(!is.numeric(D2[,2]) || !is.numeric(D2[,3])){return("The upper and lower limit must be numeric.")}
        D2 <- overlap(D2)
      }
    }
    #diag = 0 ?
    if(sum(diag(D2)) != 0){
      return("D2: The diagonal must be equal to 0.")
    }
    #symetric ?
    if(!isSymmetric(D2)){
      return("D2 is not symmetric.")
    }
  }

  if(class(D1)[1] != "dist") {
    #matrix ?
    if(min(D1) < 0){
      return("First source is a contingency table and can not contain negative value.")
    }
    if(!is.matrix(D1)){

      D1 <- as.matrix(CAdist(D1,nCP="max",graph=FALSE))
    }
    #diag = 0 ?
    if(sum(diag(D1)) != 0){
      stop("D1: The diagonal must be equal to 0.")
    }
    #symetric ?
    if(!isSymmetric(D1)){
      stop("D1 is not symmetric.")
    }
  }

  if(!dim(D1)[1] == dim(D1)[2]){return("D1 not a square matrix.")}
  if(!dim(D2)[1] == dim(D2)[2]){return("D2 not a square matrix.")}
  if(!length(D1) == length(D2)){return("Contingency table and the second source have not the same size.")}
  #dist (positive matrix) ?
  if(min(D1) < 0){return("D1 is not positive.")}
  if(min(D2) < 0){return("D2 is not positive.")}

  return("")
}

generateClone <- function(D1,D2,i,j){
  D1bis <- rbind(D1,clone=D1[i,])
  D1bis <- cbind(D1bis,clone=c(D1[i,],0))

  D2bis <- rbind(D2,clone=D2[j,])
  D2bis <- cbind(D2bis,clone=c(D2[j,],0))

  return(list(D1bis,D2bis))
}

temp_cov <- function (df) {
  indice <- function(inf1,sup1,inf2,sup2){
    if(sup1 < inf1 || inf2 > sup2){
      return(99)
    }else{
      inf1 <- as.numeric(inf1)
      sup1 <- as.numeric(sup1)
      inf2 <- as.numeric(inf2)
      sup2 <- as.numeric(sup2)
      total_overlap_inf <- min(inf1,inf2)
      total_overlap_sup <- max(sup1,sup2)
      total_overlap <- total_overlap_sup - total_overlap_inf
      within_overlap_inf <- max(inf1,inf2)
      within_overlap_sup <- min(sup1,sup2)
      within_overlap <- within_overlap_sup - within_overlap_inf
      i <- within_overlap / total_overlap
      return(i)
    }
  }
  recouvrement_indice_matrix <- function(df){
    df <- as.data.frame(df)
    overlap_matrix <- matrix(nrow = length(df[,1]), ncol = length(df[,1]))
    for(i in 1:length(df[,1])){
      for(j in 1:length(df[,1])){
        if(i != j){
          ens1 <- df[i,]
          ens2 <- df[j,]
          ind <- indice(ens1[1,2],ens1[1,3],ens2[1,2],ens2[1,3])
          overlap_matrix[i,j] <- overlap_matrix[j,i] <- ind
        }
      }
    }
    diag(overlap_matrix) <- 1
    return(overlap_matrix)
  }

  O <- recouvrement_indice_matrix(df[,1:3])

  D <- O + 1
  D <- 2 - D
  D <- D / 2

  dimnames(D) <- list(df[,1],df[,1])

  return(D)
}

let <- function(alphabet) function(i) {
  base10toA <- function(n, A) {
    stopifnot(n >= 0L)
    N <- length(A)
    j <- n %/% N
    if (j == 0L) A[n + 1L] else paste0(Recall(j - 1L, A), A[n %% N + 1L])
  }
  vapply(i-1L, base10toA, character(1L), alphabet)
}

Int2Alpha <- function(int){
  first <- trunc(int)
  l1 <- let(LETTERS)(first)
  second <- arrondi((int - trunc(int))*10,2)
  l2 <- c()
  for(i in 1:length(second)){
    if(second[i] == 0){
      l2 <- c(l2,"")
    }else{
      f <- nchar(second[i])
      if(f == 1){
        l2 <- c(l2,as.character(second[i]))
      }
      if(f >= 3){
        new <- unlist(strsplit(as.character(second[i]),".",fixed=T))
        new <- paste0(new[1],new[2])
        l2 <- c(l2,new)
      }
    }
  }
  return(paste0(l1,l2))
}

make.chars <- function(length.out, case, n.char = NULL) {
  if(is.null(n.char))
    n.char <- ceiling(log(length.out, 26))
  m <- sapply(n.char:1, function(x) {
    rep(rep(1:26, each = 26^(x-1)) , length.out = length.out)
  })
  m.char <- switch(case,
                   'lower' = letters[m],
                   'upper' = LETTERS[m]
  )
  m.char <- LETTERS[m]
  dim(m.char) <- dim(m)

  apply(m.char, 1, function(x) paste(x, collapse = ""))
}


get.letters <- function(length.out, case = 'upper'){
  max.char <- ceiling(log(length.out, 26))
  grp <- rep(1:max.char, 26^(1:max.char))[1:length.out]
  unlist(lapply(unique(grp), function(n) make.chars(length(grp[grp == n]), case = case, n.char = n)))
}

Alpha2Int <- function(alpha){
  letter <- get.letters(800)
  first <- c()
  second <- c()
  for(i in 1:length(alpha)){
    first <- c(first,str_extract(alpha[i], "[aA-zZ]+"))
    tmp <- str_extract(alpha[i], "[0-9]+")
    if(is.na(tmp)){tmp <- ""}
    second <- c(second,tmp)
  }
  n <- nchar(max(second))
  second <- as.numeric(paste0("0.",second))

  second[is.na(second)] <- 0
  for(i in 1:length(alpha)){
    first[i] <- which(letter == first[i])
  }
  first <- as.numeric(first)
  return(first+second)
}

CAdist <- function(df, nCP = NULL, graph = TRUE){
  #===============================================================================
  # GENERATE DISTANCE MATRIX WITH CORRESPONDENCE ANALYSiS (CA) (Fr: Analyse Factorielle des Correspondances)
  #
  # Authors : A.COULON
  # Last update : 23 october 2018
  #
  # Arguments:
  # df	      The pottery contingence table
  # nCP       The number of principal component to keep
  # graph     Logical for print or not the plot of the CA
  #
  #===============================================================================

  if(!is.table(df)){
    if(is.matrix(df)){df <- as.data.frame(df)}
    if(!is.data.frame(df)){stop("df is not a data frame.")}
  }

  if(FALSE %in% (df == trunc(df))){stop("The data frame must contain integers.")}
  max <- min( nrow(df)-1,ncol(df)-1)

  if(nCP == "max"){nCP <- max}
  if(nCP > max){nCP <- max}

  df.CA <- CA(df, ncp = nCP, graph=graph)
  dist <- dist(df.CA$row$coord)/max(dist(df.CA$row$coord))
  return(as.matrix(dist))
}


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

plot_EPPM <- function(x,show,permute,col_weight) {

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
    theme_grey(base_size = 15) + geom_segment(aes(x = 0, y = -.04, xend = 0, yend = .04),
                                              linetype = "blank") +
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

elbow_finder <- function(x_values, y_values) {
  # Max values to create line
  max_x_x <- max(x_values)
  max_x_y <- y_values[which.max(x_values)]
  max_y_y <- max(y_values)
  max_y_x <- x_values[which.max(y_values)]
  max_df <- data.frame(x = c(max_y_x, max_x_x), y = c(max_y_y, max_x_y))

  # Creating straight line between the max values
  fit <- lm(max_df$y ~ max_df$x)

  # Distance from point to line
  distances <- c()
  for(i in 1:length(x_values)) {
    distances <- c(distances,
      abs( coef(fit)[2] * x_values[i] - y_values[i] + coef(fit)[1] ) / sqrt( coef(fit)[2]^2 + 1^2))
  }
  return(distances)
}

corCriterion <- function(tree,D1,D2) {
  d2 <- cophenetic(tree)
  res <- abs(cor(as.dist(D1),d2) - cor(as.dist(D2),d2))
  return(res)
}

server <- function(input, output, session) {
  #------------------------------------------------------#
  #  modalDialog wait message
  #------------------------------------------------------#
  waitModal <- function() {
    modalDialog(title="Initialisation", size = c("s"),
                div(style="height:25px;"),
                h4("Wait . . . ", style="margin-left:100px;"),
                hr(style="border-color: #222222;"),
                footer = NULL
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
                helpText("Note: Drag and drop clusters to sort them.
                         Drag unused row into the trash to remove them from the seriograph.
                         You can insert 'Hiatus' by picking it up and dropping it off where you want."),
                orderInput('seq2',span(icon("sort-amount-up",lib = "font-awesome"),'Sort the cluster'),
                           items = Int2Alpha(sort(unique(values[["cluster"]]))),
                           placeholder = 'There must be at least one item', connect = "rm",
                           item_class = 'success'),
                hr(style="border-color: #b2b2b2;"),
                column(4,
                       orderInput('hiatus', span(icon("plus",lib = "font-awesome"),'Add'), items = c("Hiatus"),
                                  connect = 'seq2', as_source = TRUE, item_class = 'warning', width = '70px')
                ),
                column(8,
                       orderInput('rm', span(icon("trash",lib = "font-awesome"),'Remove from seriograph'),
                                  items = NULL,connect = "seq2", placeholder = 'Drag the items here to delete them',
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
  #------------------------------------------------------#
  #  modalDialog subdivise
  #------------------------------------------------------#
  subModal1 <- function(failed = FALSE) {
    modalDialog(title="divide a cluster:",size=c("s"),
                div(style="height:15px;"),
                h3("Subdivide a cluster",style="text-align:center;"),hr(style="border-color: #222222;"),
                h4("Select one cluster to divide!"),
                selectInput("where","Cluster to divide:",
                            choices = as.list(sort(Int2Alpha((values[["div"]])))),
                            selected = 1,
                            multiple = FALSE),
                footer = tagList(
                  modalButton("Cancel"),
                  actionButton("sub1", "Next")
                )
    )
  }
  #------------------------------------------------------#
  #  modalDialog subdivise2
  #------------------------------------------------------#
  subModal2 <- function(failed = FALSE) {#
    modalDialog(title="divide a cluster:",size=c("s"),
                div(style="height:15px;"),
                h3("Subdivide a cluster",style="text-align:center;"),hr(style="border-color: #222222;"),
                h4("How many sub-class do you want ?"),
                sliderInput("sk",label = "Number:" ,min = 2,
                    max = length(values[["cluster_ori"]][values[["cluster_ori"]] == Alpha2Int(input$where)])-1,
                    step = 1, value = 2),
                footer = tagList(
                  modalButton("Cancel"),
                  actionButton("back", "Back"),
                  actionButton("sub2", "Finish")
                )
    )
  }
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
                  actionButton("import2", icon("arrow-right",lib = "font-awesome"))
                )
    )
  }

  #------------------------------------------------------#
  #  modalDialog import
  #------------------------------------------------------#
  import2Modal <- function(failed = FALSE) {
    modalDialog(title="Step 2:",size=c("l"),
                div(style="height:15px;"),
                h3("Step 2",style="text-align:center;"),hr(style="border-color: #222222;"),
                h4("Upload file temporal data"),
                radioButtons("datatype", "Type",
                             choices = c("Stratigraphic connection" = 1,
                                         "Time range data" = 2),
                             selected = 1),

                fileInput("fileD2", span(icon("import", lib = "glyphicon"),"Choose CSV File"),
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
                           checkboxInput("headerD2", "Header (colnames)", TRUE),
                           checkboxInput("rownamesD2", "rownames", FALSE)
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
                                        choices = c("Dot" = '.',
                                                    "Comma" = ","),
                                        selected = '.')
                    )
                  )
                )),
                tableOutput("contentsD2"),
                if (failed)
                  div(tags$b("Your data must be match with the example above ! Change parameter or file.",
                             style = "color: red;")),
                footer = tagList(
                  modalButton("Cancel"),
                  actionButton("back1", icon("arrow-left",lib = "font-awesome")),
                  actionButton("import3", icon("arrow-right",lib = "font-awesome"))
                )
    )
  }

  #------------------------------------------------------#
  #  modalDialog import
  #------------------------------------------------------#
  import3Modal <- function(failed = FALSE) {
    modalDialog(title="Step 3:",size=c("l"),
                div(style="height:15px;"),
                h3("Step 3",style="text-align:center;"),hr(style="border-color: #222222;"),
                h4("Visualisation and confirm"),
                tableOutput("contentsD1"),

                tableOutput("contentsD2"),


                v <- validation(values[["D1"]],values[["D2"]],input$datatype),


                if(v == ""){
                  div(tags$b("Valid data!", style = "color: green;"))
                }else{
                  div(tags$b(paste("Invalid data!",v), style = "color: red;"))
                },

                if (failed)
                  div(tags$b("Your data must be match with the example above ! Change parameter or file.",
                             style = "color: red;")),
                footer = tagList(
                  modalButton("Cancel"),
                  actionButton("back2", icon("arrow-left",lib = "font-awesome")),
                  actionButton("Finishimport", "Finish")
                )
    )
  }

  #------------------------------------------------------#
  #  modalDialog select axis
  #------------------------------------------------------#
  pre_selectaxis <- function() {
    modalDialog(title="Select number of axis:",size=c("l"),
                div(style="height:15px;"),
                h3("Select number of axis:",style="text-align:center;"),hr(style="border-color: #222222;"),

                HTML('<p>Our first source of information is a contingency table. We need to construct a distance matrix from this table. We perform a correspondence analysis (CA) on the contingency table and then use the distances of the components (chi-square metric).</p>
                     <p>It is necessary to choose the number of axes (CA components) to be used to construct the distance matrix.</p>
                     <p>By examining the eigenvalues, it is possible to determine the number of principal axes to be considered. The eigenvalues correspond to the amount of information retained by each axis.</p>
                     <p>A heuristic method is often used: one constructs the graph of eigenvalues sorted in descending order and retains the eigenvalues (and associated principal components) preceding the first "elbow".</p>
                     <p>Another approach is based on maintaining an approximation quality of the original data, measured by the explained inertia rate (e.g. 85%). As many axes as necessary are chosen so that the sum of the corresponding eigenvalues exceeds the target inertia rate.</p>'),
                footer = tagList(
                  actionButton("nextCA", "OK")
                )
    )
  }
  selectaxis <- function() {
    modalDialog(title="Select number of axes:",size=c("l"),
                div(style="height:15px;"),
                h3("Select number of axis:",style="text-align:center;"),hr(style="border-color: #222222;"),

                plotlyOutput("axisCAplot"),
                sliderInput("axisNumber",label = "Number:" ,min = 1,
                            max = min( nrow(values[["data1"]])-1,ncol(values[["data1"]])-1),
                            step = 1, value = min( nrow(values[["data1"]])-1,ncol(values[["data1"]])-1)),

                footer = tagList(
                  modalButton("Cancel"),
                  actionButton("aCA", "OK")
                )
    )
  }
  selectaxis2 <- function() {
    modalDialog(title="Select number of axis:",size=c("l"),
                div(style="height:15px;"),
                h3("Select number of axis:",style="text-align:center;"),hr(style="border-color: #222222;"),

                plotlyOutput("axisCAplot"),
                sliderInput("axisNumber",label = "Number:" ,min = 1,
                            max = min( nrow(values[["data1"]])-1,ncol(values[["data1"]])-1),
                            step = 1, value = min( nrow(values[["data1"]])-1,ncol(values[["data1"]])-1)),

                footer = tagList(
                  actionButton("aCA", "OK")
                )
    )
  }

  values <- reactiveValues()
  values[["nCP"]] <- "max"
  values[["cluster"]] <- c()
  values[["ADD"]] <- FALSE
  values[["hiatus"]] <- NULL
  values[["sortperiodtest"]] <- FALSE
  values[["D_alpha"]] <- NULL
  values[["valid"]] <- 1

  values[["alphaplot"]] <- NULL
  values[["dendrogram"]] <- NULL
  values[["seriograph"]] <- NULL
  values[["CA"]] <- NULL
  values[["truecolor"]] <- NULL
  values[["order"]] <- NULL
  values[["table_serio"]] <- NULL
  shinyjs::hide("toggle")
  shinyjs::show("nodata")

  observeEvent(input$axe1,{
    values[["axe1"]] <- input$axe1
  })

  observeEvent(input$axe2,{
    values[["axe2"]] <- input$axe2
  })

  observeEvent(input$width1,{
    values[["width1"]] <- input$width1
  })
  observeEvent(input$height1,{
    values[["height1"]] <- input$height1
  })

  observeEvent(input$width2,{
    values[["width2"]] <- input$width2
  })
  observeEvent(input$height2,{
    values[["height2"]] <- input$height2
  })

  observeEvent(input$width3,{
    values[["width3"]] <- input$width3
  })
  observeEvent(input$height3,{
    values[["height3"]] <- input$height3
  })

  observeEvent(input$width4,{
    values[["width4"]] <- input$width4
  })
  observeEvent(input$height4,{
    values[["height4"]] <- input$height4
  })

  observeEvent(input$width5,{
    values[["width5"]] <- input$width5
  })
  observeEvent(input$height5,{
    values[["height5"]] <- input$height5
  })

  observeEvent(input$width6,{
    values[["width6"]] <- input$width6
  })
  observeEvent(input$height4,{
    values[["height6"]] <- input$height6
  })

  observeEvent(input$method,{
    if(input$method == 1){values[["method"]] <- "ward.D"}
    if(input$method == 2){values[["method"]] <- "ward.D2"}
    if(input$method == 3){values[["method"]] <- "single"}
    if(input$method == 4){values[["method"]] <- "complete"}
    if(input$method == 5){values[["method"]] <- "average"}
    if(input$method == 6){values[["method"]] <- "mcquitty"}
    #if(input$method == 7){values[["method"]] <- "median"}
    #if(input$method == 8){values[["method"]] <- "centroid"}

    values[["cluster"]] <- values[["cluster_ori"]]
    values[["ADD"]] <- FALSE
    values[["div"]] <- as.numeric(names(table(values[["cluster_ori"]]))[table(values[["cluster_ori"]])>1])
    enable("subdivise")
    values[["hiatus"]] <- NULL
    values[["sortperiodtest"]] <- FALSE

    values[["step"]] <- NULL
    values[["which"]] <- NULL
    values[["where"]] <- NULL
    values[["order"]] <- NULL
    values[["firstClustSett"]] <- TRUE
  })


  #raw stratigraphic data (Network)
  values[["data2"]] <- values[["raw"]] <- data.frame(
    ens = c("AI09","AI08","AI07","AI06","AI05","AI04","AI03",
            "AI02","AI01","AO05","AO04","AO03","AO02","AO01","APQR03","APQR02","APQR01"),
    linked_at = c("AI08,AI06","AI07","AI04","AI05","AI01","AI03","AI02","","","AO04","AO03",
            "AO02,AO01","","","APQR02","APQR01","")
  )
  #contingency table
  values[["data1"]] <- values[["cont"]] <- data.frame(
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

  output$data1 <- renderTable(rownames = TRUE, bordered = TRUE,digits = 0,{
    values[["data1"]]
  })

  output$data2 <- renderTable(bordered = TRUE,{
    values[["data2"]]
  })

  output$data1f <- renderTable(rownames = TRUE, bordered = TRUE,digits = 0,{
    data.frame(
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
  })

  output$data2f <- renderTable(bordered = TRUE,{
    data.frame(
      ens = c("AI09","AI08","AI07","AI06","AI05","AI04","AI03",
              "AI02","AI01","AO05","AO04","AO03","AO02","AO01","APQR03","APQR02","APQR01"),
      linked_at = c("AI08,AI06","AI07","AI04","AI05","AI01","AI03","AI02","","","AO04","AO03",
                    "AO02,AO01","","","APQR02","APQR01","")
    )
  })

  output$data3f <- renderTable(bordered = TRUE,digits=0,{
    data.frame(
      ens = c("AI09","AI08","AI07","AI06","AI05","AI04","AI03",
              "AI02","AI01","AO05","AO04","AO03","AO02","AO01","APQR03","APQR02","APQR01"),
      inf = c(400,420,500,590,600,550,400,750,700,750,900,800,700,800,850,1000,1050),
      sup = c(500,500,600,800,800,700,600,800,900,850,1100,950,800,900,900,1100,1100)
    )
  })

  output$axisCAplot <- renderPlotly({
    cont <- values[["data1"]]
    tmp <- CA(cont, graph=FALSE)

    x <- paste0(rownames(tmp$eig)," (",arrondi(tmp$eig[,3],1),"%)")
    Eigen_value <- tmp$eig[,1]
    df2 <- data.frame(x,Eigen_value)
    df2 <- within(df2,Principal_Component <- factor(x,levels = x))

    plot_ly(df2,
            x = ~Principal_Component,
            y = ~Eigen_value,
            type = "bar"
    )
  })

  output$timerange <- renderPlotly({
    data <- values[["data2"]]
    data <- data[order(data[,1]),]
    res <- timerange(data = data,cluster = values[["cluster"]])
    values[["order"]] <- res$order
    res$plot
  })


  output$selectAlpha <- renderPlot({
    if(values[["firstClustSett"]]){
      showModal(pre_selectaxis())
    }else{
      showModal(waitModal())
      shinyjs::show("toggle")
      shinyjs::hide("nodata")
      progress1 <- Progress$new(session, min = 0, max = 2)
      on.exit(progress1$close())
      progress1$set(message = 'Initialisation and calcul of the default alpha ..')

      raw <- values[["data2"]]
      cont <- values[["data1"]]
      if(values[["valid"]] == 1){
        D2 <- adjacency(raw)
      }else{
        D2 <- overlap(raw)
      }
      D1 <- CAdist(cont,nCP=values[["nCP"]],graph=TRUE)
      values[["CA"]] <- recordPlot()
      progress1$set(value = 1)
      sa <- hclustcompro_select_alpha(D1, D2, method = values[["method"]],resampling=FALSE, ite= 5)
      values[["alphaplot"]] <- sa$alpha.plot
      print(sa$alpha.plot)

      progress1$set(value = 2)

      updateSliderInput(session, "alpha", value = sa$alpha)
      updateSliderInput(session, "k", value = values[["nb_classes"]],min=2,
                        max = length(values[["raw"]][,1])-1,step = 1)

      removeModal()
    }
  })

  observeEvent(input$aCA,values[["firstClustSett"]] <- FALSE)

  output$dendrogramme <- renderPlot({
    raw <- values[["data2"]]
    cont <- values[["data1"]]
    if(values[["valid"]] == 1){
      D2 <- adjacency(raw)
    }else{
      D2 <- overlap(raw)
    }
    D1 <- CAdist(cont,nCP=values[["nCP"]],graph=FALSE)
    alpha <- input$alpha

    #create mixing matrix
    D1 <- as.dist(D1)
    D2 <- as.dist(D2)
    D1 <- as.matrix(D1)
    D2 <- as.matrix(D2)
    #sort labels
    D1.Tmp <- D1[sort(labels(D1)[[1]]),sort(labels(D1)[[1]])]
    D2.Tmp <- D2[sort(labels(D2)[[1]]),sort(labels(D2)[[1]])]
    D1 <- as.dist(D1.Tmp)
    D2 <- as.dist(D2.Tmp)
    #normalization
    if(max(D1) != 0){
      D1 <- D1 / max(D1)
    }
    if(max(D2) != 0){
      D2 <- D2 / max(D2)
    }
    #mixt matrix
    dist.mixt <- (alpha) * D1 + (1-alpha) * D2
    dist.mixt <- as.dist(dist.mixt)
    values[["D_alpha"]] <- dist.mixt

    #CAH
    tree <- hclust(dist.mixt, method = values[["method"]])
    values[["tree"]] <- tree

    #wss
    D <- as.matrix(dist.mixt)
    n <- length(D[,1])
    if(length(D[,1]) > 20){n <- 20}
    values[["n"]] <- n

    wss <- c()
    sil <- c()
    for ( i in 2:(n-1)){
      cutree <- cutree(tree,i)
      res <- cluster.stats(D,cutree)
      wss <- c(wss,res$within.cluster.ss)
      sil <- c(sil,res$avg.silwidth)
    }
    values[["wss_v"]] <- wss
    values[["sil_v"]] <- sil
    Within_Sum_of_Square <- rbind(seq(2,(n-1)),wss)
    color_wss <- elbow_finder(seq(2,(n-1)),wss)

    dataplot <- t(Within_Sum_of_Square)
    dataplot <- as.data.frame(dataplot)
    colnames(dataplot) <- c("Number_of_groups","Within_Sum_of_Square")
    q <- ggplot(dataplot,aes(x=Number_of_groups, y=Within_Sum_of_Square)) +
      geom_point(aes(color=color_wss),size=2) + scale_colour_gradientn(colours=c("brown1","chartreuse3")) +
      geom_line(linetype = 3) +
      scale_x_continuous(breaks = dataplot$Number_of_groups ,
                         labels = dataplot$Number_of_groups ) +
      ylim(0,max(wss)+1) +
      labs(title = "Within Sum of Square" , subtitle = "Partition evaluation") +
      theme(legend.position="none")
    values[["wss"]] <- q

    #silhouette
    Average_Sil_Width <- rbind(seq(2,(n-1)),sil)
    dataplot <- t(Average_Sil_Width)
    dataplot <- as.data.frame(dataplot)
    dataplot <- na.omit(dataplot)
    colnames(dataplot) <- c("Number_of_groups","Average_Sil_Width")
    p <- ggplot(dataplot,aes(x=Number_of_groups, y=Average_Sil_Width)) +
      geom_point(aes(color=Average_Sil_Width),size=2) +
      scale_colour_gradientn(colours=c("brown1","chartreuse3")) +
      geom_line(linetype = 3) + ylim(0,1) +
      scale_x_continuous(breaks = dataplot$Number_of_groups,
                         labels = dataplot$Number_of_groups) +
      labs(title = "Average Silhouette Width", subtitle = "Partition evaluation") +
      theme(legend.position="none")
    values[["silhouette"]] <- p

    #dendrogram
    title <- paste("Cluster Dendrogram with alpha =",alpha)
    subtitle <- paste("The mixing matrix is made of ",alpha * 100,"% of D1 and ",(1-alpha) * 100,"% of D2.")
    plot(tree,
         h=-1, main = title, sub = subtitle, xlab=paste0("dist.mixt, method: \"",values[["method"]],"\""))
    if(length(input$k) == 0){
      k <- 2
    }else{
      k <- input$k
    }
    cluster <- cutree(tree, k = k)
    order <- unique(cluster[tree$order])
    for(i in 1:length(order)){
      if(i != order[i] && i < order[i]){
        cluster[cluster == i] <- 0
        cluster[cluster == order[i]] <- i
        cluster[cluster == 0] <- order[i]
        order <- unique(cluster[tree$order])
      }
    }

    if(values[["ADD"]] == FALSE && values[["sortperiodtest"]] == FALSE){
      values[["clusterserio"]] <- values[["cluster"]] <- cluster
      values[["cluster_ori"]] <- cluster
      values[["div"]] <- as.numeric(names(table(values[["cluster_ori"]]))[table(values[["cluster_ori"]])>1])
    }
    rect.hclust(tree, k = k, border = rainbow_hcl(k, c = 80, l = 60, start = 0, end = 360*(k-1)/k))
    col <- rainbow_hcl(length(unique(cluster)), c = 80, l = 60, start = 0,
                       end = 360*(length(unique(cluster))-1)/(length(unique(cluster))))
    for(i in 1:length(unique(cluster))){
      x <- mean(which(cluster[tree$order] == sort(unique(cluster))[i]))
      mtext(paste0(LETTERS[i]), side = 1, line = .5, at = x, col = col[i], cex=1.5)
    }
    mtext("Cluster:",side = 1, line = .5, at = 0, col = "#767676", cex=1.2)
    mtext("Order:",side = 1, line = 1.5, at = 0, col = "#767676", cex=1.2)
    if(values[["sortperiodtest"]]){

    }
    #if subdivide
    if(values[["ADD"]]){
      #rectangle
      for(i in 1:length(values[["which"]])){
        #col <- rainbow_hcl(length(unique(values[["cluster"]])), c = 80, l = 60, start = 0,
        #            end = 360*(length(unique(values[["cluster"]]))-1)/(length(unique(values[["cluster"]]))))
        #index <- which(unique(values[["cluster"]][tree$order]) > values[["where"]][[i]])
        rect.hclust(tree,h = tree$height[values[["step"]][[i]]+1],which = values[["which"]][[i]],
                    border = values[["truecolor"]][values[["where"]][[i]]])#col[index])
      }
      #cluster white erase
      for(i in 1:length(unique(values[["cluster_ori"]]))){
        x <- mean(which(values[["cluster_ori"]][tree$order] == sort(unique(values[["cluster_ori"]]))[i]))
        mtext(Int2Alpha(sort(unique(values[["cluster_ori"]]))[i]), side = 1, line = .5, at = x, col = "white",
              cex=1.5)
      }
      #cluster col
      for(i in 1:length(unique(values[["cluster"]]))){
        x <- mean(which(values[["cluster"]][tree$order] == sort(unique(values[["cluster"]]))[i]))
        j <- trunc(sort(unique(values[["cluster"]]))[i])
        mtext(Int2Alpha(sort(unique(values[["cluster"]]))[i]), side = 1, line = .5, at = x,
              col = values[["truecolor"]][j], cex=1.5)
      }
      for(i in 1:length(unique(values[["clusterserio"]]))){
        #if not remove  period
        if(!nchar(sort(unique(values[["clusterserio"]]))[i]) == 5){
          x <- mean(which(values[["clusterserio"]][tree$order] == sort(unique(values[["clusterserio"]]))[i]))
          #y <- which(unique(values[["clusterserio"]][tree$order]) == sort(unique(values[["clusterserio"]]))[i])
          j <- trunc(sort(unique(values[["cluster"]]))[i])
          mtext(paste0(sort(unique(values[["clusterserio"]]))[i]), side = 1, line = 1.5, at = x,
                col = "#333333", cex=1.3)
        #if remove elt  period
        }else{
          x <- mean(which(values[["clusterserio"]][tree$order] == sort(unique(values[["clusterserio"]]))[i]))
          #y <- which(unique(values[["clusterserio"]][tree$order]) == sort(unique(values[["clusterserio"]]))[i])
          j <- trunc(sort(unique(values[["cluster"]]))[i])
          mtext("Rm", side = 1, line = 1.5, at = x, col = "#333333", cex=1.3)
        }
      }
    #not subdivide
    }else{
      for(i in 1:length(unique(values[["clusterserio"]]))){
        #if not remove  period
        if(!nchar(sort(unique(values[["clusterserio"]]))[i]) == 5){
          x <- mean(which(values[["clusterserio"]][tree$order] == sort(unique(values[["clusterserio"]]))[i]))
          #y <- which(unique(values[["clusterserio"]][tree$order]) == sort(unique(values[["clusterserio"]]))[i])
          j <- trunc(sort(unique(values[["cluster"]]))[i])
          mtext(paste0(sort(unique(values[["clusterserio"]]))[i]), side = 1, line = 1.5, at = x,
                col = "#333333", cex=1.3)
        #if remove elt  period
        }else{
          x <- mean(which(values[["clusterserio"]][tree$order] == sort(unique(values[["clusterserio"]]))[i]))
          #y <- which(unique(values[["clusterserio"]][tree$order]) == sort(unique(values[["clusterserio"]]))[i])
          j <- trunc(sort(unique(values[["cluster"]]))[i])
          mtext("Rm", side = 1, line = 1.5, at = x, col = "#333333", cex=1.3)
        }
      }
    }
    values[["dendrogram"]] <- recordPlot()
  })

  output$EPPM <- renderPlot({
    seq <- values[["clusterserio"]]
    hiatus <- values[["hiatus"]]
    cont <- values[["data1"]]
    cont <- cont[sort(labels(cont)[[1]]),]
    cont2 <- mutate(cont, Cluster = factor(seq))
    cont2 <- group_by(cont2, Cluster)
    cont2 <- summarise_all(cont2, sum)
    cont2 <- as.data.frame(cont2)
    row.names(cont2) <- as.character(unique(sort(seq)))
    cont2 <- cont2[,-1]
    remove <- which(nchar(rownames(cont2)) == 5)
    #remove <- which(rownames(cont2) >= 900)
    if(length(remove) >= 1){
      cont2 <- cont2[-remove,]
    }
    cont2 <- cont2[order(as.numeric(rownames(cont2)),decreasing = T),]
    rownames(cont2) <- as.character(paste0("Cluster.",rownames(cont2)))
    #hiatus
    if(length(hiatus) > 0 && !is.na(hiatus)){
      #insert hiatus
      label <- "Hiatus"
      insert <- values[["hiatus"]]
      insert <- sort(length(cont2[,1]) - insert)
      for (i in length(insert):1){
        if(insert[i] < length(cont2[,1]) && insert[i] >= 1 && insert[i] == trunc(insert[i])){
          cont2 <- rbind(cont2[1:insert[i],], label = rep(0,length(cont2[1,])),
                         cont2[(insert[i]+1):length(cont2[,1]),])
          if(length(insert)== 1){
            row.names(cont2)[insert[i]+1] <- label
          }else{
            row.names(cont2)[insert[i]+1] <- paste0(label,".",(length(insert)+1 - i))
          }
        }else{
          warning("Some positions are outside the range. Ex: ",insert[i])
        }
      }
    }
    show <- c("both","frequency","EPPM")[as.numeric(input$show)]
    values[["table_serio"]] <- cont2
    plot_EPPM(cont2, permute=input$permute, col_weight = input$col_weight, show = show)
    values[["seriograph"]] <- recordPlot()
  })

  output$CA <- renderScatterD3({

    updateNumericInput(session, "axe1", max = min( nrow(values[["data1"]])-1,ncol(values[["data1"]])-1))
    updateNumericInput(session, "axe2", max = min( nrow(values[["data1"]])-1,ncol(values[["data1"]])-1))

    res <- explor::prepare_results(CA(values[["data1"]], ncp = min( nrow(values[["data1"]])-1,ncol(values[["data1"]])-1)))
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

    var_data <- CA_var_data(res, xax, yax, var_sup, var_hide, var_lab_min_contrib)

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

  output$WSS <- renderPlot({
    print(values[["wss"]])
  })

  output$silhouette <- renderPlot({
    print(values[["silhouette"]])
  })

  observeEvent(input$aCA,{
    values[["nCP"]] <- input$axisNumber
    removeModal()
  })
  observeEvent(input$nextCA,{
    removeModal()
    showModal(selectaxis2())
  })
  observeEvent(input$sort,{
    showModal(sortPeriod())
  })

  observeEvent(input$tri,{
    remove <- input$rm_order
    remove <- remove[2][remove[2] != "Hiatus"]
    order <- input$seq2_order

    if(is.na(order[2])[1] || ( length(order[2]) == 1 && order[2][1] == "Hiatus" )){
      removeModal()
    }else{
      hiatus <- which(order[2] == "Hiatus")
      for (i in 1:length(hiatus)){
        hiatus[i] <- hiatus[i] - 1 * i
      }

      values[["hiatus"]] <- hiatus
      order <- order[2][order[2] != "Hiatus"]
      remove <- remove[!is.na(remove)]

      if(length(remove) >= 1){
        remove <- Alpha2Int(remove)
      }

      split <- Alpha2Int(order)

      if(length(split)==0){
        removeModal()
        showModal(sortPeriod(TRUE))
      }else{
        cluster <- values[["cluster"]]
        seq <- cluster

        if(length(remove) >= 1){
          for(i in 1:length(remove)){
            removeseq <- which(seq %in% remove[i])
            seq[removeseq] <- seq[removeseq] + as.numeric(paste0("0.00",i))
          }
        }
        c <- 0
        d <- 1
        pelt <- 0
        for (i in 1:length(split)) {
          elt <- as.numeric(split[i])
          felt <- as.numeric(split[i+1])
          if(elt != as.integer(elt)){
            if(as.integer(elt) != as.integer(pelt)){
              c <- c + 1
              d <- 1
            }
            seq[cluster == elt] <- as.numeric(paste0(c,".",d))
            d <- d + 1
          }else{
            c <- c + 1
            seq[cluster == elt] <- c
            d <- 1
          }
          pelt <- elt
        }
        values[["clusterserio"]] <- seq
        values[["sortperiodtest"]] <- TRUE
        removeModal()
      }
    }
  })

  observeEvent(input$subdivise,{
    showModal(subModal1())
  })

  observeEvent(input$sub1,{
    removeModal()
    showModal(subModal2())
  })

  observeEvent(input$sub2,{
    removeModal()
    ok <- TRUE
    cluster <- values[["cluster"]]
    tree <- values[["tree"]]
    where <- Alpha2Int(input$where)
    values[["div"]] <- values[["div"]][values[["div"]] != where]
    if(length(values[["div"]]) == 0){
      disable("subdivise")
    }
    sk <- input$sk
    step <- length(cluster)-1
    while(ok){
      clusters <- cutree(tree,h = tree$height[step])
      if(length(unique(clusters[cluster == where])) == sk){
        ok <- FALSE
      }else{
        step <- step - 1
      }
    }
    #reorder value
    order <- unique(clusters[tree$order])
    for(i in 1:length(order)){
      if(i != order[i] && i < order[i]){
        clusters[clusters == i] <- 0
        clusters[clusters == order[i]] <- i
        clusters[clusters == 0] <- order[i]
        order <- unique(clusters[tree$order])
      }
    }
    nb <- which(values[["cluster_ori"]] == Alpha2Int(input$where))
    nb <- clusters[nb]
    which <- sort(unique(nb))
    #rename value
    resume <- table(clusters[cluster == where])
    c <- 1
    for(i in as.numeric(names(resume))){
      clusters[clusters == i] <- c
      c <- c+1
    }
    c <- c-1
    clusters <- clusters / 10 ^nchar(c)
    for(i in 1:length(cluster)){
      ifelse(cluster[i] == where,
             cluster[i] <- where + clusters[i],
             cluster[i] <- cluster[i])
    }
    values[["clusterserio"]] <- values[["cluster"]] <- cluster
    values[["step"]][[length(values[["step"]])+1]] <- step
    values[["which"]][[length(values[["which"]])+1]] <- which
    values[["where"]][[length(values[["where"]])+1]] <- Alpha2Int(input$where)
    values[["ADD"]] <- TRUE
  })

  observeEvent(input$alpha,{
    values[["cluster"]] <- values[["cluster_ori"]]
    values[["ADD"]] <- FALSE
    values[["div"]] <- as.numeric(names(table(values[["cluster_ori"]]))[table(values[["cluster_ori"]])>1])
    enable("subdivise")
    values[["hiatus"]] <- NULL
    values[["sortperiodtest"]] <- FALSE
    values[["order"]] <- NULL
  })

  observeEvent(input$k,{
    values[["cluster"]] <- values[["cluster_ori"]]
    values[["ADD"]] <- FALSE
    values[["div"]] <- as.numeric(names(table(values[["cluster_ori"]]))[table(values[["cluster_ori"]])>1])
    enable("subdivise")
    values[["hiatus"]] <- NULL
    values[["sortperiodtest"]] <- FALSE
    values[["truecolor"]] <- rainbow_hcl(input$k, c = 80, l = 60, start = 0,
                               end = 360*(input$k-1)/input$k)
    values[["step"]] <- NULL
    values[["which"]] <- NULL
    values[["where"]] <- NULL
    values[["order"]] <- NULL
    values[["table_serio"]] <- NULL
  })

  observeEvent(input$back,{
    removeModal()
    showModal(subModal1())
  })

  observeEvent(input$reset,{
    values[["cluster"]] <- values[["cluster_ori"]]
    values[["ADD"]] <- FALSE
    values[["div"]] <- as.numeric(names(table(values[["cluster_ori"]]))[table(values[["cluster_ori"]])>1])
    enable("subdivise")
    values[["hiatus"]] <- NULL
    values[["sortperiodtest"]] <- FALSE
    values[["step"]] <- NULL
    values[["which"]] <- NULL
    values[["where"]] <- NULL
    values[["order"]] <- NULL
    values[["table_serio"]] <- NULL
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
    if(input$datatype == 2){values[["timerange"]] <- df}
    values[["D2"]] <- df
    return(head(df))
  },rownames=TRUE)

  observeEvent(input$data,{
    showModal(import1Modal())
    shinyjs::hide(id = "myBox")
  })

  observeEvent(input$back1,{
    showModal(import1Modal())
    shinyjs::hide(id = "myBox")
  })

  observeEvent(input$back2,{
    showModal(import2Modal())
    shinyjs::hide(id = "myBox")
  })

  observeEvent(input$import2,{
    showModal(import2Modal())
    shinyjs::hide(id = "myBox")
  })

  observeEvent(input$import3,{
    showModal(import3Modal())
  })

  observeEvent(input$Finishimport,{

    v <- validation(values[["D1"]],values[["D2"]],input$datatype)

    if(!v == ""){
      removeModal()
      showModal(import3Modal(TRUE))
    }else{
      removeModal()
      values[["nCP"]] <- "max"
      values[["cluster"]] <- values[["cluster_ori"]]
      values[["ADD"]] <- FALSE
      values[["div"]] <- as.numeric(names(table(values[["cluster_ori"]]))[table(values[["cluster_ori"]])>1])
      enable("subdivise")
      values[["data1"]] <- values[["D1"]]
      values[["data2"]] <- values[["D2"]]
      values[["D1"]] <- NULL
      values[["D2"]] <- NULL
      values[["PA"]] <- NULL
      values[["sortperiodtest"]] <- FALSE
      values[["clusterserio"]] <- c(1)
      values[["valid"]] <- input$datatype
      values[["firstClustSett"]] <- TRUE
      updateNumericInput(session, "axe1", value = 1)
      updateNumericInput(session, "axe2", value = 2)
    }

  })

  observeEvent(input$data,{
    shinyjs::disable("import2")
  })

  observeEvent(input$import2,{
    shinyjs::disable("import3")
  })

  observeEvent(input$fileD1,{
      shinyjs::enable("import2")
  })

  observeEvent(input$fileD2,{
      shinyjs::enable("import3")
  })

  observeEvent(input$toggle, {
    if(input$toggle %% 2 == 0){
      shinyjs::hide(id = "myBox")
    }else{
      shinyjs::show(id = "myBox")
    }
  })

  observeEvent(input$refresh, {
    session$reload()
  })

  observeEvent(input$axesCA, {
    showModal(selectaxis())
  })

  #EXPORT BUTTON
  output$alpha.pdf <- downloadHandler(
    filename = "alpha.pdf",
    content = function(file) {
      cairo_pdf(file,width=12)
      print(values[["alphaplot"]])
      dev.off()
    }
  )

  output$alpha.png <- downloadHandler(
    filename = "alpha.png",
    content = function(file) {
      png(file,
          width = values[["width1"]],
          height = values[["height1"]]
      )
      print(values[["alphaplot"]])
      dev.off()
    }
  )

  output$dendrogram.pdf <- downloadHandler(
    filename = "dendrogram.pdf",
    content = function(file) {
      cairo_pdf(file,width=12)
      print(values[["dendrogram"]])
      dev.off()
    }
  )

  output$dendrogram.png <- downloadHandler(
    filename = "dendrogram.png",
    content = function(file) {
      png(file,
          width = values[["width2"]],
          height = values[["height2"]]
      )
      print(values[["dendrogram"]])
      dev.off()
    }
  )

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
          width = values[["width3"]],
          height = values[["height3"]]
      )
      print(values[["seriograph"]])
      dev.off()
    }
  )

  output$silplot.pdf <- downloadHandler(
    filename = "silplot.pdf",
    content = function(file) {
      cairo_pdf(file,width=12)
      print(values[["silplot"]])
      dev.off()
    }
  )

  output$silplot.png <- downloadHandler(
    filename = "silplot.png",
    content = function(file) {
      png(file,
          width = values[["width4"]],
          height = values[["height4"]]
      )
      print(values[["silplot"]])
      dev.off()
    }
  )

  output$wss.pdf <- downloadHandler(
    filename = "wss.pdf",
    content = function(file) {
      cairo_pdf(file,width=7)
      print(values[["wss"]])
      dev.off()
    }
  )

  output$wss.png <- downloadHandler(
    filename = "wss.png",
    content = function(file) {
      png(file,
          width = 500,
          height = 500
      )
      print(values[["wss"]])
      dev.off()
    }
  )

  output$silhouette.pdf <- downloadHandler(
    filename = "silhouette.pdf",
    content = function(file) {
      cairo_pdf(file,width=7)
      print(values[["silhouette"]])
      dev.off()
    }
  )

  output$silhouette.png <- downloadHandler(
    filename = "silhouette.png",
    content = function(file) {
      png(file,
          width = 500,
          height = 500
      )
      print(values[["silhouette"]])
      dev.off()
    }
  )

  output$table_serio.csv <- downloadHandler(
    filename = function() {
      paste("table_serio", ".csv", sep = "")
    },
    content = function(file) {
      write.csv2(values[["table_serio"]], file, row.names = TRUE)
    }
  )


  observeEvent(input$seriograph.pdf,{
    showNotification("wait for backup in progress",duration=15)
  })

  output$detail <- renderPlotly({
    raw <- values[["data2"]]
    cont <- values[["data1"]]
    if(values[["valid"]] == 1){
      D2 <- adjacency(raw)
    }else{
      D2 <- temp_cov(raw)
    }
    D1 <- CAdist(cont,nCP=values[["nCP"]],graph=FALSE)
    method <- values[["method"]]
    acc = 2


    METHODS <- c("ward.D", "single", "complete", "average", "mcquitty",
                 "median", "centroid", "ward.D2")
    if (method == "ward") {
      message("The \"ward\" method has been renamed to \"ward.D\"; note new avaible method \"ward.D2\"")
      method <- "ward.D"
    }
    i.meth <- pmatch(method, METHODS)
    if (is.na(i.meth))
      stop("invalid clustering method", paste("", method))
    if (i.meth == -1)
      stop("ambiguous clustering method", paste("", method))

    method <- METHODS[i.meth]
    if(!is.numeric(acc)) {stop("acc must be numeric!")}
    if(acc < 1) {
      acc <- 1;
      message("The minimal value allowed for acc is 1.")
    }

    if(class(D1)[1] == "hclustcompro_cl"){
      if(is.null(D1$D1) || is.null(D1$D2)){
        stop("We don't find dissmilarities in hclustcompro_cl object.")
      }else{
        D2 <- D1$D2
        D1 <- D1$D1
      }
    }else{
      if(is.null(D1) || is.null(D2)){
        stop("We don't find dissmilarities in D1 and D2 arguments.")
      }
      #test
      if(class(D1)[1] != "dist"){
        #matrix ?
        if(!is.matrix(D1)){
          D1 <- as.matrix(CAdist(cont,nCP=values[["nCP"]],graph=FALSE))
        }
        #diag = 0 ?
        if(sum(diag(D1)) != 0){
          stop("D1: in distance matrix the diagonal is equal to 0.")
        }
        #symetric ?
        if(!isSymmetric(D1)){
          stop("D1 is not symmetric.")
        }
      }
      if(class(D2)[1] != "dist"){
        #matrix ?
        if(!is.matrix(D2)){
          if(valid){
            D2 <- adjacency(raw)
          }else{
            D2 <- temp_cov(raw)
          }
        }
        #diag = 0 ?
        if(sum(diag(D2)) != 0){
          stop("D2: in distance matrix the diagonal is equal to 0.")
        }
        #symetric ?
        if(!isSymmetric(D2)){
          stop("D2 is not symmetric.")
        }
      }
      if(!length(D1[,1]) == length(D1[1,])){stop("D1 not a square matrix.")}
      if(!length(D2[,1]) == length(D2[1,])){stop("D2 not a square matrix.")}
      #dist (positive matrix) ?
      if(min(D1) < 0){stop("D1 not a matrix dissimilarity (positive).")}
      if(min(D2) < 0){stop("D2 not a matrix distance (positive).")}
      if(median(D2) != 0 && median(D2) != 1){
        if(median(D1) == 0 || median(D1) == 1){
          message("We switch the two matrices. Stratigraphic data must be in D1 and ceramic data in D2.")
          tmp <- D1
          D1 <- D2
          D2 <- tmp
        }
      }
      #normalization
      if(max(D1) != 0){
        D1 <- D1 / max(D1)
      }
      if(max(D2) != 0){
        D2 <- D2 / max(D2)
      }
      #create mixing matrix
      D1 <- as.dist(D1)
      D2 <- as.dist(D2)
      D1 <- as.matrix(D1)
      D2 <- as.matrix(D2)
      #sort labels
      D1.Tmp <- D1[sort(labels(D1)[[1]]),sort(labels(D1)[[1]])]
      D2.Tmp <- D2[sort(labels(D2)[[1]]),sort(labels(D2)[[1]])]
      D1 <- as.matrix(D1.Tmp)
      D2 <- as.matrix(D2.Tmp)
    }

    D1save <- as.matrix(D1)
    D2save <- as.matrix(D2)

    save_lines <- NULL
    alpha_seq <- seq(0,1,0.01)

    showModal(waitModal())

    withProgress(message = 'Resampling process', value = 0, {
      nb_elt <- length(D1save[,1])
      for(i in 1:nb_elt){
        # update progress bar
        incProgress(1/nb_elt, detail = paste("Doing part", i))
        Sys.sleep(0.1)
        df <- NULL
        for(j in 1:nb_elt){
          if(j != i ){
            list <- generateClone(D1save,D2save,i,j)
            D1 <- list[[1]]
            D2 <- list[[2]]
            dist.dend <- c()
            for(alpha in alpha_seq){
              Mdist <- (alpha)*D1 + (1-alpha)*D2
              mixt.dist <- as.dist(Mdist)
              tree <- fastcluster::hclust(mixt.dist,method=method)
              new <- corCriterion(tree,D1,D2)
              dist.dend <- c(dist.dend,new)
            }
            if(is.null(df)){
              df <- data.frame(alpha = alpha_seq, height = dist.dend)
            }else{
              df <- cbind(df,dist.dend)
            }
          }
        }
        mean <- unlist(lapply(as.data.frame(t(df[-1])),mean))
        if(is.null(save_lines)){
          save_lines <- cbind(df$alpha,mean)
        }else{
          save_lines <- cbind(save_lines,mean)
        }
      }
    })
    colnames(save_lines) <- c("alpha",rownames(D1)[1:nb_elt])
    rownames(save_lines) <- c()

    removeModal()

    p <- plot_ly(width = 760, height = 400)
    for (k in 2:length(colnames(save_lines))) {
      p <- add_lines(p,x=save_lines[,1],y=save_lines[,k],name=colnames(save_lines)[k])
    }
    p
  })

  observe({
    if(2 > 1){
      output$SilPlot <- renderPlot({
        cluster <- values[["cluster"]]
        cls <- sort(unique(cluster))
        for (i in length(cls):1){
          cluster[cluster == cls[i]] <- i
        }
        Dmatrix <- values[["D_alpha"]]
        names <- labels(Dmatrix)
        text <- c("")
        for (i in 1:length(names)){
          text <- paste0(text,i,":",names[i],"\n")
        }
        sil <- silhouette(cluster,Dmatrix)
        col <- rainbow_hcl(length(unique(trunc(values[["cluster"]]))), c = 80, l = 60, start = 0,
                end = 360*(length(unique(trunc(values[["cluster"]])))-1)/(length(unique(trunc(values[["cluster"]])))))
        col <- col[trunc(sort(values[["cluster"]]))]
        par(fig=c(.2,1,0,1), new=TRUE)
        plot(sil, col = col)
        par(fig=c(0,.2,0,1), new=TRUE)
        plot.new()
        mtext(text,3,adj=0,padj=1.1)
        values[["silplot"]] <- recordPlot()
      })
    }
  })


  observeEvent(input$wiki, {
    updateTabsetPanel(session, "navbar", selected = "4")
  })
  observeEvent(input$over, {
    updateTabsetPanel(session, "tabsetperso", selected = "Overview")
  })
  observeEvent(input$clustsett, {
    updateTabsetPanel(session, "tabsetperso", selected = "Clustering settings")
  })
  observeEvent(input$resviz, {
    updateTabsetPanel(session, "tabsetperso", selected = "Results visualization")
  })

  ##dendrogramme click
  output$click_info <- renderPrint({
    if(length(input$k) == 0){
      cat("hcluscompro")
    }else{

      size <- c()

      cluster <- values[["cluster"]]
      nb <- length(unique(cluster))
      tot <- length(cluster)
      cat("The selected partition contains",nb,"clusters among",tot,"observations.\ncluster:\t")
      for (k in unique(values[["cluster"]])){
        size <- c(size,length(cluster[cluster == k]))
        cat(k,"\t")
      }
      cat("\nSize:\t\t")
      for(y in 1:nb){
        cat(size[y],"\t")
      }

    }
  })

  output$formula <- renderUI({
    withMathJax(paste0("$$ \\boldsymbol{D}_\\alpha=\\alpha \\boldsymbol{D}_1+(1-\\alpha) \\boldsymbol{D}_2 $$"))
  })


}
