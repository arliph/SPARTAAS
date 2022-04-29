CAdist <- function(df, nPC = NULL, graph = TRUE){

  #===============================================================================
  # GENERATE DISTANCE MATRIX WITH CORRESPONDENCE ANALYSiS (CA) (Fr: Analyse Factorielle des Correspondances)
  #
  # Authors : A.COULON
  # Last update : 23 october 2018
  #
  # Arguments:
  # df	      The ceramic contingence table
  # nPC       The number of principal component to keep
  # graph     Logical for print or not the plot of the CA
  #
  #===============================================================================

  if(!is.table(df)){
    if(is.matrix(df)){df <- as.data.frame(df)}
    if(!is.data.frame(df)){stop("df is not a data frame.")}
  }

  if(FALSE %in% (df == trunc(df))){stop("The data frame must contain integers.")}
  max <- min( nrow(df)-1,ncol(df)-1)

  if(!is.null(nPC)){
    if(nPC == "max"){nPC <- max}
    if(nPC > max){
      message("nPC must be less than: ",max+1,".")
      nPC <- max
    }
  }else{
    tmp <- FactoMineR::CA(df, graph=graph)

    x <- paste0(rownames(tmp$eig)," (",arrondi(tmp$eig[,3],1),"%)")
    Eigen_value <- tmp$eig[,1]
    df2 <- data.frame(x,Eigen_value)
    df2 <- within(df2,Principal_Component <- factor(x,levels = x))

    p <- plotly::plot_ly(df2,
      x = ~Principal_Component,
      y = ~Eigen_value,
      type = "bar"
    )
    print(p)
    if(interactive()){
      OK <- TRUE
      while(OK){
        nPC <- as.integer(readline("Select the number of axes:"))
        if(is.na(nPC) || nPC < 1 || nPC > max){
          message("must be a integer between 1 and ",max)
        }else{
          OK <- FALSE
        }
      }
    }else{
      stop("Need a nPC value.")
    }

  }
  df.CA <- FactoMineR::CA(df, ncp = nPC, graph = graph)
  dist <- dist(df.CA$row$coord) / max(dist(df.CA$row$coord))
  return(D=as.matrix(dist))
}

adjacency <- function(network){

  #===============================================================================
  # GENERATE DISTANCE MATRIX WITH NETWORK
  #
  # Authors : A.COULON
  # Last update : 23 october 2018
  #
  # Arguments:
  # raw	      The data frame with the stratigraphic element and the relation with the other
  #
  #===============================================================================

  #data.frame ?
  if(!is.data.frame(network)){stop("network data must be a data.frame.")}
  #number of col  ?
  if(length(network[1,])!=2){stop("network must have 2 colunms.")}
  M <- matrix(0,ncol = length(network[,1]),nrow = length(network[,1]),dimnames = list(network[,1],network[,1]))
  diag(M) <- 1
  for(i in 1:length(network[,1])){
    rel <- as.character(network[i,2])
    rel <- strsplit(rel,",")
    for(elt in rel[[1]]){
      j <- which(network[,1] == elt)
      M[i,j] <- 1
      M[j,i] <- 1
    }
  }
  D <- 1 - M
  return(D)
}

overlap <- function(temporal){

  #===============================================================================
  # GENERATE DISSIMILARITY MATRIX WITH TEMPORAL
  #
  # Authors : A.COULON
  # Last update : 23 october 2018
  #
  # Arguments:
  # temporal	      The data frame
  #
  #===============================================================================

  indice <- function(inf1,sup1,inf2,sup2){
    if(inf2 == inf1 && inf1 == sup2 && inf1 == sup1){
      return(1)
    }
    if(sup1 < inf1 || inf2 > sup2){
      return(99)
    }
    #we define the bound of the total overlap as follow : min of lower bounds and max of upper bounds
    total_overlap_inf <- min(inf1,inf2)
    total_overlap_sup <- max(sup1,sup2)

    #The overlap is the difference between the bounds
    total_overlap <- total_overlap_sup - total_overlap_inf

    #we define the bound of the within overlap as follow : max of lower bounds and min of upper bounds
    within_overlap_inf <- max(inf1,inf2)
    within_overlap_sup <- min(sup1,sup2)

    #The overlap is the difference between the bounds
    within_overlap <- within_overlap_sup - within_overlap_inf

    #the overlap index is the ratio between the within overlap and the total overlap
    i <- within_overlap / total_overlap
    return(i)
  }

  df <- as.data.frame(temporal)
  #we create a square matrix of the same size as the number of rows
  overlap_matrix <- matrix(nrow = length(df[,1]), ncol = length(df[,1]))

  #we look each couple i j and compute the overlap index for this couple and save the value in the matrix[i,j]
  for(i in 1:(length(df[,1])-1)){
    for(j in (i+1):length(df[,1])){
      #recuperation of the elements rows i and j
      ens1 <- df[i,]
      ens2 <- df[j,]
      #compute the overlap index
      ind <- indice(ens1[1,2],ens1[1,3],ens2[1,2],ens2[1,3])
      #if(ind < 0){ind <- 0}
      if(ind==99){stop("The upper and lower bound are not in the correct column (sup1 < inf1 || inf2 > sup2).")}
      #if(ind <= 0){ind <- 0}
      #if(ind > 0){ind <- 1}
      #fill the matrix
      overlap_matrix[i,j] <- overlap_matrix[j,i] <- ind
    }
  }

  #fill the diag with 1
  diag(overlap_matrix) <- 1
  #[-1,1] to [0,2] (2 mean identical time range)
  overlap_matrix <- overlap_matrix + 1
  #inversion of values (high values now mean distant elements)
  overlap_matrix <- 2 - overlap_matrix
  #[0,2] to [0,1]
  overlap_matrix <- overlap_matrix / 2
  rownames(overlap_matrix) <- colnames(overlap_matrix) <- df[,1]
  return(D=overlap_matrix)
}

perioclust <- hclustcompro <- function(D1,D2,alpha = "EstimateAlphaForMe",k = NULL,title = "notitle",method = "ward.D2",suppl_plot = TRUE) {
  #===============================================================================
  # CAH CLASSIFICATION USING TWO DIST MATRIX
  #
  # Authors : A.COULON, L.BELLANGER and P.HUSI
  # Last update : 23 october 2018
  #
  # Arguments:
  # D1	      The fisrt distance matrix. Archeological context: Stratigraphic or timerange distance.
  # D2	      The second distance matrix. Archeological context: ceramic distance.
  # alpha	    The mixing parameter
  # k         The number of clusters
  # title	    The title to print on the dendrogram. (optionnal)
  # method    The method to use in hclust (see doc)
  # suppl_plot Logical for plot the WSS and average sil plot.
  #
  #===============================================================================

  sort = TRUE
  if(is.null(labels(D1))){
    sort <- FALSE
  }
  if(is.null(labels(D2))){
    sort <- FALSE
  }

  METHODS <- c("ward.D","single","complete","average","mcquitty","median","centroid","ward.D2")
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
  cont <- NULL
  Number_of_groups <- c()
  y <- c()
  yend <- c()
  xend <- c()
  if(!"dist" %in% class(D2)){
    #matrix ?
    if(!is.matrix(D2)){
      if(length(D2[1,]) == 3){
        message("generation of a temporal dissimilarities matrix.")
        D2 <- overlap(D2)
      }else{
        message("generation of a stratigraphic dissimilarities matrix.")
        D2 <- adjacency(D2)
        #diag = 0 ?
        if(sum(diag(D2)) != 0){
          stop("D2: The diagonal must be equal to 0.")
        }
      }
    }

    #symetric ?
    if(!isSymmetric(D2)){
      stop("D2 is not symmetric.")
    }
    if(!dim(D2)[1] == dim(D2)[2]){stop("D2 not a square matrix.")}
  }

  if(!"dist" %in% class(D1)) {
    #matrix ?
    if(!is.matrix(D1)){
      if(is.data.frame(D1)){cont <- D1;EPPM <- TRUE}else{EPPM <- FALSE}
      message("generation of a ceramic distances matrix.")
      D1 <- as.matrix(CAdist(D1,nPC="max",graph=FALSE))
    }
    #diag = 0 ?
    if(sum(diag(D1)) != 0){
      stop("D1: The diagonal must be equal to 0.")
    }
    #symetric ?
    if(!isSymmetric(D1)){
      stop("D1 is not symmetric.")
    }
    if(!dim(D1)[1] == dim(D1)[2]){stop("D1 not a square matrix.")}
  }

  D1 <- as.matrix(D1)
  max <- dim(D1)[1]
  D2 <- as.matrix(D2)

  if(FALSE %in% (sort(rownames(D1)) == sort(rownames(D2)))){
    warning("The names of D1 and D2 are not the same.")
    sort <- FALSE
  }
  if(FALSE %in% (sort(colnames(D1)) == sort(colnames(D2)))){
    warning("The names of D1 and D2 are not the same.")
    sort <- FALSE
  }
  if(FALSE %in% (sort(rownames(D1)) == sort(colnames(D2)))){
    warning("The names of D1 and D2 are not the same.")
    sort <- FALSE
  }

  if(!length(D1) == length(D2)){stop("D1 and D2 have not the same size.")}
  #dist (positive matrix) ?
  if(min(D1) < 0){stop("D1 is not positive.")}
  if(min(D2) < 0){stop("D2 is not positive.")}

  #no alpha
  if(alpha == "EstimateAlphaForMe"){
    sa <- hclustcompro_select_alpha(D1, D2, method = method, resampling = FALSE)
    alpha <- sa$alpha[1]
    sa <- list(alpha.plot = NULL,conf = NULL)
  }else{
    sa <- list(alpha.plot = NULL,conf = NULL)
  }
  #alpha 0-1 ?
  if(alpha > 1 || alpha < 0){stop("alpha must be between 0 and 1.")}
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
  if (sort) {
    D1.Tmp <- D1[sort(labels(D1)[[1]]),sort(labels(D1)[[1]])]
    D2.Tmp <- D2[sort(labels(D2)[[1]]),sort(labels(D2)[[1]])]
    D1 <- as.dist(D1.Tmp)
    D2 <- as.dist(D2.Tmp)
  }

  #mixt matrix
  D_alpha <- (alpha)*D1 + (1-alpha)*D2
  message("\nThe mixing parameter is ",alpha,". The mixing matrix is made of ",arrondi((alpha)*100,2),"% of D1 and ",
      arrondi((1-alpha)*100,2),"% of D2.\n\n")
  D_alpha <- as.dist(D_alpha)

  #CAH
  tree <- SPARTAAS::hclust(d=D1,method=method,d2=D2,alpha=alpha)

  if(title == "notitle"){
    title <- paste("Dendrogram with alpha =",alpha)
  }
  subtitle <- paste("The mixing matrix is made of ",arrondi((alpha)*100,2),"% of D1 and ",
                    arrondi((1-alpha)*100,2),"% of D2.")

  D <- as.matrix(D_alpha)
  n <- length(D[,1])

  #wss sil value
  if(suppl_plot){
    if(length(D[,1]) > 20){n <- 20}
    wss <- c()
    sil <- c()
    for ( i in 2:(n-1)){
      cutree <- stats::cutree(tree,i)
      res <- fpc::cluster.stats(D,cutree)
      wss <- c(wss,res$within.cluster.ss)
      sil <- c(sil,res$avg.silwidth)
    }
    #wss
    Within_Sum_of_Square <- rbind(seq(2,(n-1)),wss)
    color_wss <- elbow_finder(seq(2,(n-1)),wss)
    dataplot <- t(Within_Sum_of_Square)
    dataplot <- as.data.frame(dataplot)
    colnames(dataplot) <- c("Number_of_groups","Within_Sum_of_Square")
    q <- ggplot2::ggplot(dataplot,aes(x=Number_of_groups, y=Within_Sum_of_Square)) +
      ggplot2::geom_point(aes(color=color_wss)) +
      ggplot2::scale_colour_gradientn(colours=c("brown1","chartreuse3")) +
      ggplot2::geom_line(linetype = 3) +
      ggplot2::scale_x_continuous(breaks = dataplot$Number_of_groups ,
                                  labels = dataplot$Number_of_groups ) +
      ylim(0,max(wss)+1) +
      labs(title = "Within Sum of Square" , subtitle = "Partition evaluation") +
      ggplot2::theme(legend.position="none")

    #silhouette
    Average_Sil_Width <- rbind(seq(2,(n-1)),sil)
    dataplot <- t(Average_Sil_Width)
    dataplot <- as.data.frame(dataplot)
    dataplot <- na.omit(dataplot)
    colnames(dataplot) <- c("Number_of_groups","Average_Sil_Width")
    p <- ggplot2::ggplot(dataplot,aes(x = Number_of_groups, y = Average_Sil_Width)) +
      ggplot2::geom_point(aes(color=Average_Sil_Width)) +
      ggplot2::scale_colour_gradientn(colours=c("brown1","chartreuse3")) +
      ggplot2::geom_line(linetype = 3) + ylim(0,1) +
      ggplot2::scale_x_continuous(breaks = dataplot$Number_of_groups,
                                  labels = dataplot$Number_of_groups) +
      labs(title = "Average Silhouette Width", subtitle = "Partition evaluation") +
      ggplot2::theme(legend.position="none")

    # Nouvelle page
    grid::grid.newpage()
    # Creer la mise en page : nrow = 2, ncol = 2
    grid::pushViewport(viewport(layout = grid.layout(2, 2)))
    # Une fonction pour definir une region dans la mise en page
    define_region <- function(row, col){
      grid::viewport(layout.pos.row = row, layout.pos.col = col)
    }

    # Arranger les graphiques dans la fenetre
    hc <- tree
    hcdata <- ggdendro::dendro_data(tree, type="rectangle")
    dendrogram <- ggplot2::ggplot() + ggplot2::ggtitle(title) + ylim(-0.2, max(tree$height)) + ylab("height") + xlab("") +
      ggplot2::geom_segment(data=segment(hcdata), aes(x=x, y=y, xend=xend, yend=yend)) +
      ggplot2::geom_text(data=label(hcdata), aes(x=x, y=y,label=label,hjust=1,angle=90,fontface="bold",cex=1.5),
                         size=3) +
      ggplot2::theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank())

    print(dendrogram, vp = define_region(1:2, 1))
    print(p, vp = define_region(1, 2))
    print(q, vp = define_region(2, 2))

    if(is.null(k)){
      if(interactive()){
        OK <- TRUE
        while(OK){
          k <- as.integer(readline("Select the number of clusters:"))
          if(is.na(k) || k < 1 || k > max){
            message("must be a integer between 1 and ",max)
          }else{
            OK <- FALSE
          }
        }
      }else{
        stop("No interactive session. Use parameter k to set clusters.")
      }
    }
  }else{
    if(is.null(k)){
      plot(tree,h = -1,main = title,sub = subtitle,xlab = paste0("D_alpha, method: \"",method,"\""))
      if(interactive()){
        OK <- TRUE
        while(OK){
          k <- as.integer(readline("Select the number of clusters:"))
          if(is.na(k) || k < 1 || k > max){
            message("must be a integer between 1 and ",max)
          }else{
            OK <- FALSE
          }
        }
      }else{
        stop("No interactive session. Use parameter k to set clusters.")
      }
    }
  }


  plot(tree,h = -1,main = title,sub = subtitle,xlab = paste0("D_alpha, method: \"",method,"\""))
  stats::rect.hclust(tree, k = k, border = colorspace::rainbow_hcl(k, c = 80, l = 60, start = 0, end = 360 * (k-1) / k))
  cluster <- stats::cutree(tree, k = k)
  height <- tree$height[length(tree$height) - k + 2]
  order <- unique(cluster[tree$order])
  for(i in 1:length(order)){
    if(i != order[i] && i < order[i]){
      cluster[cluster == i] <- 0
      cluster[cluster == order[i]] <- i
      cluster[cluster == 0] <- order[i]
      order <- unique(cluster[tree$order])
    }
  }
  col <- colorspace::rainbow_hcl(length(unique(cluster)),
                     c = 80, l = 60, start = 0,
                     end = 360*(length(unique(cluster))-1)/(length(unique(cluster))))
  for(i in 1:length(unique(cluster))){
    x <- mean(which(cluster[tree$order] == sort(unique(cluster))[i]))

    if(length(unique(cluster)) > 338){
      symb_letter <- sort(unique(cluster))[i]
    }else{
      symb_letter <- Int2Alpha(sort(unique(cluster))[i])
    }
    mtext(symb_letter, side = 1, line = .5, at = x, col = col[i], cex = 1.5)
  }
  cutree <- grDevices::recordPlot()

  return(
    structure(
      list(
        D1 = D1,
        D2 = D2,
        D_alpha = D_alpha,
        alpha = alpha,
        alpha.plot = sa$alpha.plot,
        conf_alpha = sa$conf,
        tree = tree,
        cluster = cluster,
        cutree = cutree,
        call = match.call(),
        cont = cont,
        wss=q,
        avg_sil=p
      ),
      class = c("hclustcompro_cl","list")
    )
  )
}

hclustcompro_subdivide <- function(hclustcompro_cl,cluster,nb_class){

  #test
  if("hclustcompro_cl" != class(hclustcompro_cl)[1]){stop("Not a hclustcompro object.")}
  if(length(cluster) != length(nb_class)){stop("cluster and nb_cluster must have the same length.")}
  title <- paste("Cluster Dendrogram with alpha =",hclustcompro_cl$alpha)
  subtitle <- paste("The mixing matrix is made of ",arrondi((hclustcompro_cl$alpha)*100,2),"% of D1 and ",
                    arrondi((1-hclustcompro_cl$alpha)*100,2),"% of D2.")
  xlab <- paste("D_alpha, method: \"ward.D2\"")
  if(!is.null(hclustcompro_cl$oldcluster)){
    hclustcompro_cl$cluster <- hclustcompro_cl$oldcluster
  }
  k <- length(unique(hclustcompro_cl$cluster))
  color <- colorspace::rainbow_hcl(k, c = 80, l = 60, start = 0,
                       end = 360*(k-1)/(k))
  plot(hclustcompro_cl$tree, h=-1, main = title, sub = subtitle, xlab=xlab)
  stats::rect.hclust(hclustcompro_cl$tree,k = k,border = color)
  clusters <- c()

  cmpt <- 0
  for(i in 1:length(cluster)){
    where <- cluster[i]
    sk <- nb_class[i]
    number <- length(hclustcompro_cl$cluster[hclustcompro_cl$cluster == where])
    if(!sk %in% seq(2,(number-1)) ){stop("You try to create more sub cluster than possible.")}
    #col2 <- rainbow_hcl((k+sk-1), c = 80, l = 60, start = 0, end = 360*((k+sk-1)-1)/(k+sk-1))
    #add
    col <- colorspace::rainbow_hcl(k, c = 80, l = 60, start = 0,
                                end = 360*(k-1)/(k))
    col <- col[where]
    color <- append(color, rep(col,sk-1), after = where+cmpt)
    cmpt <- cmpt + sk-1
    #end
    if(length(hclustcompro_cl$cluster[hclustcompro_cl$cluster == where]) < sk){
      stop("There is not enought element in this class to create ",sk," sub-classes.")
    }
    subdiv <- sub_div(hclustcompro_cl$tree,hclustcompro_cl$cluster,where,sk)
    step <- subdiv$step
    clusters <- fusion(clusters,subdiv$cluster)

    #for color
    nbcluster <- stats::cutree(hclustcompro_cl$tree,h=hclustcompro_cl$tree$height[step])
    order <- unique(nbcluster[hclustcompro_cl$tree$order])
    for(i in 1:length(order)){
      if(i != order[i] && i < order[i]){
        nbcluster[nbcluster == i] <- 0
        nbcluster[nbcluster == order[i]] <- i
        nbcluster[nbcluster == 0] <- order[i]
        order <- unique(nbcluster[hclustcompro_cl$tree$order])
      }
    }

    nb <- which(hclustcompro_cl$cluster == where)
    nb <- nbcluster[nb]
    which <- sort(unique(nb))
    stats::rect.hclust(hclustcompro_cl$tree,h = hclustcompro_cl$tree$height[step+1],which = which,
                border = col)  #col2[seq(where,where+sk-1)])
  }
  #col2 <- rainbow_hcl(length(unique(clusters)), c = 80, l = 60, start = 0,
  #                   end = 360*(length(unique(clusters))-1)/(length(unique(clusters))))
  for(i in 1:length(unique(clusters))){
    col <- color[i]
    x <- mean(which(clusters[hclustcompro_cl$tree$order] == sort(unique(clusters))[i]))
    # col2[i], cex=1.5)
    graphics::mtext(Int2Alpha(sort(unique(clusters))[i]), side = 1, line = .5, at = x, col = col, cex=1.5)
  }
  cutree <- grDevices::recordPlot()
  return(
    structure(
      list(
        D1 = hclustcompro_cl$D1,
        D2 = hclustcompro_cl$D2,
        D_alpha = hclustcompro_cl$D_alpha,
        alpha = hclustcompro_cl$alpha,
        tree = hclustcompro_cl$tree,
        cluster = clusters,
        oldcluster = hclustcompro_cl$cluster,
        cutree = cutree,
        call = match.call(),
        cont = hclustcompro_cl$cont
      ),
      class = c("hclustcompro_cl","list")
    )
  )
}

seriograph <- function(cont, order = NULL, insert = NULL, show = "both", permute = TRUE){

  show <- tolower(show)
  SHOWS <- c("both", "frequency", "eppm")
  i.show <- pmatch(show, SHOWS)
  if (is.na(i.show))
    stop("invalid show method", paste("", show))
  if (i.show == -1)
    stop("ambiguous show method", paste("", show))

  show <- SHOWS[i.show]
  if(show == "eppm"){
    show <- "EPPM"
  }

  if(class(cont)[1] == "hclustcompro_cl"){
    if(is.null(cont$cont)){
      stop("We don't find contingency table in hclustcompro_cl object.")
    }else{
      hclustcompro_cl <- cont
      cont <- hclustcompro_cl$cont
    }
    remove <- Int2Alpha(sort(unique(hclustcompro_cl$cluster)))
    if(is.null(order)){
      order <- Int2Alpha(sort(unique(hclustcompro_cl$cluster)))
    }
    remove <- remove[!remove %in% order]
    if(length(remove) >= 1){
      remove <- Alpha2Int(remove)
    }
    Cluster <- c(0)

    split <- Alpha2Int(order)
    seq <- hclustcompro_cl$cluster
    if(length(remove) >= 1){
      removeseq <- which(seq %in% remove)
      for(i in 1:length(removeseq)){
        seq[removeseq[i]] <- 900 + i
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
        seq[hclustcompro_cl$cluster == elt] <- as.numeric(paste0(c,".",d))
        d <- d + 1
      }else{
        c <- c + 1
        seq[hclustcompro_cl$cluster == elt] <- c
        d <- 1
      }
      pelt <- elt
    }
    message(blue("Distribution of elements in the different clusters:\n"))

    if(length(remove) >= 1){
      removeseq <- which(seq >= 900)
      seq2 <- seq[-removeseq]
      print(seq2)
    }else{
      print(seq)
    }

    cont <- cont[sort(labels(cont)[[1]]),]
    cont2 <- dplyr::mutate(cont, Cluster = factor(seq))
    cont2 <- dplyr::group_by(cont2, Cluster)
    cont2 <- dplyr::summarise_all(cont2, list(sum=sum))
    cont2 <- as.data.frame(cont2)
    rownames(cont2) <- cont2[,1]
    cont2 <- cont2[,-1]
    if(length(remove) >= 1){
      remove <- which(rownames(cont2) >= 900)
      cont2 <- cont2[-remove,]
    }
    cont2 <- cont2[order(as.numeric(rownames(cont2)),decreasing = T),]
    rownames(cont2) <- order

    if(!is.null(insert)){
      #insert hiatus
      ignore <- FALSE
      if(is.list(insert)){
        if(!is.numeric(insert[[1]])){
          if(!is.numeric(insert[[2]])){
            warning("Insert must contain integers corresponding to the position where to insert Hiatus.")
            ignore <- TRUE
          }else{
            tmp <- insert[[1]]
            insert[[1]] <- insert[[2]]
            insert[[2]] <- tmp
          }
        }
        if(!ignore){
          if(FALSE %in% (insert[[1]] == trunc(insert[[1]])) ){
            warning("Insert must contain integers corresponding to the position where to insert Hiatus.")
            ignore <- TRUE
          }else{
            if(!is.character(insert[[2]])){
              warning("Insert must contain labels as a character vector in the list.")
              ignore <- TRUE
            }
          }
        }
        if(!ignore){
          insert[[1]] <- length(cont2[,1]) - insert[[1]]
          label <- make.names(insert[[2]],unique=TRUE)
          for (i in length(insert[[1]]):1){
            after <- sort(insert[[1]])[i]
            if(after < length(cont2[,1]) && after >= 1 && after == trunc(after)){
              cont2 <- rbind(cont2[1:after,], label = rep(0,length(cont2[1,])),
                             cont2[(after+1):length(cont2[,1]),])
              row.names(cont2)[after+1] <- label[i]
            }else{
              warning("Some positions are outside the range. Ex: ",(length(cont2[,1]) - after))
            }
          }
        }
      }else{
        if(is.vector(insert)){
          if(!is.numeric(insert)){
            warning("Insert must contain integers corresponding to the position where to insert Hiatus.")
            ignore <- TRUE
          }
          if(!ignore){
            if(FALSE %in% (insert[[1]] == trunc(insert[[1]])) ){
              warning("Insert must contain integers corresponding to the position where to insert Hiatus.")
              ignore <- TRUE
            }
          }
          if(!ignore){
            insert <- length(cont2[,1]) - insert
            label <- "Hiatus"
            insert <- sort(insert)
            for (i in length(insert):1){
              if(insert[i] < length(cont2[,1]) && insert[i] >= 1 && insert[i] == trunc(insert[i])){
                cont2 <- rbind(cont2[1:insert[i],], label = rep(0,length(cont2[1,])),
                               cont2[(insert[i]+1):length(cont2[,1]),])
                if(length(insert)== 1){
                  row.names(cont2)[insert[i]+1] <- label
                }else{
                  row.names(cont2)[insert[i]+1] <- paste0(label, ".", (length(insert) + 1 - i) )
                }
              }else{
                warning("Some positions are outside the range. Ex: ",(length(cont2[,1]) - insert[i]))
              }
            }
          }
        }else{
          stop("insert must be a list of two vector (position after which to insert and label)
               or a vector (position)")
        }
      }
    }
    res <- plot_EPPM(cont2,show,permute)
    return(
      structure(
        list(
          seriograph = res$seriograph,
          contingency = res$contingency,
          frequency = res$frequency,
          ecart = res$ecart
        ),
        class = "seriograph"
      )
    )

  }else{
    #test data.frame
    if(!is.data.frame(cont)){
      cont <- as.data.frame(cont)
      if(!is.data.frame(cont)){stop("cont must be a data frame.")}
      if(!FALSE %in% (cont == trunc(cont))){stop("A contingency table must contains integers values.")}
    }
    cont2 <- cont

    #order
    if(!is.null(order)){
      new_order <- rep(NA,length(order))

      for(i in 1:length(order)){
        new_order[i] <- which(labels(cont2)[[1]] == order[i])
      }
      cont2 <- cont2[new_order,]

      #for (i in 1:length(order)){
      #  order[i] <- which(rownames(cont2) == order[i])
      #}
      #cont2 <- cont2[order(as.numeric(order),decreasing = F),]
    }

    if(!is.null(insert)){
      #insert hiatus
      ignore <- FALSE
      if(is.list(insert)){
        if(!is.numeric(insert[[1]])){
          if(!is.numeric(insert[[2]])){
            warning("Insert must contain integers corresponding to the position where to insert Hiatus.")
            ignore <- TRUE
          }else{
            tmp <- insert[[1]]
            insert[[1]] <- insert[[2]]
            insert[[2]] <- tmp
          }
        }
        if(!ignore){
          if(FALSE %in% (insert[[1]] == trunc(insert[[1]])) ){
            warning("Insert must contain integers corresponding to the position where to insert Hiatus.")
            ignore <- TRUE
          }else{
            if(!is.character(insert[[2]])){
              warning("Insert must contain labels as a character vector in the list.")
              ignore <- TRUE
            }
          }
        }
        if(!ignore){
          insert[[1]] <- length(cont2[,1]) - insert[[1]]
          label <- make.names(insert[[2]],unique=TRUE)
          for (i in length(insert[[1]]):1){
            after <- sort(insert[[1]])[i]
            if(after < length(cont2[,1]) && after >= 1 && after == trunc(after)){
              cont2 <- rbind(cont2[1:after,], label = rep(0,length(cont2[1,])),
                             cont2[(after+1):length(cont2[,1]),])
              row.names(cont2)[after+1] <- label[i]
            }else{
              warning("Some positions are outside the range. Ex: ",(length(cont2[,1]) - after))
            }
          }
        }
      }else{
        if(is.vector(insert)){
          if(!is.numeric(insert)){
            warning("Insert must contain integers corresponding to the position where to insert Hiatus.")
            ignore <- TRUE
          }
          if(!ignore){
            if(FALSE %in% (insert[[1]] == trunc(insert[[1]])) ){
              warning("Insert must contain integers corresponding to the position where to insert Hiatus.")
              ignore <- TRUE
            }
          }
          if(!ignore){
            insert <- length(cont2[,1]) - insert
            label <- "Hiatus"
            insert <- sort(insert)
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
                warning("Some positions are outside the range. Ex: ",(length(cont2[,1]) - insert[i]))
              }
            }
          }
        }else{
          stop("insert must be a list of two vector (position after which to insert and label)
               or a vector (position)")
        }
      }
    }
    res <- plot_EPPM(cont2,show,permute)
    return(
      structure(
        list(
          seriograph = res$seriograph,
          contingency = res$contingency,
          frequency = res$frequency,
          ecart = res$ecart
        ),
        class = "seriograph"
      )
    )
  }
}


print.hclustcompro_cl <- function(x, ...){

  cat("     hclustcompro | hclustcompro_cl

Alpha:",x$alpha,"

call: ")
  print(x$call)
  cat("
value:
      .. $D1          The first dissimilarities matrix
      .. $D2          The second dissimilarities matrix
      .. $D_alpha     The mixing dissimilarities matrix
      .. $alpha       The mixing parameter
      .. $tree        An object of class hclust (see ?hclust)

If you ran hclustcompro like this: classif <- hclustcompro(...)
select another number of cluster on the plot:
> plot(classif$tree, h = -1)
> rect.hclust(classif$tree, k = k) #with k the number of cluster

documentation:
> ?hclustcompro
      ")
}

hclustcompro_detail_resampling <- function(D1, D2 = NULL, acc = 2, method = "ward.D2", iter = 5){

  #===============================================================================
  # CAH CLASSIFICATION USING TWO DIST MATRIX
  #
  # Authors : A.COULON, L.BELLANGER and P.HUSI
  # Last update : 23 october 2018
  #
  # Arguments:
  # D1	                  The fisrt distance matrix. Archeological context: stratigraphic distance.
  # D2	                  The second distance matrix. Archeological context: ceramic distance.
  # acc	                  The accuracy to round.
  # method                The method to use in hclust (see doc)
  # iter                   The number of replication for each element
  #
  #===============================================================================
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
      stop("We don't find dissmilarities matrices in hclustcompro_cl object.")
    }else{
      D2 <- D1$D2
      D1 <- D1$D1
    }
  }else{
    if(is.null(D1) || is.null(D2)){
      stop("We don't find dissmilarities in D1 and D2 arguments.")
    }
    #test
    if(!"dist" %in% class(D1)){
      #matrix ?
      if(!is.matrix(D1)){
        message("generation of a ceramic distances matrix.")
        D1 <- as.matrix(CAdist(D1,nPC="max",graph=FALSE))
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
    if(!"dist" %in% class(D2)){
      #matrix ?
      if(!is.matrix(D2)){
        if(length(D2[1,]) == 3){
          message("generation of a temporal dissimilarities matrix.")
          D2 <- overlap(D2)
        }else{
          message("generation of a stratigraphic dissimilarities matrix.")
          D2 <- adjacency(D2)
          #diag = 0 ?
          if(sum(diag(D2)) != 0){
            stop("D2: The diagonal must be equal to 0.")
          }
        }
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

  save_lines <- NULL
  alpha_seq <- seq(0,1,0.01)

  D1save <- as.matrix(D1)
  D2save <- as.matrix(D2)

  #echantillonnage
  if(iter == "max"){iter <- length(D1save[1,])}
  if(!is.numeric(iter)){stop("iter must be numeric.")}
  if(iter < 1){iter <- 1; message("iter must be higher than 1.")
  warning("iter: must be between 1 and the number of elements.")}
  if(iter > length(D1save[1,])-1){
    if(iter > length(D1save[1,])){
      message("iter can't be higher than the number of element: ",length(D1save[1,]),".")
      warning("iter must be between 1 and the number of elements.")
    }
    iter <- length(D1save[1,])-1
  }

  cat("Resampling process:\n")
  #creation barre de progression(pb)
  pb <- utils::txtProgressBar(min = 0, max = length(D1save[1,]), style = 3)
  df <- NULL
  for (i in 1:length(D1save[1,])){
    for(j in sample((1:(length(D1save[,1])))[-i],iter,replace=FALSE)){
      if(iter == 1 && i == j){
        if(j == 1){j <- j + 1}else{j <- j - 1}
      }
      list <- generateClone(D1save,D2save,i,j)
      D1 <- list[[1]]
      D2 <- list[[2]]
      dist.dend <- c()
      for(alpha in alpha_seq){
        Mdist <- (alpha) * D1 + (1-alpha) * D2
        mixt.dist <- as.dist(Mdist)
        tree <- stats::hclust(mixt.dist,method=method)
        new <- corCriterion(tree,D1,D2)
        dist.dend <- c(dist.dend,new)
      }
      if(is.null(df)){
        df <- data.frame(alpha = alpha_seq, height = dist.dend)
      }else{
        df <- cbind(df,dist.dend)
      }
    }
    mean <- unlist(lapply(as.data.frame(t(df[-1])),mean))
    if(is.null(save_lines)){
      save_lines <- cbind(df$alpha,mean)
    }else{
      save_lines <- cbind(save_lines,mean)
    }
    df <- NULL
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)

  colnames(save_lines) <- c("alpha",rownames(D1save))
  rownames(save_lines) <- c()

  p <- plotly::plot_ly()
  for (k in 2:length(colnames(save_lines))) {
    p <- plotly::add_lines(p,x=save_lines[,1],y=save_lines[,k],name=colnames(save_lines)[k])
  }
  p
}

print.selectAlpha_obj <- function(x, ...){

  cat("        select Alpha | selectAlpha_obj

Estimated alpha:",x$alpha,"

call: ")
  print(x$call)
  cat("
value:
      .. $alpha             The estimate value of the parameter alpha (min CorCrit_alpha).
      .. $alpha.plot        The CorCrit for all the possible alpha.
    if resampling = TRUE
      .. $sd                The standard deviation.
      .. $conf              The confidence interval of alpha.
      .. $boxplot           The boxplot of the distribution.
      .. $values            All the potential alpha values obtained from clones.

documentation:
> ?hclustcompro_select_alpha
")
}

summary.selectAlpha_obj <- function(object, ...){
  cat("Estimate Alpha:",object$alpha,"    ( sd:",object$sd,")
Confidence interval: [",object$conf[1],";",object$conf[2],"]

Estimation calculated with",NROW(object$values),"observations.

?hclustcompro_select_alpha for details.
")
}

print.seriograph <- function(x, ...){
  cat("")
}

hclustcompro_select_alpha <- function(D1,D2,acc = 2, resampling = TRUE, method = "ward.D2",iter = 5, suppl_plot = TRUE){

  #===============================================================================
  # CAH CLASSIFICATION USING TWO DIST MATRIX
  #
  # Authors : A.COULON, L.BELLANGER and P.HUSI
  # Last update : 08 avril 2019
  #
  # Arguments:
  # D1	                  The fisrt distance matrix. Archeological context: Stratigraphic distance.
  # D2	                  The second distance matrix. Archeological context: ceramic distance.
  # acc             	    The accuracy to round.
  # resampling            Logical if we have to calculate confidence interval with resampling strategy
  # method                The method to use in hclust (see doc)
  #
  #===============================================================================

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

  corCriterion_ <- function(alpha) {
    res <- c()
    for (elt in alpha){
      Mdist <- (elt)*D1 + (1-elt)*D2
      mixt.dist <- as.dist(Mdist)
      tree <- stats::hclust(mixt.dist,method=method)
      d2 <- stats::cophenetic(tree)
      res <- c(res,abs(cor(as.dist(D1),d2) - cor(as.dist(D2),d2)))
    }
    return(res)
  }

  cor_ <- function(alpha) {
    c1 <- c()
    c2 <- c()
    for (elt in alpha){
      Mdist <- (elt)*D1 + (1-elt)*D2
      mixt.dist <- as.dist(Mdist)
      tree <- stats::hclust(mixt.dist,method=method)
      d2 <- stats::cophenetic(tree)
      c1 <- c(c1, cor(as.dist(D1),d2))
      c2 <- c(c2, cor(as.dist(D2),d2))
    }
    return(list(c1,c2))
  }

  corCriterion_clone <- function(alpha,a,b) {
    res <- c()
    D1bis <- rbind(D1,clone=D1[a,])
    D1bis <- cbind(D1bis,clone=c(D1[a,],0))
    D2bis <- rbind(D2,clone=D2[b,])
    D2bis <- cbind(D2bis,clone=c(D2[b,],0))
    for (elt in alpha){
      Mdist <- (elt)*D1bis + (1-elt)*D2bis
      mixt.dist <- as.dist(Mdist)
      tree <- stats::hclust(mixt.dist,method=method)
      d2 <- stats::cophenetic(tree)
      res <- c(res,abs(cor(as.dist(D1bis),d2) - cor(as.dist(D2bis),d2)))
    }
    return(res)
  }

  if(!is.numeric(acc)) {stop("acc must be numeric!")}
  if(acc < 1) {
    acc <- 1;
    message("The minimal value allowed for acc is 1.")
  }
  seq_acc=0.001
  if(acc==1){seq_acc=0.1}
  if(acc==2){seq_acc=0.01}
  if(acc==3){seq_acc=0.001}
  #test
  if(!"dist" %in% class(D1)){
    #matrix ?
    if(!is.matrix(D1)){
      message("generation of distances matrix.")
      D1 <- as.matrix(CAdist(D1,nPC="max",graph=FALSE))
    }
    #diag = 0 ?
    if(sum(diag(D1)) != 0){
      stop("D1: in distance matrix the diagonal is equal to 0.")
    }
    #symetric ?
    if(!isSymmetric(D1)){
      stop("D1 is not symmetric.")
    }
    if(!length(D1[,1]) == length(D1[1,])){stop("D1 not a square matrix.")}
  }
  if(!"dist" %in% class(D2)){
    #matrix ?
    if(!is.matrix(D2)){
      if(length(D2[1,]) == 3){
        message("generation of a temporal dissimilarities matrix.")
        D2 <- overlap(D2)
      }else{
        message("generation of a stratigraphic dissimilarities matrix.")
        D2 <- adjacency(D2)
        #diag = 0 ?
        if(sum(diag(D2)) != 0){
          stop("D2: The diagonal must be equal to 0.")
        }
      }
    }
    #symetric ?
    if(!isSymmetric(D2)){
      stop("D2 is not symmetric.")
    }
    if(!length(D2[,1]) == length(D2[1,])){stop("D2 not a square matrix.")}
  }

  D1 <- as.matrix(D1)
  D2 <- as.matrix(D2)
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

  #resampling check
  if(iter == "max"){iter <- length(D1[1,])}
  if(!is.numeric(iter)){stop("iter must be numeric.")}
  if(!is.logical(resampling)){stop("resampling must be TRUE or FALSE.")}
  if(iter < 1){iter <- 1; message("iter must be higher than 1.")
  warning("iter: must be between 1 and the number of elements.")}
  if(iter > length(D1[1,])-1){
    if(iter > length(D1[1,])){
      message("iter can't be higher than the number of element: ",length(D1[1,]),".")
      warning("iter: must be between 1 and the number of elements.")
    }
    iter <- length(D1[1,])-1
  }

  if(suppl_plot){
    #estime alpha
    corCrit <- c()
    for(i in seq(0,1,seq_acc)){corCrit <- c(corCrit,corCriterion_(i))}
    index <- which(corCrit == min(corCrit))
    alpha <- seq(0,1,seq_acc)[index]

    add_text_resamp <- ""
    if(resampling){
      #echantillonnage
      cat("\nResampling process:\n")
      pb <- utils::txtProgressBar(min = 0, max = length(D1[1,]), style = 3)
      set <- c()
      res <- c()
      for (i in 1:length(D1[1,])){
        for(j in sample((1:(length(D1[,1])))[-i],iter,replace=FALSE)){
          if(i!=j){
            tmp <- stats::optimize(corCriterion_clone,lower=0,upper=1,a=i,b=j)
            set <- c(set,tmp$minimum)
          }
        }
        utils::setTxtProgressBar(pb, i)
      }

      res <- set
      sd <- arrondi(sd(res),acc)
      box <- quantile(res,c(0,0.025,0.25,0.5,0.75,0.975,1))
      conf <- c(arrondi(box[2],acc),arrondi(box[6],acc))
      values <- arrondi(res,acc)

      BoxPlotValue <- boxplot(
        res,
        horizontal = TRUE,
        range = 0,
        main="boxplot: Distribution of all the potentials values of alpha.",
        xlab=expression(alpha),
        xaxt="n"
      )
      axis(
        side = 1,
        at = c(box[1],box[3],box[4],box[5],box[7]),
        labels = arrondi(c(box[1],box[3],box[4],box[5],box[7]),acc)
      )
      points(arrondi(mean(res),acc),1,col="red")
      text(arrondi(mean(res),acc),1.3,
           bquote(hat(alpha)^"*" == .(arrondi(mean(res),acc))))
      boxplot <- grDevices::recordPlot()

      add_text_resamp <- paste(" The IC95% is calculated with",iter,
                               "clones out of",length(D1[1,])-1,"available.")
    }


    plot(
      seq(0,1,seq_acc),corCriterion_(seq(0,1,seq_acc)),type="l",
      xlab = expression(alpha),
      ylab = expression(CorCrit[alpha]),
      ylim=c(0,1),
      axes = TRUE,
      main = "Select alpha",
      sub = paste0("method: ",method,".",add_text_resamp),
      xaxp = c(0,1,5),
      mgp= c(2, 1,0),
      col = "slateblue2",
      xaxt= "n",
      lwd=1
    )
    correlation_curve <- cor_(seq(0,1,seq_acc))
    lines(seq(0,1,seq_acc),correlation_curve[[1]],lty=3,col="grey60")
    lines(seq(0,1,seq_acc),correlation_curve[[2]],lty=3,col="grey60")
    axis(
      side = 1,
      at = c(.1,.2,.3,.4,.5,.6,.7,.8,.9),
      labels = rep("",9),
      col = "grey"
    )
    axis(
      side = 1,
      at = arrondi(alpha,acc),
      labels = rep("",length(alpha)),
      col = "slateblue2"
    )
    if(resampling){
      axis(
        side = 1,
        at = c(0,arrondi(box[2],acc),arrondi(box[6],acc),1),
        labels = c(0,arrondi(box[2],acc),arrondi(box[6],acc),1)
      )
      abline(v = arrondi(box[6],acc),col = "gray50",lty=3,lwd=2)
      abline(v = arrondi(box[2],acc),col = "gray50",lty=3,lwd=2)
      close(pb)
    }else{
      axis(
        side = 1,
        at = c(0,1),
        labels = c(0,1)
      )
    }

    for(i in 1:length(alpha)){
      graphics::mtext(arrondi(alpha[i],acc),side = 1,line = .4,at = arrondi(alpha[i],acc),col = "slateblue2")
    }
    alpha.plot <- grDevices::recordPlot()

    if(resampling){
      return(structure(
        list(alpha = arrondi(alpha,acc),
             alpha.plot = alpha.plot,
             sd = sd,
             conf = conf,
             boxplot = boxplot,
             values = values
        ),
        class = c("selectAlpha_obj","list")
      )
      )
    }else{
      return(structure(
        list(alpha = arrondi(alpha,acc),
             alpha.plot = alpha.plot,
             more = "Need resampling = TRUE!"
        ),
        class = c("selectAlpha_obj","list")
      )
      )
    }
  }else{

    #estime alpha
    #alpha <- stats::optimize(corCriterion_,lower=0,upper=1)$minimum
    corCrit <- c()
    for(i in seq(0,1,seq_acc)){corCrit <- c(corCrit,corCriterion_(i))}
    index <- which(corCrit == min(corCrit))
    alpha <- seq(0,1,seq_acc)[index]

    if(resampling){
      pb <- utils::txtProgressBar(min = 0, max = length(D1[1,]), style = 3)
      set <- c()
      res <- c()
      for (i in 1:length(D1[1,])){
        for(j in sample((1:(length(D1[,1])))[-i],iter,replace=FALSE)){
          if(i!=j){
            tmp <- stats::optimize(corCriterion_clone,lower=0,upper=1,a=i,b=j)
            set <- c(set,tmp$minimum)
          }
        }
        utils::setTxtProgressBar(pb, i)
      }
      conf <- arrondi(quantile(set,c(0.025,0.975)),acc)
      close(pb)
      return(structure(
        list(alpha = arrondi(alpha,acc),
             conf = conf,
             more = "Need suppl_plot = FALSE!"
        ),
        class = c("selectAlpha_obj","list")
      )
      )
    }else{
      return(structure(
        list(alpha = arrondi(alpha,acc),
             more = "Need suppl_plot = FALSE!"
        ),
        class = c("selectAlpha_obj","list")
      )
      )
    }


  }


}

timerange <- function(data, cluster = NULL, add = NULL, density = NULL, color = NULL, reorder = FALSE){

  #test
  if(is.null(cluster)){cluster <- rep(1,dim(data)[1])}
  if(!is.data.frame(data)){stop("data must be a data.frame.")}
  if(!is.null(add)){
    if(!is.data.frame(add)){stop("add must be a data.frame (even if there is only one column).")}
  }
  if(!is.vector(cluster)){stop("cluster must be a vector.")}
  if(!is.null(density)){
    if(!is.vector(density)){stop("cluster must be a vector.")}
  }

  #call fonction overlap plot
  df <- cbind(data,cluster)
  res <- overlap_plot(df,add,density,color,reorder_color = reorder)

  return(structure(list(order = res$order, plot = res$plot, density = res$density, cluster = res$reorder_cluster), class = c("temp_obj", "list")))
}

histogram <- function(df, cluster = NULL, order = NULL, k = 4, reorder_color = FALSE){
  ################################################################################################
  # df : Borne inf, Borne sup, Observation (ici GT), Nombre d'Observation (ici NTI), Numéro de groupe
  # k : les k observations les plus représentees
  ################################################################################################

  Cluster <- c()
  reorder <- TRUE

  # df <- data_NTI[,1:3]
  # cluster <- data_NTI[,4]
  # order <- order_NTI
  # k = 4

  #|
  #| check inputs
  #|

  if(is.null(cluster)){cluster <- rep(1,dim(df)[1])}

  if(length(df[1,]) == 2){
    df <- data.frame(
      niv1 <- df[,1],
      niv2 <- df[,1],
      Observation_number <- df[,2]
    )
  }

  if(k < 1){k<-1}
  df <- cbind(df,cluster)
  names(df) = c("niv1", "niv2", "Observation_number", "Cluster")
  ncol <- length(unique(trunc(df$Cluster)))
  df$Cluster = as.factor(df$Cluster)

  if(is.null(order)){
    order <- unique(df$Cluster)
    reorder <- FALSE
  }

  ### on ordonne
  df_order = {}
  for(i in order){df_order = rbind(df_order, df[df$Cluster == i,])}
  df = df_order

  if(reorder){
    ## Re numerotation (ancien au recent)
    order_tmp <- as.numeric(order)
    Cluster_tmp <- as.numeric(as.character(df$Cluster))

    new <- c()
    count <- 1
    for (i in 1:length(order_tmp)){
      if(order_tmp[i] == trunc(order_tmp[i])){
        new <- c(new,count)
        count <- count + 1
      }
      if(order_tmp[i] != trunc(order_tmp[i])){
        t <- order_tmp[i]
        if(i != length(order_tmp)){
          tun <- order_tmp[i+1]
        }else{
          tun <- 0
        }
        if(trunc(t) == trunc(tun)){
          new <- c(new,count+order_tmp[i]-trunc(t))
        }else{
          new <- c(new,count+order_tmp[i]-trunc(t))
          count <- count + 1
        }
      }
    }
    new <- sort(new)
    Cluster_tmp <- Cluster_tmp * 10
    for(i in 1:length(new)){
      Cluster_tmp[Cluster_tmp == 10*order_tmp[i]] <- new[i]
    }
    df$Cluster <- as.factor(Cluster_tmp)
  }


  if(sort(unique(as.character(df$niv1) == as.character(df$niv2)))){
    #### On définit les couleurs

    col <- rainbow_hcl(ncol,
                       c = 80, l = 60, start = 0,
                       end = 360*(ncol-1)/ncol)

    plot <- list()
    for(i in unique(df$Cluster)){
      ## On recupere les donnees correspondant aux groupes i
      dg = df[df$Cluster == i,]

      ## On regroupe les données en fonction de Observations et on compte le nombre de Observation_number
      ## On calcule la proportion
      dg = dg %>% group_by(niv1) %>% summarise(Observation_number = sum(Observation_number),
                                               Cluster = first(Cluster)) %>% arrange(-Observation_number)
      dg$sum_Observations = sum(dg$Observation_number)
      dg$Proportion = dg$Observation_number / dg$sum_Observations

      data <- dg[order(dg$Proportion,decreasing=T),][1:k,]

      data$niv1 = factor(data$niv1, levels = data$niv1)
      col_index <- trunc(as.numeric(i))

      p1 <- plot_ly(data, x = data$niv1, y = data$Proportion, type = "bar",hoverinfo = 'text',
                    textposition = 'auto',showlegend = F,
                    hovertext = paste("</br>", data$niv1,
                                 "</br> Proportion :",arrondi(data$Proportion,3),
                                 "</br> Number of observation :",data$Observation_number,
                                 "</br> Cluster :", data$Cluster
                    ),
                    marker = list(
                      color = col[col_index]
                    )
      ) %>%
        layout(
          title = "",
          yaxis = list(title = "Proportion", showline = FALSE, dtick = .1,range= c(0,1)),
          xaxis = list(title = paste("Class",(as.numeric(as.character(data$Cluster)))),tickangle = 90)
        )

      plot[[i]] <- p1
    }

    plot2 <- subplot(plot,shareY = TRUE,shareX=TRUE,margin = 0.003)
    plot2
    return(structure(list(plot = plot2), class = c("prop_obj", "list")))

  }else{
    plot <- list()
    for(i in unique(df$Cluster)){
      #### On définit les couleurs
      col <- rainbow_hcl(length(unique(trunc(as.numeric(cluster)))),
                         c = 80, l = 60, start = 0,
                         end = 360*(length(unique(trunc(as.numeric(cluster))))-1) /
                           (length(unique(trunc(as.numeric(cluster))))))

      ## On recupere les donnees correspondant aux groupes i
      dcluster = df[df$Cluster == i,]
      sum_Observations = sum(dcluster$Observation_number)

      nb <- dcluster %>% group_by(niv1) %>% summarise(n = length(unique(niv2)),sum = sum(Observation_number))
      max <- length(nb$n)
      if(k > max){k_tmp <- max}else{k_tmp <- k}
      nb <- nb[order(nb$sum,decreasing=TRUE),][1:k_tmp,]

      clustplot <- list()
      if(k > 10){k <- 10}
      for(l in 1:k){
        if(l > k_tmp){
          p1 <- plot_ly(x= "",y = 0, type = "bar",hoverinfo = 'text',
                        showlegend = F,
                        marker = list(
                          color = col[1],
                          line = list(
                            color = 'black',
                            width = 1
                          )
                        )
          )
        }else{
          test <- dcluster[dcluster$niv1 == as.character(nb$niv1)[l],] %>% group_by(niv2) %>%
            summarise(nb=sum(Observation_number))
          test$nb <- test$nb / sum_Observations
          test$niv2 = factor(test$niv2, levels = test$niv2)
          test <- test[order(test$nb,decreasing=TRUE),]

          p1 <- plot_ly(test, x = as.character(nb$niv1)[l], y = test$nb, type = "bar", hoverinfo = 'text',
                        textposition = 'auto',showlegend = F,
                        hovertext = paste("</br>", as.character(nb$niv1)[l],
                                     " -- ", test$niv2,
                                     "</br> Proportion :",arrondi(test$nb,3),
                                     "</br> Number of observation :",arrondi(test$nb * sum_Observations,0),
                                     "</br> Cluster :", i
                        ),
                        marker = list(
                          color = col[as.numeric(i)],
                          line = list(
                            color = 'white',
                            width = 1
                          )
                        )
          )
          p1 <- p1 %>% layout(yaxis = list(title = 'Proportion'), barmode = 'stack',
                              title = "",
                              xaxis = list(title = "", tickangle = -45),
                              yaxis = list(title = "Proportion", showline = FALSE, dtick = .1,range= c(0,1)))

        }
        clustplot[[l]] <- p1
        plot[[i]] <- subplot(clustplot,shareY = TRUE,shareX=TRUE,margin = 0)
        }

    }
    if(k*length(unique(df$Cluster)) > 99){nr <- 2;margin <- c(0.004,0.004,.05,.05)}else{nr <- 1;margin <- 0.004}
    plot2 <- subplot(plot, shareY = TRUE, shareX = FALSE, margin = margin,nrows=nr)
    plot2
    return(structure(list(plot = plot2), class = c("prop_obj", "list")))
  }
}

#overload hclust
hclust <- function(d, method = "complete", members = NULL, d2 = NULL, alpha = NULL){
  if (!is.null(d2)) {
    if(!length(d) == length(d2)){stop("d and d2 have not the same size.")}
    if (is.null(alpha)) {
      sa <- hclustcompro_select_alpha(d, d2, method = method, resampling = FALSE)
      alpha <- sa$alpha[1]
    }
    alpha <- as.numeric(alpha)
    if(!(alpha>=0 & alpha<=1)){
      warning("alpha must be between 0 and 1.")
      sa <- hclustcompro_select_alpha(d, d2, method = method, resampling = FALSE)
      alpha <- sa$alpha[1]
    }
    #normalization
    if(max(d) != 0){
      d <- d / max(d)
    }
    if(max(d2) != 0){
      d2 <- d2 / max(d2)
    }
    d <- as.dist(alpha * d + (1 - alpha) * d2)
  }
  stats::hclust(d, method, members)
}
