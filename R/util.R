generateClone <- function(D1,D2,i,j){
  D1bis <- rbind(D1,clone=D1[i,])
  D1bis <- cbind(D1bis,clone=c(D1[i,],0))

  D2bis <- rbind(D2,clone=D2[j,])
  D2bis <- cbind(D2bis,clone=c(D2[j,],0))

  return(list(D1bis,D2bis))
}

corCriterion <- function(tree,D1,D2) {
  d2 <- cophenetic(tree)
  res <- abs(cor(as.dist(D1),d2) - cor(as.dist(D2),d2))
  #res2 <- list()
  #res2$un <- cor(as.dist(D1),d2)
  #res2$de <- cor(as.dist(D2),d2)
  #return(res2)
  return(res)
}

search <- function(x){
  return(which(x == min(x)))
}

fusion <- function(c1, c2) {
  c3 <- c()
  for (i in 1:length(c2)) {
    c3[i] <- max(c1[i], c2[i])
  }
  names(c3) <- names(c2)
  return(c3)
}

sub_div <- function(tree,cluster,where,sk) {
  ok <- TRUE
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
  return(list(cluster = cluster,step = step))
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
    if(!(sum(x[,i]) == arrondi(sum(x[,i])))){stop("The data frame must contains integers.")}
  }
  ens = 0
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
  data <- gather(data,key = "type", value = "frequency",-.data$ens,-.data$rowsum, factor_key = TRUE)
  data <- mutate(data,frequency = .data$frequency / .data$rowsum)
  #remove col rowsum
  data <- data[-2]
  data$frequency[is.na(data$frequency)] <- 0

  #EPPM
  ecart <- as.data.frame(EPM(x[-length(x[1,])]))
  save_ecart <- ecart
  ecart <- cbind(ecart,Weight=c(0))
  ecart[ecart < 0] <- 0
  ecart <- mutate(ecart,ens = labels)
  ecart <- gather(ecart,key = "type", value = "EPPM", -.data$ens, factor_key = TRUE)
  ecart$EPPM[ecart$type == "Weight"] <- 0

  # Join data and EPPM
  data <- inner_join(data,ecart,by = c("ens", "type"))
  data <- mutate(data,frequency = .data$frequency - .data$EPPM)
  data <- gather(data,key = "legend", value = "frequency", -.data$ens,
                 -.data$type, factor_key = TRUE)
  #if insert hiatus in df replace NA by 0 (NA are generating by divising by 0 (rowsum of a hiatus is 0))
  data[is.na(data)] <- 0
  #center (div /2 et copy for one part - and the other +)
  data <- rbind(data,data)
  data <- mutate(data,frequency = .data$frequency * c(rep(.5, nrow(data)/2), rep(-.5, nrow(data)/2)))

  #color for Weight
  tmp <- data$frequency[data$type == "Weight"]
  levels(data$legend) <- c(levels(data$legend),"1","2","3","4","5","6","7","8","9","10","Weight")
  data$legend[data$type == "Weight"] <- "Weight"

  breaks <- c("frequency","EPPM","Weight")

  #define the color
  if(show == "frequency"){
    color <- c("frequency"="#bebebe","EPPM"="#bebebe","Weight"="#7f7f7f")
    #remove EPPM from legend
    breaks <- c("frequency","Weight")
    default_col <- "#bebebe"
  }else{
    if(show == "EPPM"){
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
  contingency <- x[,1:(length(x[1,])-1)]
  contingency <- x[,-length(x[1,])]
  frequency <- contingency
  rowsum <- rowSums(contingency)
  for (k in 1:length(contingency[,1])){
    frequency[k,] <- frequency[k,] / rowsum[k]
  }

  cat("\nContingency table:\n")
  print(contingency)


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
      subtitle = "EPPM (Ecart Positif au Pourcentage Moyen) Positive Deviation from the Average Percentage",
      caption = "Warnings: The frequencies are calculated independently for each element (row).
You can see the relative number of observations in the Weight column"
    ) + geom_segment(aes(x = 0, y = -.04, xend = 0, yend = .04),linetype = "blank") +
    theme(legend.position = "bottom",
          panel.spacing = unit(.2, "lines"),
          strip.text.x = element_text(angle = 90),
          axis.text.x = element_text(hjust=1,margin = margin(t = -13)),
          axis.ticks.length=unit(.5, "cm"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.caption = element_text(face = "italic", hjust = 0)

    )
  plot(p)
  return(list(
      seriograph = p,
      contingency = contingency,
      frequency = frequency,
      ecart = save_ecart
    )
  )
}

overlap_plot <- function(df,add=NULL,density = NULL, color = NULL, reorder_color = FALSE){
  if(!is.null(color)){
    change_color <- TRUE
    save_color <- as.character(color)
  }else{
    change_color <- FALSE
  }

  #fonction pour faire des lignes
  line <- function(x0, x1, y0, y1, color = "#666666",width = 2){
    l = {}
    for(i in 1:length(x0)){
      l[[i]] = list(type = "line", y0 = y0[i], y1 = y1[i],
                    x0 = x0[i], x1 = x1[i], line = list(color = color, width = width))
    }
    return(l)
  }

  Lower_bound <- NULL
  Upper_bound <- NULL
  Observations <- NULL
  lower_median <- NULL

  Cluster <- df[,4]

  if(is.numeric(Cluster)){
    if(TRUE %in% (Cluster != trunc(Cluster))){

      class_subdivide <- which(Cluster != trunc(Cluster))
      nb_cust <- as.vector(table(trunc(unique(Cluster[class_subdivide]))))
      class_subdivide <- unique(trunc(Cluster[Cluster != trunc(Cluster)]))

      Cluster. <- Cluster
      Cluster <- trunc(Cluster)

      df[,4] <- trunc(Cluster)
      df <- cbind(df,Cluster.)

      names(df) = c("Observations", "Lower_bound", "Upper_bound", "Cluster", "Cluster.")
      df$Observations <- as.factor(df$Observations)
      df$Cluster = as.factor(df$Cluster)
      df$Cluster. = as.factor(df$Cluster.)
      if(!is.null(density)){
        df <- cbind(df,density)
      }
      if(change_color){
        df <- cbind(df,save_color)
      }
      if(!is.null(add)){
        if(dim(add)[2] == 1){
          n <- names(add)
          eval(parse(
            text=paste0("add <- data.frame(",n,
                        " = add[order(df$Cluster,df$Cluster.,df$Lower_bound,df$Upper_bound),])")))
        }else{
          add = as.data.frame(add)[order(df$Cluster,df$Cluster.,df$Lower_bound,df$Upper_bound),]
        }
      }
      df = df[order(df$Cluster,df$Cluster.,df$Lower_bound,df$Upper_bound),]
      #add = add[order(df$Cluster,df$Cluster.,df$Lower_bound,df$Upper_bound),]

      ### on permute les groupes en fonction de la médiane basse des groupes
      order = df  %>% group_by(Cluster) %>%
        summarise(median = median(Lower_bound)) %>%
        arrange(median)

      if(!is.null(add)){add <- cbind(add,Cluster=df$Cluster)}

      df_fin = {}
      add_fin = {}
      for(i in order$Cluster){
        add_order = {}
        df_order = {}
        new_tab = filter(df,Cluster == i)
        if(!is.null(add)){ new_add = filter(add,Cluster == i)}
        if(!new_tab$Cluster[1] == new_tab$Cluster.[1]){
          order_y = new_tab %>% group_by(Cluster.) %>%
            summarise(median = median(Lower_bound)) %>% arrange(median)
          for(i in order_y$Cluster.){
            if(!is.null(add)){ add_order = rbind(add_order,new_add[new_tab$Cluster. == i,]) }
            df_order = rbind(df_order,new_tab[new_tab$Cluster. == i,])
          }
        }else{
          if(!is.null(add)){add_order = new_add[new_tab$Cluster == i,]}
          df_order = new_tab[new_tab$Cluster == i,]
        }
        if(!is.null(add)){add_fin = rbind(add_fin,add_order)}
        df_fin = rbind(df_fin,df_order)
      }
      if(!is.null(add)){
        if(dim(add)[2] == 2){
          add <- data.frame(add_fin)
          index <- which(names(add) == "Cluster")
          print(n)
          eval(parse(text = paste0("add <- data.frame(",n," = add[,-index])")))
        }else{
          add = add_fin
          index <- which(names(add) == "Cluster")
          add <- data.frame(add[,-index])
        }
      }
      df = df_fin

      order = factor(unique(df_fin$Cluster.), levels = sort(unique(df_fin$Cluster.)))

      df$Cluster = factor(df$Cluster,levels = unique(df$Cluster))
      df$Cluster. = factor(df$Cluster.,levels = unique(df$Cluster.))

      ### on récupère les premières et dernières observation de chaque groupe
      intersection = df %>% group_by(Cluster) %>%
        summarise(First_observations = first(Observations), Last_observations = last(Observations))
      ### On récupère la position du changement de groupe
      intersection_place = which(df$Observations %in% intersection$First_observations)[-1] - 1.5


      ### on récupère les premières et dernières observation de chaque groupe et sous groupes
      intersection_SG = df %>% group_by(Cluster.) %>%
        summarise(First_observations = first(Observations), Last_observations = last(Observations))


      ### On calcule un l'intervalle de confiance pour les dates des groupes et sous groupes
      intervalle_inf = intervalle_sup = {}
      for(i in unique(df$Cluster.)){
        intervalle_inf[i] = arrondi(quantile(df[df$Cluster. == i,]$Lower_bound,.10),0)
        intervalle_sup[i] = arrondi(quantile(df[df$Cluster. == i,]$Upper_bound,.90),0)
      }

      intersection_SG$borne_inf = intervalle_inf
      intersection_SG$borne_sup = intervalle_sup


      ### On récupère la position du changement de groupe et sous groupe
      intersection_SG_place = which(df$Observations %in% intersection_SG$First_observations)[-1] - 1.5
      intersection_SG_x_min = c(-.5,intersection_SG_place)
      intersection_SG_x_max = c(intersection_SG_place, length(levels(df$Observation))-.5)

      intersection_SG$x_min = intersection_SG_x_min
      intersection_SG$x_max = intersection_SG_x_max

      #df$Cluster = factor(df$Cluster,levels = unique(df$Cluster))
      #df$Cluster. = factor(df$Cluster.,levels = unique(df$Cluster.))
      df$Observations = factor(df$Observations,levels = unique(df$Observations))


      if(reorder_color){
        ## Re numerotation (ancien au recent)
        order_tmp <- as.numeric(levels(df$Cluster.))
        Cluster_tmp <- as.numeric(as.character(df$Cluster.))

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
        df$Cluster. <- as.factor(Cluster_tmp)
        df$Cluster <- as.factor(trunc(Cluster_tmp))

        class_subdivide <- which(Cluster_tmp != trunc(Cluster_tmp))
        nb_cust <- as.vector(table(trunc(unique(Cluster_tmp[class_subdivide]))))
        class_subdivide <- unique(trunc(Cluster_tmp[Cluster_tmp != trunc(Cluster_tmp)]))

        ### On définit les couleurs
        if(change_color){
          color <- rainbow_hcl(length(unique(df$save_color)),
                               c = 80, l = 60, start = 0,
                               end = 360*(length(unique(df$save_color))-1)/(length(unique(df$save_color))))
        }else{
          color <- rainbow_hcl(length(unique(df$Cluster)),
                               c = 80, l = 60, start = 0,
                               end = 360*(length(unique(df$Cluster))-1)/(length(unique(df$Cluster))))
          new_color = {}
          for(i in 1:length(color)){
            if(i %in% class_subdivide){
              a = which(class_subdivide == i)
              color_i = rep(color[i],nb_cust[a])
            }else{
              color_i = color[i]
            }
            new_color = c(new_color,color_i)
          }
          color <- new_color
        }
      }else{
        Cluster_tmp <- as.numeric(as.character(df$Cluster.))

        class_subdivide <- which(Cluster_tmp != trunc(Cluster_tmp))
        nb_cust <- as.vector(table(trunc(unique(Cluster_tmp[class_subdivide]))))
        class_subdivide <- unique(trunc(Cluster_tmp[Cluster_tmp != trunc(Cluster_tmp)]))

        #### On définit les couleurs
        if(change_color){
          color <- rainbow_hcl(length(unique(df$save_color)),
                               c = 80, l = 60, start = 0,
                               end = 360*(length(unique(df$save_color))-1)/(length(unique(df$save_color))))
        }else{
          color <- rainbow_hcl(length(unique(df$Cluster)),
                               c = 80, l = 60, start = 0,
                               end = 360*(length(unique(df$Cluster))-1)/(length(unique(df$Cluster))))

          new_color = {}
          for(i in 1:length(color)){
            if(i %in% class_subdivide){
              a = which(class_subdivide == i)
              color_i = rep(color[i],nb_cust[a])
            }else{
              color_i = color[i]
            }
            new_color = c(new_color,color_i)
          }

          color <- new_color[as.numeric(order)]
        }
      }

      #graphique
      #fonction pour faire des lignes verticales
      vline <- function(x = 0, color = "#666666", dash = "", width = 2){
        l = {}
        for(i in 1:length(x)){
          l[[i]] = list(type = "line", y0 = 0, y1 = 1, yref = "paper",
                        x0 = x[i], x1 = x[i], line = list(color = color, width = width, dash = dash))
        }
        return(l)
      }
      a = vline(c(-.5,intersection_place,length(levels(df$Observation))-.5))
      b = vline(intersection_SG_place, dash = 'dot', width = 1)

      c = line(x0 = intersection_SG$x_min, x1 = intersection_SG$x_max, y0 = intersection_SG$borne_inf,
               y1 = intersection_SG$borne_inf)
      d = line(x0 = intersection_SG$x_min, x1 = intersection_SG$x_max, y0 = intersection_SG$borne_sup,
               y1 = intersection_SG$borne_sup)


      #on définit les axes
      ay <- list(
        title = "Years",
        ticks = "outside",
        dtick = 200)

      if(!is.null(add)){
        #add_text
        col_add_text <- NULL

        for(j in 1:dim(add)[1]){
          add_text <- NULL
          for(i in 1:dim(add)[2]){
            if(dim(add)[2] == 1){
              add_text <- paste0("</br> ",names(add)[i]," : ", as.character(add[,i][j]))
            }else{
              add_text <- paste0(add_text,"</br> ",names(add)[i]," : ", as.character(add[,i][j]))
            }
          }
          col_add_text <- c(col_add_text,add_text)
        }

        df$add_text <- col_add_text
      }else{
        df$add_text <- ""
      }

      if(change_color){
        save_color <- df$save_color
        name <- as.character(df$save_color)
      }else{
        save_color <- df$Cluster.
        name <- Int2Alpha(as.numeric(as.character(save_color)))
      }

      #on affiche le plot
      p = plot_ly(df) %>%
        add_segments(x = df$Observations, xend = df$Observations,
                     y = df$Lower_bound, yend = df$Lower_bound + ((df$Upper_bound - df$Lower_bound)/2),
                     showlegend = T, color = save_color, colors = color,
                     name=paste("Cluster", name),
                     legendgroup = paste("Cluster", name),
                     hoverinfo = 'text',
                     text = paste("</br> Observation :", df$Observations,
                                  "</br> Lower bound :", df$Lower_bound,
                                  "</br> Upper bound :", df$Upper_bound,
                                  "</br> Cluster :", Int2Alpha(as.numeric(as.character(df$Cluster.))),
                                  df$add_text
                     )) %>%
        add_segments(x = df$Observations, xend = df$Observations,
                     y =  df$Lower_bound + ((df$Upper_bound - df$Lower_bound)/2), yend = df$Upper_bound,
                     showlegend = F, color = save_color, colors = color,
                     legendgroup = paste("Cluster", name),
                     hoverinfo = 'text',
                     text = paste("</br> Observation :", df$Observations,
                                  "</br> Lower bound :", df$Lower_bound,
                                  "</br> Upper bound :", df$Upper_bound,
                                  "</br> Cluster :", Int2Alpha(as.numeric(as.character(df$Cluster.))),
                                  df$add_text
                     )) %>%
        add_markers(x = intersection_SG$First_observations, y = intersection_SG$borne_inf,
                    color = I('#666666'), marker = list(size = 0.01),
                    hoverinfo = "text", showlegend = F,
                    text = paste("</br> Interval : [", intersection_SG$borne_inf,
                                 ", ", intersection_SG$borne_sup, "]", sep = "")) %>%
        add_markers(x = intersection_SG$First_observations, y = intersection_SG$borne_sup,
                    color = I('#666666'), marker = list(size = 0.01),
                    hoverinfo = "text", showlegend = F,
                    text = paste("</br> Interval : [", intersection_SG$borne_inf,
                                 ", ", intersection_SG$borne_sup, "]", sep = "")) %>%
        add_markers(x = intersection_SG$Last_observations, y = intersection_SG$borne_inf,
                    color = I('#666666'), marker = list(size = 0.01),
                    hoverinfo = "text", showlegend = F,
                    text = paste("</br> Interval : [", intersection_SG$borne_inf,
                                 ", ", intersection_SG$borne_sup, "]", sep = "")) %>%
        add_markers(x = intersection_SG$Last_observations, y = intersection_SG$borne_sup,
                    color = I('#666666'), marker = list(size = 0.01),
                    hoverinfo = "text", showlegend = F,
                    text = paste("</br> Interval : [", intersection_SG$borne_inf,
                                 ", ", intersection_SG$borne_sup, "]", sep = ""))

      p <- p %>%
        layout(shapes = c(a,b,c,d),
               title = "Timerange",
               yaxis = list(title = "Years",zeroline = FALSE),
               xaxis = list(title = "Observations")
        )
      plot <- p
      if(!is.null(density)){
        save_dens <- df$density
        df$density <- df$density / sum(df$density)

        line <- function(x0, x1, y0, y1, color = "#666666",width = 2){
          l = {}
          for(i in 1:length(x0)){
            l[[i]] = list(type = "line", y0 = y0[i], y1 = y1[i], yref = "y2",
                          x0 = x0[i], x1 = x1[i], line = list(color = color, width = width))
          }
          return(l)
        }

        l <- line(x0 = 0:(length(df$density)-2),
                  y0 = df$density[1:(length(df$density)-1)],
                  x1 = 1:(length(df$density)-1),
                  y1 = df$density[2:length(df$density)],
                  color = "#E16A86",width = 2.5)

        q <- plot_ly(df) %>%
          add_markers(x = df$Observation, y = df$density, showlegend = F,
                      hoverinfo = 'text', color = I("#E16A86"), marker = list(size = 2),
                      hovertext = paste("</br> Observation :", df$Observations,
                                        "</br> Proportion :", arrondi(df$density,3),
                                        "</br> Variable :", arrondi(save_dens,3),
                                        df$add_text)
          ) %>%
          layout(shapes=c(a,b,l))


        plot <- subplot(
          p,q,nrows = 2, heights = c(0.7,0.3),shareX = TRUE, shareY = FALSE,titleY = TRUE
        ) %>%
          layout(title = "Periodization of observations by clusters",
                 yaxis = list(title="Years",zeroline = FALSE),
                 yaxis2 = list(title="Density"),
                 xaxis = list(title = "Observations")
          )
        density <- df$density
        names(density) <- df$Observation
      }

      df$Cluster <- df$Cluster.
      order = df %>% group_by(Cluster) %>%
        summarise(lower_median = median(Lower_bound),upper_median = median(Upper_bound)) %>%
        arrange(lower_median)
      #order$Cluster <- new
      return(structure(list(order = order, plot = plot, density = density, reorder_cluster = save_color),
                       class = c("temp_obj", "list")))
    }
  }

  names(df) = c("Observations", "Lower_bound", "Upper_bound", "Cluster")
  df$Cluster = as.factor(df$Cluster)
  if(!is.null(density)){
    df <- cbind(df,density)
  }
  if(change_color){
    df <- cbind(df,save_color)
  }
  if(!is.null(add)){
    if(dim(add)[2] == 1){
      n <- names(add)
      eval(parse(text=paste0("add <- data.frame(",n," = add[order(df$Lower_bound, df$Upper_bound),])")))
    }else{
      add = as.data.frame(add)[order(df$Lower_bound, df$Upper_bound),]

    }
  }

  df = df[order(df$Lower_bound, df$Upper_bound),]

  ### on ordonne les groupes en fonction de la médiane basse des groupes
  order = df  %>% group_by(Cluster) %>%
    summarise(lower_median = median(Lower_bound),upper_median = median(Upper_bound)) %>% arrange(lower_median)

  df_order = {}
  add_order = {}
  for(i in order$Cluster){
    if(!is.null(add)){
      if(dim(add)[2] == 1){
        add_order = c(add_order, as.character(add[df$Cluster == i,]))
      }else{
        add_order = rbind(add_order, add[df$Cluster == i,])
      }
    }
    df_order = rbind(df_order, df[df$Cluster == i,])
  }
  if(!is.null(add)){
    if(dim(add)[2] == 1){
      eval(parse(text=paste0("add <- data.frame(",n," = add_order)")))
    }else{
      add = as.data.frame(add_order)
    }
  }
  df = df_order
  df$Observations = factor(df$Observations, unique(df$Observations))


  ### on récupère les premières et dernières observations de chaque groupe
  intersection = df %>% group_by(Cluster) %>% summarise(First_observations = first(Observations),
                                                        Last_observations = last(Observations))

  ### On ordonne
  intersection_order = {}
  for(i in order$Cluster){intersection_order = rbind(intersection_order,
                                                     intersection[intersection$Cluster == i,])}
  intersection = intersection_order
  intersection$Cluster = factor(intersection$Cluster, levels = intersection$Cluster)

  ### On calcule un intervalle de confiance pour les dates
  intervalle_inf = intervalle_sup = {}
  for(i in unique(df$Cluster)){
    intervalle_inf[i] = arrondi(quantile(df[df$Cluster == i,]$Lower_bound,.10),0)
    intervalle_sup[i] = arrondi(quantile(df[df$Cluster == i,]$Upper_bound,.90),0)
  }

  intersection$borne_inf = intervalle_inf
  intersection$borne_sup = intervalle_sup

  ### On récupère la position du changement de groupe
  intersection_place = which(df$Observations %in% intersection$First_observations)[-1] - 1.5
  intersection_x_min = c(-.5,intersection_place)
  intersection_x_max = c(intersection_place, length(levels(df$Observation))-.5)

  nb_cluster <- length(unique(df$Cluster))
  if(nb_cluster > 1){
    intersection$x_min = intersection_x_min
    intersection$x_max = intersection_x_max
  }else{
    intersection$x_min = c(-.5)
    intersection$x_max = c(length(levels(df$Observation))-.5)
  }

  #### On définit les couleurs pour le graph
  if(change_color){
    color <- rainbow_hcl(length(unique(df$save_color)),
                         c = 80, l = 60, start = 0,
                         end = 360*(length(unique(df$save_color))-1)/(length(unique(df$save_color))))
  }else{
    color <- rainbow_hcl(length(unique(df$Cluster)),
                         c = 80, l = 60, start = 0,
                         end = 360*(length(unique(df$Cluster))-1)/(length(unique(df$Cluster))))
  }

  #graphique
  #fonction pour faire des lignes verticales
  vlinebis <- function(x = 0, color = "#666666"){
    l = {}
    for(i in 1:length(x)){
      l[[i]] = list(type = "line", y0 = 0, y1 = 1, yref = "paper",
                    x0 = x[i], x1 = x[i], line = list(color = color))
    }
    return(l)
  }

  #on définit les axes
  ay <- list(
    title = "Years",
    ticks = "outside",
    dtick = 200)

  #add_text
  if(!is.null(add)){
    col_add_text <- NULL
    for(j in 1:dim(add)[1]){
      add_text <- NULL
      for(i in 1:dim(add)[2]){
        if(dim(add)[2] == 1){
          add_text <- paste0("</br> ",names(add)[i]," : ", as.character(add[,i][j]))
        }else{
          add_text <- paste0(add_text,"</br> ",names(add)[i]," : ", as.character(add[,i][j]))
        }
      }
      col_add_text <- c(col_add_text,add_text)
    }
    df$add_text <- col_add_text
  }else{
    df$add_text <- ""
  }

  if(reorder_color){
    #re order
    order_tmp <- as.numeric(unique(df$Cluster))
    Cluster_tmp <- as.numeric(df$Cluster)
    for(i in 1:length(order_tmp)){
      if(i != order_tmp[i] && i < order_tmp[i]){
        Cluster_tmp[Cluster_tmp == i] <- 0
        Cluster_tmp[Cluster_tmp == order_tmp[i]] <- i
        Cluster_tmp[Cluster_tmp == 0] <- order_tmp[i]
        order_tmp <- unique(Cluster_tmp)
      }
    }
    df$Cluster <- as.factor(Cluster_tmp)
    order$Cluster <- as.factor(unique(Cluster_tmp))
  }


  if(change_color){
    save_color <- df$save_color
    name <- as.character(df$save_color)
  }else{
    save_color <- df$Cluster
    name <- Int2Alpha(as.numeric(as.character(df$Cluster)))
  }

  #on affiche le plot
  if(nb_cluster > 1){

    a = vlinebis(c(-.5,intersection_place,length(levels(df$Observation))-.5))
    b = line(x0 = intersection$x_min, x1 = intersection$x_max,
             y0 = intersection$borne_inf, y1 = intersection$borne_inf)
    c = line(x0 = intersection$x_min, x1 = intersection$x_max,
             y0 = intersection$borne_sup, y1 = intersection$borne_sup)

    Cl <- save_color
    p = plot_ly(df) %>%
      add_segments(x = df$Observations, xend = df$Observations,
                   y = df$Lower_bound, yend = df$Lower_bound + ((df$Upper_bound - df$Lower_bound)/2),
                   showlegend = T, color = save_color, colors = color, name=paste("Cluster", name),
                   legendgroup = paste("Cluster", name),
                   hoverinfo = 'text',
                   text = paste("</br> Observation :", df$Observations,
                                "</br> Lower bound :", df$Lower_bound,
                                "</br> Upper bound :", df$Upper_bound,
                                "</br> Cluster :", Int2Alpha(as.numeric(as.character(df$Cluster))),
                                df$add_text
                   )) %>%
      add_segments(x = df$Observations, xend = df$Observations,
                   y = df$Lower_bound + ((df$Upper_bound - df$Lower_bound)/2), yend = df$Upper_bound,
                   showlegend = F, color = save_color, colors = color,
                   legendgroup = paste("Cluster", name),
                   hoverinfo = 'text',
                   text = paste("</br> Observation :", df$Observations,
                                "</br> Lower bound :", df$Lower_bound,
                                "</br> Upper bound :", df$Upper_bound,
                                "</br> Cluster :", Int2Alpha(as.numeric(as.character(df$Cluster))),
                                df$add_text
                   )) %>%
      add_markers(x = intersection$First_observations, y = intersection$borne_inf,
                  color = I('#666666'), marker = list(size = 0.01),
                  hoverinfo = "text", showlegend = F,
                  text = paste("</br> Interval : [", intersection$borne_inf,
                               ", ", intersection$borne_sup, "]", sep = "")) %>%
      add_markers(x = intersection$First_observations, y = intersection$borne_sup,
                  color = I('#666666'), marker = list(size = 0.01),
                  hoverinfo = "text", showlegend = F,
                  text = paste("</br> Interval : [", intersection$borne_inf,
                               ", ", intersection$borne_sup, "]", sep = "")) %>%
      add_markers(x = intersection$Last_observations, y = intersection$borne_inf,
                  color = I('#666666'), marker = list(size = 0.01),
                  hoverinfo = "text", showlegend = F,
                  text = paste("</br> Interval : [", intersection$borne_inf,
                               ", ", intersection$borne_sup, "]", sep = "")) %>%
      add_markers(x = intersection$Last_observations, y = intersection$borne_sup,
                  color = I('#666666'), marker = list(size = 0.01),
                  hoverinfo = "text", showlegend = F,
                  text = paste("</br> Interval : [", intersection$borne_inf,
                               ", ", intersection$borne_sup, "]", sep = "")) %>%
      layout(shapes = c(a,b,c))
  }else{

    a = vlinebis(c(-.5,length(levels(df$Observation))-.5))

    Cl <- NULL
    p = plot_ly(df) %>%
      add_segments(x = df$Observations, xend = df$Observations,
                   y = df$Lower_bound, yend = df$Lower_bound + ((df$Upper_bound - df$Lower_bound)/2),
                   showlegend = T, color = save_color, colors = color, name=paste("Cluster", name),
                   legendgroup = paste("Cluster", name),
                   hoverinfo = 'text',
                   text = paste("</br> Observation :", df$Observations,
                                "</br> Lower bound :", df$Lower_bound,
                                "</br> Upper bound :", df$Upper_bound,
                                df$add_text
                   )) %>%
      add_segments(x = df$Observations, xend = df$Observations,
                   y = df$Lower_bound + ((df$Upper_bound - df$Lower_bound)/2), yend = df$Upper_bound,
                   showlegend = F, color = save_color, colors = color,
                   legendgroup = paste("Cluster", name),
                   hoverinfo = 'text',
                   text = paste("</br> Observation :", df$Observations,
                                "</br> Lower bound :", df$Lower_bound,
                                "</br> Upper bound :", df$Upper_bound,
                                df$add_text
                   )) %>%
      layout(shapes = c(a))
  }

  plot <- p
  if(!is.null(density)){
    save_dens <- df$density
    df$density <- df$density / sum(df$density)

    line <- function(x0, x1, y0, y1, color = "#666666",width = 2){
      l = {}
      for(i in 1:length(x0)){
        l[[i]] = list(type = "line", y0 = y0[i], y1 = y1[i], yref = "y2",
                      x0 = x0[i], x1 = x1[i], line = list(color = color, width = width))
      }
      return(l)
    }

    l <- line(x0 = 0:(length(df$density)-2),
              y0 = df$density[1:(length(df$density)-1)],
              x1 = 1:(length(df$density)-1),
              y1 = df$density[2:length(df$density)],
              color = "#E16A86",width = 2.5)

    q <- plot_ly(df) %>%
      add_markers(x = df$Observation, y = df$density, showlegend = F,
                  hoverinfo = 'text', color = I("#E16A86"), marker = list(size = 2),
                  hovertext = paste("</br> Observation :", df$Observations,
                                    "</br> Proportion :", arrondi(df$density,3),
                                    "</br> Variable :", arrondi(save_dens,3),
                                    df$add_text)
                  ) %>%
      layout(shapes = na.omit(c(a,l)))

    plot <- subplot(
      p,q, nrows = 2, heights = c(0.7,0.3), shareX = TRUE, shareY = FALSE
      ) %>%
      layout(title = "Timerange",
             yaxis = list(title = "Years",zeroline = FALSE),
             yaxis2 = list(title="Density"),
             xaxis = list(title = "Observations")
      )
    density <- df$density
    names(density) <- df$Observation
  }

  return(structure(list(order = order, density = density, plot = plot, reorder_cluster = Cl),
                   class = c("temp_obj", "list")))
}

arrondi<-function (x, acc = 0)
{
  .local <- function(x, acc) {
    x <- x * (10^acc)
    ifelse(abs(x%%1 - 0.5) < .Machine$double.eps^0.5, ceiling(x)/(10^acc),
           round(x)/(10^acc))
  }
  if (is.data.frame(x))
    return(data.frame(lapply(x, .local, acc)))
  .local(x, acc)
}


Int2Alpha <- function(int){
  #====================================#
  # convert integer into letter
  let <- function(alphabet) function(i) {
    base10toA <- function(n, A) {
      stopifnot(n >= 0L)
      N <- length(A)
      j <- n %/% N
      if (j == 0L) A[n + 1L] else paste0(Recall(j - 1L, A), A[n %% N + 1L])
    }
    vapply(i-1L, base10toA, character(1L), alphabet)
  }

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

Alpha2Int <- function(alpha){
  #====================================#
  # convert letter into integer
  get.letters <- function(length.out, case = 'upper'){
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
    max.char <- ceiling(log(length.out, 26))
    grp <- rep(1:max.char, 26^(1:max.char))[1:length.out]
    unlist(lapply(unique(grp), function(n) make.chars(length(grp[grp == n]), case = case, n.char = n)))
  }

  alpha <- toupper(alpha)
  letter <- get.letters(800)
  first <- c()
  second <- c()
  for(i in 1:length(alpha)){
    first <- c(first,str_extract(alpha[i], "[A-Z]+"))
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

elbow_finder <- function(x_values, y_values) {
  # Max values to create line
  max_x_x <- max(x_values)
  max_x_y <- y_values[which.max(x_values)]
  max_y_y <- max(y_values)
  max_y_x <- x_values[which.max(y_values)]
  max_df <- data.frame(x = c(max_y_x, max_x_x), y = c(max_y_y, max_x_y))

  # Creating straight line (between max and min)
  fit <- lm(max_df$y ~ max_df$x)

  # Distance from point to line
  distances <- c()
  for(i in 1:length(x_values)) {
    distances <- c(
      distances,
      abs( coef(fit)[2] * x_values[i] - y_values[i] + coef(fit)[1] ) / sqrt( coef(fit)[2]^2 + 1^2 )
    )
  }

  return(distances)
}

convert <- function(x){
  res <- c()
  for(i in 1:length(x)){
    split <- strsplit(as.character(x[i]), split = "")
    alpha <- split[[1]][length(split[[1]])]
    if(split[[1]][1] == "-"){stop("We can't convert ",as.character(x[i]),
                                  " to date. Enter the date as a number.")}
    if(alpha == "a"){unit <- "00"}
    if(alpha == "b"){unit <- "25"}
    if(alpha == "c"){unit <- "50"}
    if(alpha == "d"){unit <- "75"}
    if(alpha != "a" && alpha != "b" && alpha != "c" && alpha != "d"){unit <- "00"}
    if(length(split[[1]]) == 3){
      res <- c(res,as.numeric(paste0(split[[1]][1],split[[1]][2],unit)))
    }else{
      if(alpha != "a" && alpha != "b" && alpha != "c" && alpha != "d"){
        if(length(split[[1]]) == 2){
          res <- c(res,as.numeric(paste0(split[[1]][1],split[[1]][2],unit)))
        }
        if(length(split[[1]]) == 1){
          res <- c(res,as.numeric(paste0(split[[1]][1],unit)))
        }
      }else{
        res <- c(res,as.numeric(paste0(split[[1]][1],unit)))
      }
    }
  }
  return(res-100)
}
