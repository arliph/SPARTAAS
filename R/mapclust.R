# ClustCutMap!
#
#

estimateDensityMulti <- function(data){
  cat("Estimation of multidimensionnal density ..")
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
    xg <- tapply(ww*zz*xx,g,sum,na.rm=TRUE)/tapply(ww*zz,g,sum,na.rm=TRUE)
    yg <- tapply(ww*zz*yy,g,sum,na.rm=TRUE)/tapply(ww*zz,g,sum,na.rm=TRUE)
    for(i in 1:length(xx)) {
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
    return(list(n=g,mat=cbind(n=n,xg=arrondi(xg,2),yg=arrondi(yg,2),zg=arrondi(pb*100,2)),
                InertiaWithin=InertiaWithin,
                InertiaWithintot=sum(InertiaWithin),
                InertiaBetween=inertiaTot - sum(InertiaWithin),
                inertiaTot=inertiaTot,
                nsp=res))
  }

  # prepare data
  Xsta <- x
  Ysta <- y
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
  cat("Estimation of univariate density ..")
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

mapclust <-
  function (coord, var, label = "Nolabel", iter = 20, Plabel = TRUE, lonlat = TRUE, positive_var = FALSE, n = NULL)
  {
    #===================================================================================
    #add two param: transfo_Z and k
    # DENDROGRAMME CUT AND MAP
    #
    # Author : A. Coulon (Universite de Nantes)
    # Last update : 17 April 2018
    #
    # Arguments:
    # x	      	      The vector with x-coordinate (latitude)
    # y	              The vector with y-coordinates (longitude)
    # var	      	    The regionalised variable of interest (vector if univariate analysis, data.frame or matrix max: 6dimension)
    # label(optional)	The names of the samples or an id (factor)
    # Plabel          logical parameter for activate or not the Print of labels on dendrogram
    #
    # Information:
    # dlim is the select minimum distance from sample to patch centre: to
    # identify patches (units are those of coordinates)
    #
    #===================================================================================

    cuttree2 <- function(classification, index, coul)
    {
      #===================================================================================
      # COUPE DENDROGRAMME
      #
      # Author : A. Coulon (Universite de Nantes)
      # Last update : 17 April 2018
      #
      # Arguments:
      # classification	  The object return by mapclust
      # index             The index of the cut dlim in the classification$dlim vector
      #
      #===================================================================================

      cuttree <- classification$dendrogram_ggplot
      cuttree <- cuttree + geom_hline(yintercept = classification$dlim[index] + 0.03 , col = "black" , lwd = 1)
      if(classification$Plabel){
        cuttree <- cuttree + suppressWarnings(theme(axis.text.x = element_text(color =  coul,
                                                          hjust = 1, vjust = 0.25, size = 9, angle = 90)))
      }else{
        cuttree <- cuttree + theme(axis.text.x = element_blank())
      }

      return(cuttree = cuttree)
    }

    index2 <- function(classification, historique, nb_grp){

        V_NB_GRP <- c()
        index <- NULL
        for(i in 1:length(classification$dlim)){
          V_NB_GRP <- c(V_NB_GRP,length(unique(historique[,i])))
          if(nb_grp == length(unique(historique[,i]))){
            cutdlim <- classification$dlim[i]
            index <- i
          }
        }

        if(is.null(index)){
          message("There is not this number of group in the partition. We cut to the closest. possible number:")
          print(unique(V_NB_GRP))
          diff <- V_NB_GRP - nb_grp
          index <- which.min(abs(diff))
          cutdlim <- classification$dlim[index]
        }
      return(index)
    }

    #mapclust-----------------------------------------------------
    if(iter <= 9){stop("You must run at least 10 iterations.")}
    if(iter > 101){message("You try to run a high number of iterations. This will take a long time.")}
    if(iter > 201){message("Only run if you are sure you need so many iterations.")}
    if(iter >= 1000){stop("The number of iterations is too high. (<1000)")}#rm?
    if(dim(coord)[2] != 2){stop("coord: must be a data.frame or a matrix with two columns!")}
    if(!is.numeric(coord[,1])){stop("coord: must be numeric!")}
    if(!is.numeric(coord[,2])){stop("coord: must be numeric!")}
    if(!is.null(n) && !is.numeric(n)){stop("n: must be numeric!")}
    if(!is.null(n) && n <= 1){stop("n: must be higher than 1!")}
    if(!is.logical(positive_var)){stop("positive_var: must be logical!")}
    if(dim(as.data.frame(var))[2] > 1){
      for(i in 1:dim(as.data.frame(var))[2]){
        if(!is.numeric(var[,i])){stop("var: must be numeric!")}
      }
    }else{
      if(!is.numeric(var)){stop("var: must be numeric!")}
    }
    if(!is.character(label) && !is.factor(label)){stop("label: must be character (factor)!")}
    X <- coord[,1]
    Y <- coord[,2]
    binf = 0
    bsup = max(dist(cbind(X,Y)))
    pas = (max(dist(cbind(X,Y)))/iter)
    #test param

    ref <- length(X)
    if(length(Y) != ref){stop("'x' and 'y' have not the same length")}

    if(dim(as.data.frame(var))[2] > 1){
      if(length(var[,1]) != ref){stop("'coord' and 'var' have not the same length")}
    }else{
      if(length(var) != ref){stop("'coord' and 'var' have not the same length")}
    }
    if(length(label) != ref && label[1] != "Nolabel"){stop("'x' and 'label' have not the same length")}

    save_label <- label

    if(dim(as.data.frame(var))[2] > 6){
      if (rstudioapi::hasFun("showDialog")) {
        showDialog("warning","var can't have more than 6 dimension",url="")
      }else{
        message("var can't have more than 6 dimension\n")
      }
      stop("var can't have more than 6 dimension")
    }

    res <- NULL
    if(dim(as.data.frame(var))[2] > 1){
      #size uniform
      Size <- rep(1,length(var[,1]))
      multi = TRUE
    }else{
      Size <- var
      multi = FALSE
    }
    xg <- c()

    ctrl = positive_var

    if(!ctrl){
      if(dim(as.data.frame(var))[2] > 1){
        Z <- estimateDensityMulti(var)
      }else{
        Z <- transfo_z(var)
      }

    }else{

      if(dim(as.data.frame(var))[2] > 1){

          message("var is multidimensionnal. You must transform it!
                  var variables change forced!\n")

        Z <- estimateDensityMulti(var)
      }else{
        if(min(var)<0){

            message("There are negatives values in var. The data are not like a count, probability or density.
                    var change forced.\n")

          Z <- transfo_z(var)
        }else{
          Z <- var
          cat("Use var as density ... Done")
        }
      }
    }

    data <- data.frame(X,Y,Z)

    if (bsup <= binf) {stop("The upper bound (bsup) value must be higher than the lower bound (binf).")}

    slabel <- as.character(label)

    if (Plabel) {
      if ( label[1] == "Nolabel") {
        #Plabel <- FALSE
        label <- c(1:length(X))
      }else{
        #si l'utilisateur renseigne un nom
        label <- as.character(label)
      }
    }

    id <- c(1:length(X))
    data=cbind(data,id)
    data=as.data.frame(data)

    if( bsup > max(dist(data[,1:2])) ) {
      cat("The maximal value of dlim is: ",arrondi(max(dist(data[,1:2])),2))
      cat("the upper bound (bsup) have been changed in order to not run excessively high dlim.")
      bsup <- max(dist(data[,1:2]))
    }

    #remise a zero des objects
    names <- c()
    #initialisation matrice historique des groupes
    historique <- matrix(ncol=0, nrow=length(X))

    cat("\nInitialisation and calcul of the distance matrix ..")
    DiMatrix <- dist_m(data)
    cat(". Done\n")
    cat("Running Progress:\n")
    #creation barre de progression(pb)
    pb <- txtProgressBar( min = 0, max = arrondi((bsup - binf) / pas,0) , style = 3 )

    pre_k <- 0
    k <- 1
    pass <- 0
    grp <- c(1)
    data <- cbind(data,grp)
    dlim <- bsup
    #running spatialpatches (WOILLEZ et al., 2009)
    Pch <- spatialpatches2(X , Y , Z , w = rep(1 , length(Y)) , Lim.D = dlim)
    grp1 <- Pch$n

    names <- c(arrondi(dlim,2))
    historique <- cbind(historique,grp1)

    data$grp=grp

    if (!(length(unique(data$grp)) <= 1)) {
      cat("\nThe maximal value of dlim is not enough to create one patch at the top of the tree.
          The groups may be merged too early. Try higher bsup.\n")
    }


    while ( dlim >= binf) {
      # update progress bar
      setTxtProgressBar(pb, k)

      #if k have not been increment
      if(k == pre_k){
        dlim <- nextdlim
      }else{
        dlim = bsup - k * pas
        #save value of k for this loop. need for check next loop time
        pre_k <- k
      }

      if (arrondi(dlim,2) > 0.01) {

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
            if(eval(parse(text=paste0("length(data",a,"[,1]) > 1")))){
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
              if(diff2_3 > diff1_2 && diff1_2 > diff1_3){
                eval(parse(text=paste0("data",a,"$grp <- ifelse(data",a,"$grp==",a,"2,",a,"2,",a,"1)")))
              }else{
                if(diff2_3 > diff1_2 && diff1_2 < diff1_3){
                  eval(parse(text=paste0("data",a,"$grp <- ifelse(data",a,"$grp==",a,"3,",a,"2,",a,"1)")))
                }else{
                  eval(parse(text=paste0("data",a,"$grp <- ifelse(data",a,"$grp==",a,"1,",a,"1,",a,"2)")))
                }
              }
              nextdlim <- dlim - (pas/3)
              #pass take value 1 in order to re run one time the loop witjhout change k and replace dlim by nextdlim.
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
          namegrp <- c()
          namegrp <- as.factor(grp_tmp$grp)
          levels(namegrp) <- c(1:length(levels(namegrp)))
          data$grp <- as.numeric(namegrp)
        }else{
          # only one patch in all so we keep going with all the data
          Pch <- spatialpatches2(X , Y , Z , w = rep(1 , length(Y)) , Lim.D = dlim)
          grp <- Pch$n
          data$grp <- grp

          if( length(unique(data$grp)) >= 3 ){
            # if 3 or more groups
            XX <- data[,1]
            YY <- data[,2]
            ZZ <- data[,3]

            WW <- rep(1 , length(data[,2]))
            xg <- tapply(WW*ZZ*XX,data$grp,sum,na.rm=T)/tapply(WW*ZZ,data$grp,sum,na.rm=T)
            yg <- tapply(WW*ZZ*YY,data$grp,sum,na.rm=T)/tapply(WW*ZZ,data$grp,sum,na.rm=T)
            diff1_2 <- abs(xg[1]-xg[2])
            diff1_3 <- abs(xg[1]-xg[3])
            diff2_3 <- abs(xg[2]-xg[3])
            if(diff2_3 > diff1_2 && diff2_3 > diff1_3){
              data$grp <- ifelse(data$grp==2,2,1)
            }else{
              data$grp <- ifelse(data$grp==1,1,2)
            }
          }

        }

        if(arrondi(dlim,2) > 0){
          j <- length(historique[1,])
          #compare last and new compo
          v <- data$grp==historique[,j]
          somme <- sum(v %in% FALSE)
          if(somme != 0){
            #if data$grp is different than the last compo hist
            #add the dlim to name and the hist
            names <- c(names , arrondi(dlim,2))
            historique <- cbind(historique , data$grp)
          }else{
            if(arrondi(dlim,2) < names[length(names)]){
              #IF no difference : replace the last dlim name by the new
              names[length(names)] <- arrondi(dlim,2)
            }
          }
        }
      }
      if(pass == 1){
        #if pass egal to 1 don't increment K and re-run one time the while loop. pass egal 0 for next time of loop
        pass <- 0
      }else{
        k = k + 1
      }
    }#finwhile

    close(pb)
    dlim <- as.numeric(colnames(historique)[length(colnames(historique))])

    historique <- cbind(grp1,historique)
    names <- c(arrondi(bsup,2),names)
    historique <- as.matrix(historique)
    colnames(historique) <- names
    res$dlim <- names


    #silhouette
    Average_Sil_Width <- c()
    for(i in 1:length(historique[1,])){
      cluster <- historique[,i]
      if(length(unique(cluster))>1){
        sil <- silhouette(cluster,DiMatrix)
        avg <- mean(sil[,3])
        Average_Sil_Width <- c(Average_Sil_Width,avg)
      }
    }
    number <- apply(historique,2,unique)
    number <- lapply(number,length)

    Number_of_groups <- as.vector(unlist(number))[number>1]
    Average_Sil_Width <- rbind(Number_of_groups,Average_Sil_Width)

    historique <- cbind(id,historique)
    hist<-c()

    cpy_historique <- historique
    cpy_historique<-cbind(label,cpy_historique)

    #On trie l'historique par la colonne du dlim minimal
    historique <- historique[order(historique[,length(historique[1,])],decreasing=FALSE), ]
    cpy_historique <- cpy_historique[order(as.numeric(cpy_historique[,length(cpy_historique[1,])]),decreasing=F), ]


    label_ini_grp <- c()
    label_tmp <- c()
    label_avant <- c()
    label_grp <- c()

    #for(b in 2:(length(cpy_historique[,1]))){
    #  c=b-1
    #  if(!(FALSE %in% (historique[b,2:length(historique[1,])]==historique[c,2:length(historique[1,])]))){
    #    #Doublon
    #    label_tmp <- as.character(cpy_historique[b,1])
    #    label_avant <- as.character(cpy_historique[c,1])
    #    if(is.null(label_grp)){
    #      label_grp <- c(label_avant,label_tmp)
    #    }else{
    #      label_grp <- c(label_grp,label_avant,label_tmp)
    #    }
    #  }else{
    #    #Pas doublon
    #    hist<-rbind(hist,historique[c,])
    #    label_grp <- c(label_grp,as.character(cpy_historique[c,1]))
    #    label_grp <- paste(unique(label_grp), collapse = ".")
    #    label_ini_grp<-c(label_ini_grp,label_grp)
    #    label_grp <- c()
    #  }
    #}

    #add the last elt
    #hist <- rbind(hist,historique[length(historique[,1]),])
    #label_grp <- paste(unique(as.character(cpy_historique[length(historique[,1]),1])), collapse = ".")
    #label_ini_grp <- c(label_ini_grp,label_grp)
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

    #We cut the label colunms (Not use here because we work on the entire cluster the label and id of
    #each members of one cluster are still save in data)
    hist <- as.matrix(hist)
    hist <- hist[,2:length(historique[1,])]


    ecart <- ecartmatrice(hist)
    #if un seul groupe pas de div gerer ce cas

    ##Les segments
    #creation colonne x
    xcoord <- c()
    for (b in 1:length(hist[,1])) {
      xcoord <- c(xcoord , b)
    }
    hist <- as.matrix(hist)
    #creation colonne y
    ycoord <- rep(binf, length(hist[,1]))

    #creation d'une matrice ou on stockera les coordonnees pour les segments du drendrogramme
    coord <- cbind(x=xcoord, y=ycoord)

    #initialisation des segments
    segment <- c()

    #------------------------
    tmp_id <- -seq(1,length(hist[,1]))
    dend <- list(order=rep(NA,length(hist[,1]))
                 ,order.lab=rep(NA,length(hist[,1]))
                 ,height=rep(NA,length(hist[,1])-1)
                 ,dc=NA
                 ,merge=matrix(NA,nrow=length(hist[,1])-1,ncol=2)
    )
    #------------------------

    #boucle se repetant (nombre de grp initiaux -1)
    #Pour 6 grps il faut (6-1)*3=15 segments. n grp il faut (n-1)*3 segments.
    for (a in 1:(length(hist[,1])-1)) {

      #On recherche min (premiere fusion)
      mini <- min(ecart)
      #On recherche les indices
      index_vect <- which.min(ecart)
      index <- arrayInd(index_vect,c(length(ecart[,1]),length(ecart[,1])))
      grp1 <- index[1]
      grp2 <- index[2]

      #------------------------
      dend$height[a] <- res$dlim[length(res$dlim)-mini]
      dend$merge[a,1] <- tmp_id[grp1]
      dend$merge[a,2] <- tmp_id[grp2]
      tmp_id[grp1] <- tmp_id[grp2] <- a

      #ecart[,grp2] <- Inf
      #ecart[grp2,] <- Inf
      #------------------------

      #On recupere dlim dans l'objet res$dlim qui contient tous les dlim utilises. length(...)-mini pour
      #partir de la fin et remonter. mini est le nombre de dlim passer vers la gauche.

      if( length(res$dlim)-mini != 0){
        dlim_tmp <- res$dlim[length(res$dlim)-mini]

        #On construit 3 segments (deux qui monte et un qui fusionne branche)

        # segments verticaux
        # x = l'abscisse du tab coord pour grp correspondant  y = coord[grp_,2]    xend = x    yend = dlim_tmp
        segment <- rbind(segment , c(coord[grp1,1] ,  coord[grp1,2] ,  coord[grp1,1] ,    dlim_tmp))
        segment <- rbind(segment , c(coord[grp2,1] ,  coord[grp2,2] ,  coord[grp2,1] ,    dlim_tmp))

        # segment horizontal
        #                        x = coord[grp1,1]  y = dlim_tmp   xend = coord[grp2,1]  yend = dlim_tmp
        segment <- rbind(segment,  c(coord[grp1,1] ,    dlim_tmp ,        coord[grp2,1] ,    dlim_tmp))

        ## On sauvegarde la hauteur (y / dlim) a laquelle les grps sont rendus
        #MAJ valeur de y (dlim) grp1
        coord[grp1,2] <- dlim_tmp
        #MAJ ligne grp1, valeur en X du nouveau grp (fusion) moyenne des deux
        coord[grp1,1] <- (coord[grp1,1]+coord[grp2,1])/2
        #On marque la ligne du grp2 en mettant x=0, pour ne plus la prendre en compte dans prochaine
        # fusion puisque ce grp a fusionne avec grp1
        coord[grp2,1] <- 0
        #On elimine la ligne/colonne du grp2 dans les valeurs (historique grp)
        ecart[,grp2] <- Inf
        ecart[grp2,] <- Inf
      }
    }

    #------------------------
    dend$order <- seq(1,length(hist[,1]))
    dend$order.lab <- label_ini_grp
    dend$dc <- 0

    class(dend) <- c("diana","twins")
    #------------------------

    groupes=1
    dlim=1
    xend=1
    yend=1
    #On rajoute un dernier segment pour la racine de l'arbre
    segment <- rbind(segment , c(max(coord[,1]) , max(coord[,2]) , max(coord[,1]) , (max(res$dlim)+pas)))
    segment <- data.frame(segment)
    #On prepare tableau segment
    names(segment) <- c("groupes" , "dlim" , "xend" , "yend")
    row.names(segment) <- NULL


    #Plabel=F
    #On creer un tableau pour les noms (implant au bout des branches y = 0)
    label <- matrix(ncol = 3, nrow = length(hist[,1]))
    for (a in 1:length(hist[,1])) {
      label[a,1] <- a
      label[a,2] <- 0
      label[a,3] <- label_ini_grp[a]
    }
    label <- data.frame(label)

    plot <- NULL
    plot1 <- NULL

    #On dessine les segments et on ajout les noms
    auto.haut <- max(res$dlim) + pas / 10
    tmp <- 1
    for (b in 2:length(hist[,1])) {
      tmp <- c(tmp , b)
    }

    #Dendrogram
    plot1 <- ggplot(segment , col='white') +
      geom_segment(data = segment , aes(x = groupes, y = dlim, xend = xend, yend = yend)) +
      ylim(0 , max(res$dlim) + pas)
    if (Plabel) {
      plot1 <- plot1 + scale_x_continuous(breaks = tmp ,
                                          labels = as.character(label[,3]) )
      plot1 <- plot1 + theme(axis.text.x = element_text(face = "bold" , hjust = 1 , vjust = 0.25 , size = 9 ,
                                                        angle = 90))
    }
    plot1 <- plot1 + labs(title = "Hierarchical Clustering" , subtitle = "Mapclust")

    plot1 <- plot1 + theme(
      plot.margin = margin(0.5, 1, 0.5, 1, "cm")
    )
    res$dendrogram_ggplot <- plot1

    groupes=1
    dlim=1
    xend=1
    yend=1
    #Dendrogramme affiche avec indication pour le clic de cut.
    plot <- ggplot(segment , col = 'white') +
      geom_segment(data = segment,aes(x = groupes, y = dlim, xend = xend, yend = yend)) +
      ylim(0 , max(res$dlim) + pas)
    #Title subtitle caption
    plot <- plot + labs(title = "Hierarchical Clustering", subtitle = "Mapclust")
    #Axes x
    plot <- plot + scale_x_continuous(breaks = c())

    if (Plabel) {
      plot <- plot + geom_text(data = label, aes(x = as.numeric(paste(label[,1])),
                                                 y = as.numeric(paste(label[,2])),
                                                 label = as.character(label[,3])),
                               fontface = "bold", nudge_x = -0.12, colour = "gray32",
                               hjust = 0, vjust = 0.25, size = 3, angle = 90)
    }
    plot <- plot + theme(
      panel.grid = element_line(colour = "grey"),
      plot.margin = margin(0.5, 1, 0.5, 1, "cm")
    )
    #print(plot)
    minisil <- min(Average_Sil_Width)
    if(minisil > 0 ){
      minisil <- 0
    }else{
      minisil <- -1
    }
    #plot eval sil
    dataplot <- t(Average_Sil_Width)
    dataplot <- as.data.frame(dataplot)
    dataplot <- na.omit(dataplot)
    colnames(dataplot) <- c("Number_of_groups","Average_Sil_Width")
    p <- ggplot(dataplot,aes(x=Number_of_groups, y=Average_Sil_Width)) +
      geom_point(aes(color=Average_Sil_Width)) + scale_colour_gradientn(colours=c("brown1","chartreuse3")) +
      geom_line(linetype = 3) +
      scale_x_continuous(breaks = dataplot$Number_of_groups,
                         labels = dataplot$Number_of_groups) +
      ylim(minisil,1) +
      labs(title = "Average Silhouette Width", subtitle = "Partition evaluation") +
      theme(legend.position="none")

    historique<-historique[order(historique[,1],decreasing=F), ]
    Sum_of_Square <- WSS(historique[,-1],data)
    Itot <- calctot(data)
    if(dim(as.data.frame(var))[2] > 1){
      Moran <- "Impossible to calculate Moran index for multivariate data"
    }else{
      Moran <- IndiceMoran(historique[,-1],data)
    }

    wss <- c()

    for (i in Sum_of_Square) {
      if(length(i)>1){
        wss <- c(wss,sum(i))
      }
    }

    names(Sum_of_Square) <- res$dlim

    Within_Sum_of_Square <- rbind(Number_of_groups,wss)
    color_wss <- elbow_finder(Number_of_groups,wss)
    #plot eval WSS
    dataplot <- t(Within_Sum_of_Square)
    dataplot <- as.data.frame(dataplot)
    colnames(dataplot) <- c("Number_of_groups","Within_Sum_of_Square")
    q <- ggplot(dataplot,aes(x=Number_of_groups, y=Within_Sum_of_Square)) +
      geom_point(aes(color=color_wss)) + scale_colour_gradientn(colours=c("brown1","chartreuse3")) +
      geom_line(linetype = 3) +
      scale_x_continuous(breaks = dataplot$Number_of_groups ,
                         labels = dataplot$Number_of_groups ) +
      ylim(0,Itot) +
      labs(title = "Within Sum of Square" , subtitle = "Partition evaluation") +
      theme(legend.position="none")

    # Nouvelle page
    grid.newpage()
    # Creer la mise en page : nrow = 2, ncol = 2
    pushViewport(viewport(layout = grid.layout(2, 2)))
    # Une fonction pour definir une region dans la mise en page
    define_region <- function(row, col){
      viewport(layout.pos.row = row, layout.pos.col = col)
    }
    # Arranger les graphiques dans la fenetre
    print(plot, vp = define_region(1:2, 1))
    print(p, vp = define_region(1, 2))
    print(q, vp = define_region(2, 2))

    if(is.null(n)){
      if (rstudioapi::hasFun("showQuestion")) {
        showDialog("Cut Tree", "You must click on the dendrogram (plot window) to select where the tree should be cut.
               You will not be able to execute anything else until you have selected a dlim.",url="")
      }


      message("
            |=======================================================|
            | You must click on the dendrogram (plot window)        |
            | to select where the tree should be cut.               |
            | You will not be able to execute anything else until   |
            | you have selected a dlim.                             |
            |=======================================================|
            ")
      OK = FALSE
      while(OK == FALSE){
        #coupe arbre
        point <- data.frame(0,-10)
        loc <- as.numeric(grid.locator("npc"))
        # On adapte les coord (de la fenetre) a l axes y.
        point <- data.frame( loc[1], ( (loc[2]-0.1055195)/0.7516234  ) * ( max(res$dlim) + pas) )
        names(point) <- c("x","dlim")

        if (is.na.data.frame(point[1,2])){point <- data.frame(0,-10)}
        if (point[1,2] == -10){point[1,2] <- median(res$dlim)+0.1}

        #dlim max parmi ceux inferieur au dlim pointe par user
        diff <- res$dlim-point[1,2]
        if((TRUE %in% (diff <= 0))){
          #return index of the first value <= 0
          index <- which.max(diff <= 0)
        }else{
          index <- which.min(diff)
        }

        cutdlim <- res$dlim[index]
        nb <- length(unique(historique[,index+1]))

        if (nb <= 1){
          message("The number of patchs must be greater than 1. Try again")
        }else{
          OK = TRUE
        }

      }

      if (cutdlim <= min(res$dlim)) {
        cutdlim <- min(res$dlim)
      }
      cat("Select dlim ...",cutdlim)
    }else{
      index <- index2(res,historique,n-1)

      cutdlim <- res$dlim[index]
      nb <- length(unique(historique[,index+1]))

      if (cutdlim <= min(res$dlim)) {
        cutdlim <- min(res$dlim)
      }
      cat("Select dlim ...",cutdlim)
    }


    cat("\nCut the dendrogram ..")
    res$Plabel <- Plabel

    tmp<-data.frame(lon=X, lat=Y)
    if(!lonlat){
      tmp$lon <- (tmp$lon - mean(tmp$lon)) / max(tmp$lon, tmp$lat)
      tmp$lat <- (tmp$lat - mean(tmp$lat)) / max(tmp$lon, tmp$lat)
    }

    lineType <- ifelse(Size<0,"2, 5","1")

    #leaflet
    df <- data.frame(
      lat <- tmp$lat,
      lng <- tmp$lon,
      size <- 10*(0.1+(abs(Size)/max(abs(Size))))+5,
      color <- as.factor(historique[,index+1])
    )
    if(slabel[1] != "Nolabel"){
      df <- cbind(df,name = as.character(slabel))
    }
    if(multi){df$size <- rep(7,length(df$size))}
    df$category <- as.factor(historique[,index+1])

    # Make up some random levels. (TODO: Better example)
    couleur <- c("slateblue2", "blue", "chartreuse1","darkolivegreen2",
                 "#222222", "darkmagenta", "mediumpurple4",
                 "magenta","cyan3", "plum4","darkorange1", "pink1", "red",
                 "lightgoldenrod", "violetred3", "brown4", "tan1","yellow", "gray50",
                 "#967272","deepskyblue4","salmon", "cyan","palegreen", "green3")
    #plot(1:25,1:25,pch =15,cex=3,col=couleur)
    factpal <- colorFactor(couleur, df$category)
    coul <- factpal(as.factor(hist[,index]))

    cut <- cuttree2(res,index,coul)


    if(lonlat){
      esri <- structure(c("Esri.WorldStreetMap", "Esri.DeLorme", "Esri.WorldTopoMap",
                          "Esri.WorldImagery", "Esri.WorldTerrain", "Esri.WorldShadedRelief",
                          "Esri.WorldPhysical", "Esri.OceanBasemap", "Esri.NatGeoWorldMap",
                          "Esri.WorldGrayCanvas"),
                        .Names = c("Esri.WorldStreetMap", "Esri.DeLorme", "Esri.WorldTopoMap",
                                   "Esri.WorldImagery", "Esri.WorldTerrain", "Esri.WorldShadedRelief",
                                   "Esri.WorldPhysical", "Esri.OceanBasemap", "Esri.NatGeoWorldMap",
                                   "Esri.WorldGrayCanvas"))
      m = leaflet(df)
      for (provider in esri) {
        m <- m %>% addProviderTiles(provider, group = provider)
      }
      if(slabel[1] != "Nolabel"){
        m <- m %>% addCircleMarkers(lng = ~lng, lat = ~lat, popup = ~htmlEscape(name), radius = ~size, fillColor =  ~factpal(category),
                                    fill = TRUE, fillOpacity = 0.8, opacity=1, stroke = TRUE,weight = 1.5, color="black", dashArray=lineType
        )
      }else{
        m <- m %>% addCircleMarkers(lng = ~lng, lat = ~lat, radius = ~size, fillColor =  ~factpal(category),
                                    fill = TRUE, fillOpacity = 0.8, opacity=1, stroke = TRUE,weight = 1.5, color="black", dashArray=lineType
        )
      }

      m <- m %>%
        addLegend("bottomleft", pal = factpal, values = as.factor(historique[,index+1]),
                  title = "Cluster",
                  opacity = 1
        ) %>%
        addMiniMap(toggleDisplay = TRUE,
                   position = "bottomright"
        ) %>%
        addLayersControl(baseGroups = names(esri),
                         options = layersControlOptions(collapsed = TRUE),
                         position = "topright"
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
      if(!multi){
        m <- m %>%
          addLegendCustom(colors = c("black", "black", "black"),
                          labels = c(arrondi(min(abs(Size)),2),
                                     arrondi(min(abs(Size)) + (max(abs(Size))-min(abs(Size)))/2  ,2),
                                     arrondi(max(abs(Size)),2)),
                          sizes = c(10*(0.1+( min(abs(Size)) / max(abs(Size))))+5,
                                    10*(0.1+((min(abs(Size)) + (max(abs(Size))-min(abs(Size)))/2) /max(abs(Size))))+5,
                                    10*(0.1+( max(abs(Size)) / max(abs(Size))))+5),
                          title = "Size"
          )
      }
      if(min(Size)<0){
        m <- m %>%
          addLegendCustom(colors = c("white", "white"), labels = c("line: positive","dotted-line: negative"),
                          sizes = 1, title="Sign"
          )
      }
    }else{
      m <- leaflet(df)
      m <- m %>% addCircleMarkers(lng = ~lng, lat = ~lat, radius = ~size, fillColor =  ~factpal(category),
                                  fill = TRUE, fillOpacity = 0.8, opacity=1, stroke = TRUE,weight = 1.5, color="black", dashArray=lineType
      ) %>%
        addLegend("bottomleft", pal = factpal, values = as.factor(historique[,index+1]),
                  title = "Cluster",
                  opacity = 1
        )
    }


    cluster <- as.numeric(historique[,index+1])

    sil <- silhouette(cluster,DiMatrix)
    coull <- factpal(sort(cluster))
    if(slabel[1] != "Nolabel"){
      row.names(sil) <- slabel
    }
    plot(sil,col = coull, nmax = 500, cex.names=.8, max.strlen = 80,
         main="Silhouette plot (x = cluster, dist = Weighted Eucidiean Distance)")
    silPlot <- recordPlot()
    plot.new()
    print(cut)
    cat(". Done")
    print(m)

    out <- structure(list(X=X,
                          Y=Y,
                          var=var,
                          label=save_label,
                          density=Z,
                          cluster=cluster,
                          fullhist=historique,
                          dlim=res$dlim,
                          historique=hist,
                          Plabel=Plabel,
                          InertiaWithin = Sum_of_Square,
                          Moran=Moran,
                          dendrogram=as.hclust(dend),
                          dendrogram_ggplot=res$dendrogram_ggplot,
                          DiMatrix=DiMatrix,
                          cuttree=cut,
                          map=m,
                          silhouette=silPlot,
                          silhouetteData=sil,
                          AveSilData=Average_Sil_Width[2,],
                          AveSilPlot=p,
                          WSSPlot=q,
                          cutdlim=cutdlim,
                          lonlat=lonlat,
                          call=match.call()),class=c("mapclust_cl","list"))
    cat("\n\n\n
    Divisive hierarchical clustering

call: ")
    print(out$call)
    cat(sep = "","\npartition:
      number of cluster: ",length(unique(out$cluster)),"\n      Average_Sil = ",arrondi(mean(out$silhouetteData[,3]),2),", WSS = ",arrondi(sum(out$InertiaWithin[[which(out$dlim == out$cutdlim)]]),2),", cutdlim = ",out$cutdlim,
        "\nplot:
      .. $dendrogram      The global dendrogramm
      .. $cuttree         The cut dendrogramm
      .. $map             The map of the select partition
      .. $AveSilPlot      The average silhouette width plot (for each partition)
      .. $WSSPlot         The Within Sum of Square plot (for each partition)
      .. $silhouette      The silhouette of the select partition
value:
      .. $X               The x-coordinates data you used
      .. $Y               The y-coordinates data you used
      .. $var             The regionalized variable.s data you used
      .. $label           The label vector you used
      .. $density         The estimate density based on var. Equal to var if you used a unidimensionnal density variable
      .. $cluster         The vector of cluster of the select partition
      .. $fullhist        The composition cluster for each observation
      .. $hist            The composition cluster without duplicates (matches to leaf of the dendrogramm)
      .. $dlim            The vector of the different limit distances
      .. $cutdlim         The select dlim for the cut
      .. $InertiaWithin   The Sum of Square for each patch of each partition
      .. $DiMatrix        The matrix of Weighted Euclidiean distances
      .. $silhouetteData  The silhouette data of the select partition
      .. $AveSilData      The average silhouette value for each partition
      .. $Moran           The Moran index for each groups for each partitions
      .. $lonlat          Logical parameter if your coordinates are in longitude latitude format or not.

If you ran mapclust like this: classif <- mapclust(...)
    See summary of the silhouette data by running:
    > summary(classif$silhouetteData)

You can cut the dendrogram for a new dlim with: NewCut <- mapclust_cut_tree(...)
    See documentation by running: \n    > ?mapclust_cut_tree

documentation: \n> ?mapclust
        ")

    return(out)

  }

print.mapclust_cl <- function(x, ...){
  cat("
      Divisive hierarchical clustering | mapclust_cl

call: ")
  print(x$call)
  cat(sep = "","\npartition:
      number of cluster: ",length(unique(x$cluster)),"\n      Average_Sil = ",arrondi(mean(x$silhouetteData[,3]),2),", WSS = ",arrondi(sum(x$InertiaWithin[[which(x$dlim == x$cutdlim)]]),2),", cutdlim = ",x$cutdlim,
      "\nplot:
      .. $dendrogram      The global dendrogramm
      .. $cuttree         The cut dendrogramm
      .. $map             The map of the select partition
      .. $AveSilPlot      The average silhouette width plot (for each partition)
      .. $WSSPlot         The Within Sum of Square plot (for each partition)
      .. $silhouette      The silhouette of the select partition
value:
      .. $X               The x-coordinates data you used
      .. $Y               The y-coordinates data you used
      .. $var             The regionalized variable.s data you used
      .. $label           The label vector you used
      .. $density         The estimate density based on var. Equal to var if you used a unidimensionnal density variable
      .. $cluster         The vector of cluster of the select partition
      .. $fullhist        The composition cluster for each observation
      .. $hist            The composition cluster without duplicates (matches to leaf of the dendrogramm)
      .. $dlim            The vector of the different limit distances
      .. $cutdlim         The select dlim for the cut
      .. $InertiaWithin   The Sum of Square for each patch of each partition
      .. $DiMatrix        The matrix of Weighted Euclidiean distances
      .. $silhouetteData  The silhouette data of the select partition
      .. $AveSilData      The average silhouette value for each partition
      .. $Moran           The Moran index for each groups for each partitions
      .. $lonlat          Logical parameter if your coordinates are in longitude latitude format or not.

If you ran mapclust like this: classif <- mapclust(...)
    See summary of the silhouette data by running:
    > summary(classif$silhouetteData)

You can cut the dendrogram for a new dlim with: NewCut <- mapclust_cut_tree(...)
    See documentation by running: \n    > ?mapclust_cut_tree

documentation: \n> ?mapclust
      ")
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

"mapclust_cut_tree" <- function(classification, nb_grp = NA, dlim = NA)
{
  #===================================================================================
  # CUT DENDROGRAMM
  #
  # Author : A. Coulon (Universite de Nantes)
  # Last update : 09 October 2018
  #
  # Arguments:
  # classification	  The object return by mapclust
  # dlim              The value of dlim where you want to cut the dendrogramm
  # Plabel            Logical parameter for activate or not the print of labels on dendrogram
  #
  #===================================================================================
  lonlat <- classification$lonlat
  if(!is.na(dlim) && !is.na(nb_grp)){
    if(nb_grp - arrondi(nb_grp) == 0){
      message("you have entered both parameters.")
      message("The parameter dlim will be ignore!")
      dlim <- NA
    }else{
      message("you have entered a non-integer group number: ",nb_grp)
      message("The parameter nb_grp will be ignore!")
      nb_grp <- NA
    }
  }
  if(!nb_grp - arrondi(nb_grp) == 0 && is.na(dlim)){
    dlim <- nb_grp
    nb_grp <- NA
    message("you have entered a non-integer group number: ",dlim,". We consider this value as a dlim.")
  }
  if(!is.na(dlim)){
    if(dlim < 0){
      message("Dlim must be positive We change it to the closest to 0.")
    }
    if(dlim > max(classification$dlim)){
      message("Your dlim is too high.")
    }
    #On recupere les variables de l'objet classification et des autres arguments
    diff <- classification$dlim-dlim
    if((TRUE %in% (diff <= 0))){
      #return index of the first value <= 0
      index <- which.max(diff <= 0)
    }else{
      index <- which.min(diff)
    }

    cutdlim <- classification$dlim[index]
  }else{
    if(!is.na(nb_grp)){
      V_NB_GRP <- c()
      index <- NULL
      for(i in 1:length(classification$dlim)){
        V_NB_GRP <- c(V_NB_GRP,length(unique(classification$historique[,i])))
        if(nb_grp == length(unique(classification$historique[,i]))){
          cutdlim <- classification$dlim[i]
          index <- i
        }
      }

      if(is.null(index)){
        message("There is not this number of group in the partition. We cut to the closest. possible number:")
        print(unique(V_NB_GRP)[-1])
        diff <- V_NB_GRP - nb_grp
        index <- which.min(abs(diff))
        cutdlim <- classification$dlim[index]
      }
    }else{
      stop("You must enter a dlim or a number of group in order to cut the dendrogram.")
    }
  }
  N <- length(unique(classification$fullhist[,index+1]))
  change <- 0
  while (N == 1){
    change <- 1
    index <- index + 1
    cutdlim <- classification$dlim[index]
    N <- length(unique(classification$fullhist[,index+1]))
  }
  if (change == 1){message("You can't create only one cluster. We use the next dlim which create at least two cluster : ",cutdlim)}

  if(dim(as.data.frame(classification$var))[2] > 1){
    Size <- classification$var[,1]
  }else{
    Size <- classification$var
  }
  cat("Running ..")
  #classification$dendrogramme$theme correspond au metadonnees du dendrogramme (concernant le theme)
  #c'est une liste et si il y a du texte alors la liste est de longueur deux

  Plabel = classification$Plabel

  multi <- dim(as.data.frame(classification$var))[2] > 1
  tmp <- data.frame(lon=classification$X,lat=classification$Y)

  lineType <- ifelse(Size < 0,"2, 5","1")
  #leaflet
  if(!lonlat){
    tmp$lon <- (tmp$lon - mean(tmp$lon)) / max(tmp$lon, tmp$lat)
    tmp$lat <- (tmp$lat - mean(tmp$lat)) / max(tmp$lon, tmp$lat)
  }
  lat <- tmp$lat
  lng <- tmp$lon
  size <- 10*(0.1+(abs(Size)/max(abs(Size))))+5
  color <- as.factor(classification$fullhist[,index+1])
  df <- data.frame(
    lat,
    lng,
    size,
    color
  )

  if(multi){df$size <- rep(7,length(df$size))}
  # Make up some random levels. (TODO: Better example)
  couleur <- c("slateblue2", "blue", "chartreuse1","darkolivegreen2",
               "#222222", "darkmagenta", "mediumpurple4",
               "magenta","cyan3", "plum4","darkorange1", "pink1", "red",
               "lightgoldenrod", "violetred3", "brown4", "tan1","yellow", "gray50",
               "#967272","deepskyblue4","salmon", "cyan","palegreen", "green3")
  #plot(1:25,1:25,pch =15,cex=3,col=couleur)
  df$category <- as.factor(classification$fullhist[,index+1])
  factpal <- colorFactor(couleur, df$category)
  coul <- factpal(classification$historique[,index])

  if(lonlat){
    esri <- structure(c("Esri.WorldStreetMap", "Esri.DeLorme", "Esri.WorldTopoMap",
                        "Esri.WorldImagery", "Esri.WorldTerrain", "Esri.WorldShadedRelief",
                        "Esri.WorldPhysical", "Esri.OceanBasemap", "Esri.NatGeoWorldMap",
                        "Esri.WorldGrayCanvas"),
                      .Names = c("Esri.WorldStreetMap", "Esri.DeLorme", "Esri.WorldTopoMap",
                                 "Esri.WorldImagery", "Esri.WorldTerrain", "Esri.WorldShadedRelief",
                                 "Esri.WorldPhysical", "Esri.OceanBasemap", "Esri.NatGeoWorldMap",
                                 "Esri.WorldGrayCanvas"))
    m = leaflet(df)
    for (provider in esri) {
      m <- m %>% addProviderTiles(provider, group = provider)
    }

    m <- m %>% addCircleMarkers(lng = ~lng, lat = ~lat, radius = ~size, fillColor =  ~factpal(category),
                                fill = TRUE, fillOpacity = 0.8, opacity=1, stroke = TRUE,weight = 1.5, color="black", dashArray=lineType
    ) %>%
      addLegend("bottomleft", pal = factpal, values = as.factor(classification$fullhist[,index+1]),
                title = "Cluster",
                opacity = 1
      ) %>%
      addMiniMap(toggleDisplay = TRUE,
                 position = "bottomright"
      ) %>%
      addLayersControl(baseGroups = names(esri),
                       options = layersControlOptions(collapsed = TRUE),
                       position = "topright"
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
    if(!multi){
      m <- m %>%
        addLegendCustom(colors = c("black", "black", "black"),
                        labels = c(arrondi(min(abs(Size)),2),
                                   arrondi(min(abs(Size)) + (max(abs(Size))-min(abs(Size)))/2  ,2),
                                   arrondi(max(abs(Size)),2)),
                        sizes = c(10*(0.1+( min(abs(Size)) / max(abs(Size))))+5,
                                  10*(0.1+((min(abs(Size)) + (max(abs(Size))-min(abs(Size)))/2) /max(abs(Size))))+5,
                                  10*(0.1+( max(abs(Size)) / max(abs(Size))))+5),
                        title = "Size"
        )
    }
    if(min(Size)<0){
      m <- m %>%
        addLegendCustom(colors = c("white", "white"), labels = c("line: positive","dotted-line: negative"),
                        sizes = 1, title="Sign"
        )
    }
  }else{
    m <- leaflet(df)
    m <- m %>% addCircleMarkers(lng = ~lng, lat = ~lat, radius = ~size, fillColor =  ~factpal(category),
                                fill = TRUE, fillOpacity = 0.8, opacity=1, stroke = TRUE,weight = 1.5, color="black", dashArray=lineType
    ) %>%
      addLegend("bottomleft", pal = factpal, values = as.factor(classification$fullhist[,index+1]),
                title = "Cluster",
                opacity = 1
      )
  }

  cluster <- as.numeric(classification$fullhist[,index+1])

  cuttree <- classification$dendrogram_ggplot
  cuttree <- cuttree + geom_hline(yintercept = cutdlim + 0.03 , col = "black" , lwd = 1)
  print(Plabel)
  if(Plabel){
    cuttree <- cuttree + suppressWarnings(theme(axis.text.x = element_text(color = coul,
                                                          hjust = 1, vjust = 0.25, size = 9, angle = 90)))
  }else{
    cuttree <- cuttree + theme(axis.text.x = element_blank())
  }

  sil <- silhouette(cluster,classification$DiMatrix)
  coull <- factpal(sort(cluster))
  plot(sil,col = coull)
  silPlot <- recordPlot()
  plot.new()
  # Nouvelle page
  grid.newpage()
  # Creer la mise en page : nrow = 2, ncol = 2
  pushViewport(viewport(layout = grid.layout(2, 2)))
  # Une fonction pour definir une region dans la mise en page
  define_region <- function(row, col){
    viewport(layout.pos.row = row, layout.pos.col = col)
  }
  # Arranger les graphiques dans la fenetre
  print(cuttree, vp = define_region(1:2, 1))
  print(classification$AveSilPlot, vp = define_region(1, 2))
  print(classification$WSSPlot, vp = define_region(2, 2))
  print(m)

    out <- structure(list(X=classification$X,
                        Y=classification$Y,
                        var=classification$var,
                        label=classification$label,
                        density=classification$density,
                        cluster=cluster,
                        fullhist=classification$fullhist,
                        dlim=classification$dlim,
                        historique=classification$historique,
                        Plabel=classification$Plabel,
                        InertiaWithin=classification$InertiaWithin,
                        Moran=classification$Moran,
                        dendrogram=classification$dendrogram,
                        dendrogram_ggplot=classification$dendrogram_ggplot,
                        DiMatrix=classification$DiMatrix,
                        cuttree=cuttree,
                        map=m,
                        silhouette=silPlot,
                        silhouetteData=sil,
                        AveSilData=classification$AveSilData,
                        AveSilPlot=classification$AveSilPlot,
                        WSSPlot=classification$WSSPlot,
                        cutdlim=cutdlim,
                        lonlat=lonlat,
                        call=match.call()),class=c("mapclust_cl","list"))
    cat(". Done")
    cat("\n\n\n
      Divisive hierarchical clustering

call: ")
    print(out$call)
    cat(sep = "","\npartition:
      number of cluster: ",length(unique(out$cluster)),"\n      Average_Sil = ",arrondi(mean(out$silhouetteData[,3]),2),", WSS = ",arrondi(sum(out$InertiaWithin[[which(out$dlim == out$cutdlim)]]),2),", cutdlim = ",out$cutdlim,
        "\nplot:
      .. $dendrogram      The global dendrogramm
      .. $cuttree         The cut dendrogramm
      .. $map             The map of the select partition
      .. $AveSilPlot      The average silhouette width plot (for each partition)
      .. $WSSPlot         The Within Sum of Square plot (for each partition)
      .. $silhouette      The silhouette of the select partition
value:
      .. $X               The x-coordinates data you used
      .. $Y               The y-coordinates data you used
      .. $var             The regionalized variable(s) data you used
      .. $label           The label vector you used
      .. $density         The estimate density based on var. Equal to var if you used a unidimensionnal density variable
      .. $cluster         The vector of cluster of the select partition
      .. $fullhist        The composition cluster for each observation
      .. $hist            The composition cluster without duplicates (matches to leaf of the dendrogramm)
      .. $dlim            The vector of the different limit distances
      .. $cutdlim         The select dlim for the cut
      .. $InertiaWithin   The Sum of Square for each patch of each partition
      .. $DiMatrix        The matrix of Weighted Euclidiean distances
      .. $silhouetteData  The silhouette data of the select partition
      .. $AveSilData      The average silhouette value for each partition
      .. $Moran           The Moran index for each groups for each partitions
      .. $lonlat          Logical parameter if your coordinates are in longitude latitude format or not.

If you ran mapclust like this: classif <- mapclust(...)
    See summary of the data by running:
    > summary(classif$silhouetteData)

You can cut the dendrogram for a new dlim with: NewCut <- mapclust_cut_tree(...)
    See documentation by running:
    > ?mapclust_cut_tree

documentation: \n> ?mapclust
        ")
  return(out)
}


summary.mapclust_cl <- function(object, ...){
  plot(object$cuttree)
  cat("\n mapclust: Divise hierarchical clustering")
  if(dim(as.data.frame(object$var))[2] == 1){
    cat("\n -------------------------------------
  dlim\t  n\t WSS\t sil\tMoran
 -------------------------------------
")
  }else{
    cat("\n -------------------------------
  dlim\t  n\t WSS\t sil
 -------------------------------
")
  }
  j <- 1
  for (i in 1:length(object$dlim)){
    dlim <- object$dlim[i]
    number_grp <- length(unique(object$historique[,i]))
    if(number_grp == 1){
      sil <- NA
    }else{
      sil <- object$AveSilData[j]
      j <- j+1
    }
    WSS <- sum(object$InertiaWithin[[i]])
    if(dim(as.data.frame(object$var))[2] == 1){
      Moran <- mean(na.omit(unlist(object$Moran[[i]][1,])))
    }


    cat(" ",dlim,"\t")
    cat(" ",number_grp,"\t")

    cat(arrondi(WSS,2),"\t")
    cat(arrondi(sil,2),"\t")
    aste <- FALSE
    if(dim(as.data.frame(object$var))[2] == 1){
      cat(arrondi(Moran,2))
      if(anyNA(object$Moran[[i]][1,]) == TRUE){
        nbTOT <- length(object$Moran[[i]][1,])
        nbNA <- length(is.na(object$Moran[[i]][1,])[is.na(object$Moran[[i]][1,])==TRUE])
        aste <- TRUE
        cat("\t* ",nbNA,"NA /",nbTOT)
      }
    }
    cat("\n")
  }
  if(dim(as.data.frame(object$var))[2] == 1){
    cat(" -------------------------------------")
  }else{
    cat(" -------------------------------")
  }
  if(aste){
    cat("\n* Some groups in this partition have less than 3 observations and we can't calculate the Moran index for these groups. There is at least one of these groups in the partition. In this case the mean value is calculated only with the groups with more than 2 observations.")
  }
  nb <- length(unique(object$cluster))
  tot <- length(object$X)
  X_mean <- tapply(object$X,object$cluster,mean)
  Y_mean <- tapply(object$Y,object$cluster,mean)
  size <- c()
  cat("\n\n")
  cat("The selected partition contains",nb,"clusters among",tot,"observations.\ncluster:\t")
  for (k in 1:nb){
    size <- c(size,length(object$cluster[object$cluster == k]))
    cat(k,"\t")
  }
  cat("\nSize:\t\t")
  for(y in 1:nb){
    cat(size[y],"\t")
  }
  cat("\nX_mean:\t\t")
  for(y in 1:nb){
    cat(arrondi(X_mean[y],2),"\t")
  }
  cat("\nY_mean:\t\t")
  for(y in 1:nb){
    cat(arrondi(Y_mean[y],2),"\t")
  }
  if(dim(as.data.frame(object$var))[2] == 1){
    var_mean <- tapply(object$var,object$cluster,mean)
    cat("\nvar_mean:\t")
    for(y in 1:nb){
      cat(arrondi(var_mean[y],2),"\t")
    }
  }else{
    for(i in 1:dim(as.data.frame(object$var))[2]){
      var_mean <- tapply(object$var[,i],object$cluster,mean)
      cat("\nvar",i,"_mean:\t",sep="")
      for(y in 1:nb){
        cat(arrondi(var_mean[y],2),"\t")
      }
    }
  }
  object$map
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
                   abs( coef(fit)[2]*x_values[i] - y_values[i] + coef(fit)[1] ) / sqrt( coef(fit)[2]^2 + 1^2 ))
  }

  # Max distance point
  x_max_dist <- x_values[which.max(distances)]
  y_max_dist <- y_values[which.max(distances)]
  return(distances)
  #return(c(x_max_dist, y_max_dist))
}

