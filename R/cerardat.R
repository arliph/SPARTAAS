press <- function(linear.model) {
  # calculate the predictive residuals
  pr <- residuals(linear.model)/(1-lm.influence(linear.model)$hat)
  # calculate the PRESS
  PRESS <- mean(pr^2)
  return(PRESS)
}

cerardat = function(df, row.sup, date, nf = NULL, confidence = 0.95, graph = T){

  ## test
  #
  # df is data.frame
  if(!is.data.frame(df)){stop("df is not a data frame.")}
  #
  # row.sup is vector
  if(!is.vector(row.sup)){stop("row.sup is not a vector.")}
  #
  # row.sup is interger
  #if(!is.integer(row.sup)){stop("row.sup is not integer.")}
  #
  # date is vector
  if(!is.vector(date)){stop("date is not a vector.")}
  #
  # date is na or interger
  #if(!is.integer(date)){stop("date is not integer.")}
  #
  # nf is integerdate
  #if(!is.integer(nf)){stop("nf is not integer.")}
  #
  # max(row.sup) < 280 size df
  if(max(row.sup) > dim(df)[1]){stop("row.sup must contain integers between 1 and the number of lines in df.")}
  # min(row.sup) > 0
  if(min(row.sup) <= 0){stop("row.sup must contain integers between 1 and the number of lines in df.")}
  #
  # length(date) = nrow(df)
  if(!length(date) == nrow(df)){stop("date must have the same number of observations as the number of lines in df. Complete with NA if necessary.")}
  #
  # nf > 0
  if(!is.null(nf)){
    if(nf < 0){stop("nf must be positive.")}
  }

  #
  # date is na or interger
  #
  # confidence < 1
  if(confidence > 0.99){stop("Confidence must be between 0 and 0.99.")}
  # confidence >= 0
  if(confidence < 0){stop("Confidence must be between 0 and 0.99.")}
  #
  ##

  lwr = 0
  upr = 0

  #vector for ref data (every not in row.sup)
  row.ref = which(!(1:length(df[,1]) %in% row.sup))

  #todo: while convergence

  #rm row(ens) < 5
  #si row < 5
  if(sum(rowSums(df)<5)!=0){
    warning(paste0("The sums of rows ",capture.output(cat(row.names(df)[which(rowSums(df)<5)]))," are less than 5. They were suppressed from the analysis."))
    #retire index de la row dans row.sup
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
  }



  #nettoyage des GT ie GT<5
  GT_rm_ref = which(colSums(df[row.ref,])<5)

  if(sum(colSums(df[row.ref,])<5)!=0){
    warning(paste0("The sums of columns (GT) ",capture.output(cat(names(df)[which(colSums(df[row.ref,])<5)]))," are less than 5. They were suppressed from the analysis."))
  }
  #?????????????????
  #GT_rm_sup = which(colSums(df[-c(GT_rm_ref),row.sup])<5)

  #data reference
  if(sum(GT_rm_ref)==0){
    data_ref = df[row.ref,]
  }else{
    data_ref = df[row.ref,-c(GT_rm_ref)]
  }

  if(sum(GT_rm_ref)==0){
    data_sup = df[row.sup,]
  }else{
    data_sup = df[row.sup,-c(GT_rm_ref)]
  }

  #DOF Degrees of freedom
  DOF = min(ncol(data_ref)-1,nrow(data_ref)-1)

  ## test
  #
  # nf <= DOF

  if(is.null(nf)){
    #AFC avec  tous les axes
    #ca.NMI.init <- ade4::dudi.coa(data_ref,scannf=FALSE,nf=DOF)

    #nb axe inertie > 60%
    #k=min(which(ade4::inertia.dudi(ca.NMI.init, row.inertia=F, col.inertia=T)$tot.inertia[,3]>60))
    k=cerardat_estim_nf(df, row.sup, date)$nf[1]
  }else{
    k=nf
  }

  #Correspondance analysis
  ca.NMI <- ade4::dudi.coa(data_ref,scannf=FALSE,nf=k)


  res.CA = FactoMineR::CA(data_ref, ncp = k, graph = graph)

  #passer par FactoMineR::CA()
  obs_ca_eval = res.CA$row$cos2
  #obs_ca_eval = t(rowSums(res.CA$col$cos2))

  #df with coord ens (col) de reference + date monnaie
  DATA_REF = cbind(ca.NMI$li, date=date[row.ref])
  #here

  #genrate col coord for sup data
  data_sup_proj_col = ade4::suprow(ca.NMI,data_sup)$lisup
  #df with coord GT (col) supplementaire + date monnaie
  DATA_SUP = cbind(data_sup_proj_col, date=date[row.sup])

  #final data
  DATA_TOT = rbind(DATA_REF, DATA_SUP)
  DATA_TOT = DATA_TOT[sort(c(row.ref,row.sup)),]

  #linear regression
  #DATA_REF_lm = MASS::steapAIC(lm(date~.,data=DATA_REF, na.action=na.omit),trace=0)
  DATA_REF_lm = lm(date~.,data=DATA_REF, na.action=na.omit)

  R_adj = arrondi(summary(DATA_REF_lm)$adj.r.squared,3)
  R_sq = arrondi(summary(DATA_REF_lm)$r.squared,3)
  sigma = arrondi(summary(DATA_REF_lm)$sigma,2)

  if(graph == 1){
    plot(DATA_REF_lm,which=c(1))
    plot(DATA_REF_lm,which=c(2))
    MASS::boxcox(DATA_REF_lm,data=DATA_REF)
  }


  ## residus student
  mod1.rstudent<-as.data.frame(which(abs(rstudent(DATA_REF_lm))>2))
  dimnames(mod1.rstudent)[[2]]<-c("abs(rstudent)")

  # H0 Normalite des residus : Shapiro-Wilks Test
  Shapiro=shapiro.test(rstudent(DATA_REF_lm))
  if (Shapiro$p.value<0.05) {warning("The Shapiro-Wilks test indicates a problem with the normality of the residuals.",sep="")}
  # H0 : pas Autocorrelation D-W  test
  DW=lmtest::dwtest(date~.,data=DATA_REF,alternative = "two.sided")
  if (DW$p.value<0.05) {warning("The Durbin-Watson test indicates a first order autocorrelation problem.",sep="")}
  # H0 : Homoscedasticite, B-P test
  BP=lmtest::bptest(date~.,data=DATA_REF)
  if (BP$p.value<0.05) {warning("The Breusch-Pagan test indicates a heteroskedasticity problem.",sep="")}

  #eval model
  mod1.diagGl = data.frame(
    cbind(
      R_adj, R_sq, sigma,
      arrondi(shapiro.test(rstudent(DATA_REF_lm))$p.value,3),
      arrondi(lmtest::dwtest(date~.,data=DATA_REF, alternative = "two.sided")$p.value, 3),
      arrondi(lmtest::bptest(date~.,data=DATA_REF)$p.value,3)
    )
  )
  dimnames(mod1.diagGl)[[2]]<-c("R2_aj", "R2", "sigma", "Shapiro p-value", "D-W p-value", "B-P p-value")

  #prediction date ensemble
  predict_obj_row = predict(DATA_REF_lm, newdata=DATA_TOT, se.fit=TRUE, interval="confidence", level=confidence)

  #prediction date GT
  dimnames(ca.NMI$co)[[2]] = dimnames(ca.NMI$li)[[2]]
  predict_obj_col = predict(DATA_REF_lm, newdata=ca.NMI$co, se.fit=TRUE, interval="confidence", level=confidence)

  ##### 95% des dates des GT
  date_gt = arrondi(predict_obj_col$fit,0)
  dimnames(date_gt)[[2]]<-c("Fit_dateEv","lower","upper")

  #matrix proportion GT
  if(sum(GT_rm_ref)==0){
    cont = df[sort(c(row.ref,row.sup)),]
  }else{
    cont = df[sort(c(row.ref,row.sup)),-c(GT_rm_ref)]
  }

  cont.gt = cont

  for(j in 1:nrow(cont))
  {
    cont.gt[j,]<-cont[j,]/as.matrix(rowSums(cont))[j]
  }

  median.norMix <- function(x) nor1mix::qnorMix(1/2,x)

  median_tab = as.data.frame(matrix(0, nrow(cont),3))
  dimnames(median_tab)[[1]] = dimnames(cont)[[1]]
  dimnames(median_tab)[[2]] = c("Median_dateAc","lower","upper")

  date_gt_sd = arrondi(cbind(as.data.frame(predict_obj_col$fit)[,1], as.data.frame(predict_obj_col$se.fit)),1)

  for(i in 1:nrow(cont))
  {
    ex = nor1mix::norMix(mu = date_gt_sd[,1],w=unlist(as.vector(cont.gt[i,])),sigma= date_gt_sd[,2])
    median_tab[i,] = cbind(median.norMix(ex), nor1mix::qnorMix(0.025,ex), nor1mix::qnorMix(0.975,ex))
  }

  date_predict = arrondi(cbind(
    date[sort(c(row.ref,row.sup))],
    predict_obj_row$fit,
    median_tab
  ),0)
  dimnames(date_predict)[[1]]<-dimnames(predict_obj_row$fit)[[1]]
  dimnames(date_predict)[[2]]<-c("date","Fit_dateEv","lower_Ev","upper_Ev","Median_dateAc","lower_Ac","upper_Ac")

  tmp = date_predict[row.ref,]
  tmp = tmp[!is.na(tmp$date),]

  df = data.frame(
    names = rownames(tmp),
    date = tmp$date,
    lwr = tmp$lower_Ev,
    upr = tmp$upper_Ev,
    fit = tmp$Fit_dateEv
  )

  #generation du graph
  check_ref = ggplot(df) +
    geom_point(aes(x = names, y = date),colour="red",shape=3, size=3) +
    geom_errorbar( aes(x=names, ymin=lwr, ymax=upr), width=.15, colour="#3f0b18", alpha=0.9, linewidth =1)+
    ylab("Date (year)") + xlab("Site") + ggtitle("Comparison of estimated dates (black) with actual dates (red)") +
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
    )

  tmp2 = date_predict[row.sup,]
  tmp2 = tmp2[!is.na(tmp2$date),]



  df2 = data.frame(
    names = rownames(tmp2),
    date = tmp2$date,
    lwr = tmp2$lower_Ev,
    upr = tmp2$upper_Ev,
    fit = tmp2$Fit_dateEv
  )

  #generation du graph
  check_sup = ggplot(df2) +
    geom_point(aes(x = names, y = date),colour="red",shape=3, size=3) +
    geom_errorbar( aes(x=names, ymin=lwr, ymax=upr), width=.15, colour="#3f0b18", alpha=0.9, linewidth=1)+
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
    )

  print(summary(DATA_REF_lm))

  return(
    structure(
      list(
        prediction          = date_predict,
        date_gt             = date_gt,
        lm                  = DATA_REF_lm,
        k                   = k,
        predict_obj_row     = predict_obj_row,
        predict_obj_col     = predict_obj_col,
        cont_gt             = cont.gt,
        statistical_summary = mod1.diagGl,
        obs_ca_eval         = obs_ca_eval,
        check_ref           = check_ref,
        check_sup           = check_sup,
        Shapiro_Wilks       = Shapiro,
        Durbin_Watson       = DW,
        Breusch_Pagan       = BP,
        row.sup             = row.sup,
        call                = match.call()
      ),
      class = c("cerardat_obj","list")
    )
  )

}

cerardat_estim_nf <- function(df, row.sup, date){

  # df is data.frame
  if(!is.data.frame(df)){stop("df is not a data frame.")}
  #
  # row.sup is vector
  if(!is.vector(row.sup)){stop("row.sup is not a vector.")}
  #
  # row.sup is interger
  #if(!is.integer(row.sup)){stop("row.sup is not integer.")}
  #
  # date is vector
  if(!is.vector(date)){stop("date is not a vector.")}
  #
  # date is na or interger
  #if(!is.integer(date)){stop("date is not integer.")}
  #
  # nf is integerdate
  #if(!is.integer(nf)){stop("nf is not integer.")}
  #
  # max(row.sup) < 280 size df
  if(max(row.sup) > dim(df)[1]){stop("row.sup must contain integers between 1 and the number of lines in df.")}
  # min(row.sup) > 0
  if(min(row.sup) <= 0){stop("row.sup must contain integers between 1 and the number of lines in df.")}
  #
  # length(date) = 280
  if(!length(date) == nrow(df)){stop("date must have the same number of observations as the number of lines in df. Complete with NA if necessary.")}
  #

  nf = 0

  #vector for ref data (every not in col.sup)
  row.ref = which(!(1:length(df[,1]) %in% row.sup))

  #rm row(ens) < 5
  #si row < 5
  if(sum(rowSums(df)<5)!=0){
    warning(paste0("The sums of rows ",capture.output(cat(row.names(df)[which(rowSums(df)<5)]))," are less than 5. They were suppressed from the analysis."))
    #retire index de la row dans row.sup
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
  }



  #nettoyage des GT ie GT<5
  GT_rm_ref = which(colSums(df[row.ref,])<5)

  if(sum(colSums(df[row.ref,])<5)!=0){
    warning(paste0("The sums of columns (GT) ",capture.output(cat(names(df)[which(colSums(df[row.ref,])<5)]))," are less than 5. They were suppressed from the analysis."))
  }
  #?????????????????
  #GT_rm_sup = which(colSums(df[-c(GT_rm_ref),row.sup])<5)

  #data reference
  if(sum(GT_rm_ref)==0){
    data_ref = df[row.ref,]
  }else{
    data_ref = df[row.ref,-c(GT_rm_ref)]
  }

  if(sum(GT_rm_ref)==0){
    data_sup = df[row.sup,]
  }else{
    data_sup = df[row.sup,-c(GT_rm_ref)]
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


  SW <- c()
  DW <- c()
  BP <- c()

  for(i in 1:max){
    if(i == 1){
      formula <- paste0(formula,' Axis',i)
    }else{
      formula <- paste0(formula,' + Axis',i)
    }
    #model <- caret::train(as.formula(formula), data = DATA_REF, method = "lm", trControl = ctrl, na.action=na.omit)
    #RMSE <- c(RMSE,model$results$RMSE)
    #R_sq <- c(R_sq,model$results$Rsquared)

    lm = lm(as.formula(formula),data = DATA_REF)
    R_sq <- c(R_sq,summary(lm)$adj.r.squared)
    MSE <- c(MSE,mean(lm$residuals^2))
    PRESS <- c(PRESS,press(lm))

    SW = c(SW,arrondi(shapiro.test(rstudent(lm))$p.value,3))
    DW = c(DW,arrondi(lmtest::dwtest(as.formula(formula),data = DATA_REF, alternative = "two.sided")$p.value,3))
    BP = c(BP,arrondi(lmtest::bptest(as.formula(formula),data = DATA_REF)$p.value,3))

  }

  data2 = data.frame(
    nf = 1:max,
    MSE = MSE,
    PRESS = PRESS,
    R_sq = R_sq
  )

  data3 = data.frame(
    nf = rep(1:max,3),
    p.value=c(SW,DW,BP),
    test=c(rep("Shapiro-Wilks",max),rep("Durbin-Watson",max),rep("Breusch-Pagan",max))
  )

  pTest = ggplot(data3) +
    ylab("p.value") + xlab("Number of component in lm()") +
    ggtitle("hypothesis tests p.value") +
    geom_point(shape=1,aes(x = nf,y = p.value, col=test)) +
    geom_line(aes(x = nf,y = p.value, col=test)) +
    geom_line(aes(x = nf,y = 0.05)) +
    theme_classic()

  pMSE = ggplot(data2) +
    ylab("MSE") + xlab("Number of component in lm()") +
    ggtitle("Mean Squared Error (MSE)") +
    geom_point(shape=1,aes(x = nf,y = MSE)) +
    #geom_point(shape=16,cex=1.1,aes(x = which(MSE == min(MSE)),y = min(MSE))) +
    geom_line(aes(x = nf,y = MSE)) +
    #geom_ribbon(alpha=0.5,aes(x = nf,y = MSE,ymin=MSE-MSE_sd, ymax=MSE+MSE_sd)) +
    #geom_errorbar(aes(x = nf, ymin=MSE-MSE_sd, ymax=MSE+MSE_sd), width=.2,
    #              position=position_dodge(.9)) +
    #geom_vline(xintercept = which(MSE == min(MSE)), linetype = "dotted",
    #           color = "grey", size=1) +
    theme_classic()

  pPRESS = ggplot(data2) +
    ylab("PRESS") + xlab("Number of component in lm()") +
    ggtitle("PRediction Error Sum Of Squares (PRESS)") +
    geom_point(shape=1,aes(x = nf,y = PRESS)) +
    geom_point(shape=16,cex=1.1,aes(x = which(PRESS == min(PRESS)),y = min(PRESS))) +
    geom_line(aes(x = nf,y = PRESS)) +
    scale_y_log10(labels = scientific) +
    annotation_logticks(sides="l") +
    geom_vline(xintercept = which(PRESS == min(PRESS)), linetype = "dotted",
               color = "grey", size=1) +
    theme_classic()

  min_R = min(data2$R_sq)
  if(min(data2$R_sq) > 0)min_R= 0
  if(min(data2$R_sq) > 0.25)min_R= 0.25
  if(min(data2$R_sq) > 0.5)min_R= 0.5
  if(min(data2$R_sq) > 0.75)min_R= 0.75

  pR_sq = ggplot(data2) +
    ylab("adj.R_squared") + xlab("Number of component in lm()") +
    ggtitle("adj.R_squared") + ylim(min_R,1) +
    geom_point(shape=1,aes(x = nf,y = R_sq)) +
    geom_line(aes(x = nf,y = R_sq)) +
    theme_classic()

  return(
    structure(
      list(
        nf= which(data2$PRESS == min(data2$PRESS)),
        MSE = pMSE,
        PRESS = pPRESS,
        hypothesis = pTest,
        adj.R_sq = pR_sq,
        data = data2,
        data_stat = data3
      ),
      class = c("list")
    )
  )

}

print.cerardat_obj <- function(x, ...){
  ## test
  #
  # is cerardat
  if(class(x)[1] != "cerardat_obj"){stop("x is not a cerardat object.")}

  cat("cerardat
call: ")
  print(x$call)
  cat("\n")
  print(x$statistical_summary)
  cat("\n")
  head(x$date_predict)
}

plot.cerardat_obj = function(x,
                             which = NULL,
                             col1=rgb(0.93,0.23,0.23,0.5),
                             col2="black",
                             xlim=NULL,
                             ylim=NULL,...){
  #ylim=c(0,0.03)
  ## test
  #
  # is cerardat
  if(class(x)[1] != "cerardat_obj"){stop("x is not a cerardat object.")}



  if(is.null(xlim)){
    ecart = max(x$date_gt)-min(x$date_gt)
    xlim = c(min(x$date_gt-ecart*0.07),max(x$date_gt+ecart*0.07))
  }


  AUTO_TICK = trunc((xlim[2]-xlim[1])/100)

  if(is.null(which)){
    which = 1:nrow(x$cont_gt)
  }


  if(length(which) > 1){
    pb <- utils::txtProgressBar(min = 0, max = length(which), style = 3)
  }

  ens_date_sd = arrondi(cbind(
    Pred=data.frame(x$predict_obj_row$fit)[,1],
    data.frame(E_T=x$predict_obj_row$se.fit)
  ),1)

  GT_date_sd = arrondi(cbind(
    Pred=data.frame(x$predict_obj_col$fit)[,1],
    data.frame(E_T=x$predict_obj_col$se.fit)
  ),1)

  if(is.null(ylim)){
    tmp_ = c()
    for(i in 1:ncol(x$cont_gt) )
    {
      date_accumulation <- nor1mix::norMix(mu = GT_date_sd[,1], w = unlist(as.vector(x$cont_gt[i,])), sigma= GT_date_sd[,2])
      date_accumulation_density = nor1mix::dnorMixL(date_accumulation, xlim=xlim)

      tmp_ = c(tmp_,max(date_accumulation_density$y))
    }
    ylim=c(0,max(tmp_))
  }

  row.ref = which(!(1:length(ens_date_sd$Pred) %in% x$row.sup))

  for(i in which)
  {
    date_accumulation <- nor1mix::norMix(mu = GT_date_sd[,1], w = unlist(as.vector(x$cont_gt[i,])), sigma= GT_date_sd[,2])
    date_accumulation_density = nor1mix::dnorMixL(date_accumulation, xlim=xlim)

    date_evenement = nor1mix::norMix(ens_date_sd[i,1],sigma=ens_date_sd[i,2])
    date_evenement_density = nor1mix::dnorMixL(date_evenement,xlim=xlim)

    ymax = max(ylim) + 0.003
    date_accumulation_density$y[date_accumulation_density$y > ymax] = ymax
    date_evenement_density$y[date_evenement_density$y > ymax] = ymax


    date_accumulation_density$y[1] = 0
    date_evenement_density$y[1] = 0

    date_accumulation_density$y[length(date_accumulation_density$y)] = 0
    date_evenement_density$y[length(date_accumulation_density$y)] = 0
    sub=""


    if(i %in% x$row.sup){
      sub=""
    }else{
      sub=paste("Quality of row representation (cos2) in correspondence analysis:",arrondi(sum(x$obs_ca_eval[which(row.ref == i),1:x$k]),2))
    }


    plot(date_accumulation_density,xlim=xlim, ylim=ylim,xlab="Date",col="black",
         ylab=dimnames(x$cont_gt)[[1]][i],type= "l",xaxt="n",
         main="model dateEv (red) and dateAc (black)",
         sub=sub,...)


    graphics::axis(side=1,col="black",at=pretty(seq(xlim[1],xlim[2]),AUTO_TICK))

    graphics::polygon(date_accumulation_density,col=col2,border=col2)
    graphics::polygon(date_evenement_density,col=col1,border=NA)

    if(length(which) > 1)
      utils::setTxtProgressBar(pb, i)
  }

}

extract_results = function(cerardat,
                            width = 480,
                            height = 480,
                            path = "figures",
                            col1 = rgb(0.93,0.23,0.23,0.5),
                            col2 = "black",
                            xlim = NULL,
                            ylim = NULL){
  #ylim=c(0,0.03)
  ## test
  #
  # is cerardat
  if(class(cerardat)[1] != "cerardat_obj"){stop("cerardat is not a cerardat object.")}

  if(is.null(xlim)){
    ecart = max(cerardat$date_gt)-min(cerardat$date_gt)
    xlim = c(min(cerardat$date_gt-ecart*0.07),max(cerardat$date_gt+ecart*0.07))
  }

  AUTO_TICK = trunc((xlim[2]-xlim[1])/100)

  pb <- utils::txtProgressBar(min = 0, max = length(cerardat$cont_gt[1,]), style = 3)

  if(!dir.exists(path))
    dir.create(path)

  ens_date_sd = arrondi(cbind(
    Pred=data.frame(cerardat$predict_obj_row$fit)[,1],
    data.frame(E_T=cerardat$predict_obj_row$se.fit)
  ),1)

  GT_date_sd = arrondi(cbind(
    Pred = data.frame(cerardat$predict_obj_col$fit)[,1],
    data.frame(E_T=cerardat$predict_obj_col$se.fit)
  ),1)


  if(is.null(ylim)){
    tmp_ = c()
    for(i in 1:nrow(cerardat$cont_gt) )
    {
      date_accumulation <- nor1mix::norMix(mu = GT_date_sd[,1], w = unlist(as.vector(cerardat$cont_gt[i,])), sigma= GT_date_sd[,2])
      date_accumulation_density = nor1mix::dnorMixL(date_accumulation, xlim=xlim)

      tmp_ = c(tmp_,max(date_accumulation_density$y))
    }
    ylim=c(0,max(tmp_))
  }

  for(i in 1:nrow(cerardat$cont_gt) )
  {
    date_accumulation = nor1mix::norMix(mu=GT_date_sd[,1],w=unlist(as.vector(cerardat$cont_gt[i,])),sigma=GT_date_sd[,2])
    date_accumulation_density = nor1mix::dnorMixL(date_accumulation,xlim=xlim)

    date_evenement = nor1mix::norMix(ens_date_sd[i,1],sigma=ens_date_sd[i,2])
    date_evenement_density = nor1mix::dnorMixL(date_evenement,xlim=xlim)

    ymax = max(ylim) + 0.003
    date_accumulation_density$y[date_accumulation_density$y > ymax] = ymax
    date_evenement_density$y[date_evenement_density$y > ymax] = ymax

    date_accumulation_density$y[1] = 0
    date_evenement_density$y[1] = 0

    date_accumulation_density$y[length(date_accumulation_density$y)] = 0
    date_evenement_density$y[length(date_accumulation_density$y)] = 0

    grDevices::tiff(file=path.expand(paste(path,"/",
                    dimnames(cerardat$cont_gt)[[1]][i],".jpeg",sep="")),
                    width = width, height = height,units = "px", pointsize = 12)


    plot(date_accumulation_density,xlim=xlim, ylim=ylim,xlab="Date",col="black",
           ylab=dimnames(cerardat$cont_gt)[[1]][i],type= "l",xaxt="n",
           main="model dateEv (red) and dateAc (black)",
           sub=paste("Quality of row representation (cos2) in correspondence analysis: ",arrondi(cerardat$obs_ca_eval[i],2)))


    graphics::axis(side=1,col="black",at=pretty(seq(xlim[1],xlim[2]),AUTO_TICK))

    graphics::polygon(date_accumulation_density,col=col2,border=col2)
    graphics::polygon(date_evenement_density,col=col1,border=NA)

    grDevices::graphics.off()

    utils::setTxtProgressBar(pb,i)
  }
  print(paste0(getwd(),"/",path))
}

