# function to obtain p-value from correlations
# library(Hmisc) # for rcorr
# library(MASS) # for mvrnorm
# library(ggplot2) # for plotting
# library(ggpubr) # for ggarrange
# library(fGarch) # for skewed distribution
# library(WRS2) # for wincor and pbcor
# function to plot the data
# function to plot the data
report <- function(ListOut){
  dat <- ListOut[[1]]
  Scatter1 <- ListOut[[2]]
  POP <- ListOut[[3]]
  TXT <- ListOut[[4]]
  TiTxt <- TXT[[1]]
  InfoTxt <- TXT[[2]]
  Rval <- TXT[[3]]
  CountSig <- sum(dat$Sig)
  Psig <- CountSig/length(dat$Sig)
  dat$Sig = factor(1-dat$Sig)
  Fcolor <- c("orange","lightblue")
  Rsig<-0
  if (sum(dat$Sig==1)<length(dat$Sig)){Rsig <- mean(dat$R[dat$Sig==0])}
  
  Rall <- mean(dat$R)
  if (Rsig>Rval){Rpos <- c("left","right")}else{Rpos <- c("right","left")}
  
  PP1 <- ggplot(dat, aes(x=R, fill=Sig)) + 
    geom_histogram(aes(x=R,fill=Sig),binwidth=.05,boundary = 0,color="black") +
    scale_fill_manual(values = Fcolor) +
    coord_cartesian(xlim = c(-1, 1))+
    geom_vline(xintercept=Rsig, linetype="dashed", color = "darkorange", size = 1) +
    geom_label(aes(x = Rsig, y = Inf, label = "Sig."), 
               hjust = Rpos[1], 
               vjust = "top", 
               colour = "darkorange", 
               label.size = NA, 
               size = 3,
               fill = NA)+
    geom_vline(xintercept=Rval, linetype="dashed", color = "darkblue", size = 1) +
    geom_label(aes(x = Rval, y=Inf, label = "Actual"),
               hjust = Rpos[2], 
               vjust = "top", 
               colour = "darkblue", 
               label.size = NA, 
               size = 3,
               fill = NA)+
    theme(legend.position='none') +
    scale_y_continuous(expand = c(0.15,0)) +
    labs(title="",x="correlation coefficient" , y="Count")
  
  
  
  PP2 <- ggplot(dat, aes(x=p, fill=Sig)) + 
    geom_histogram(aes(x=p,fill=Sig),binwidth=.05,boundary = 0.05,color="black") +
    scale_fill_manual(values = Fcolor) +
    geom_hline(yintercept=CountSig, linetype="dashed", color = "black") +
    geom_label(aes(x = 0, y = CountSig, label = "H0"), 
               hjust = "left", 
               vjust = "bottom", 
               colour = "black", 
               label.size = NA, 
               size = 3,
               fill = NA)+
    coord_cartesian(xlim = c(0, 1))+
    theme(legend.position='none')+
    scale_y_continuous(expand = c(0.1,0)) +
    labs(title="",x="p-value" , y="Count")
  
  d=data.frame(x1=c(0.33,0.33), x2=c(0.66,0.66), y1=c(0,Psig), y2=c(Psig,1), t=c('a','b'))
  PP3 <- ggplot() +  
    geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=t), 
              color="black") +
    scale_fill_manual(values = Fcolor, 
                      name="",
                      labels=c("p<=0.05", "p>0.05")) +
    geom_hline(yintercept=0.05, linetype="dashed", color = "red") +
    geom_label(aes(x = 0.1, y = 0.05, label = "Type 1 \nError"), 
               hjust = "center", 
               vjust = "bottom", 
               colour = "red", 
               label.size = NA, 
               size = 3,
               fill = NA)+
    geom_hline(yintercept=0.8, linetype="dashed", color = "black") +
    geom_label(aes(x = 0.1, y = 0.8, label = "Type 2 \nError"), 
               hjust = "center", 
               vjust = "center", 
               colour = "black", 
               label.size = NA, 
               size = 3,
               fill = NA)+
    geom_hline(yintercept=0, color = "black") +
    coord_cartesian(xlim = c(0, 1),ylim = c(0, 1))+
    labs(title="",y="\n \n \nProportion of simulations")+
    scale_y_continuous(expand = c(0,0)) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  
  
  cols <- c("royalblue3", "orangered2" )
  PP4 <- ggplot(Scatter1, aes(x=D1, y=D2, color=Zf)) +
    geom_point()+
    scale_color_manual(name = "Population(s)", 
                       values = cols,
                       labels = POP)+
    labs(title = "",x="Random variable 1" , y="Random variable 2")
  
  PGA <- ggarrange(PP4,PP1,PP2,PP3,ncol = 2, nrow = 2,labels = c("A","B","C","D"))
  annotate_figure(PGA,
                  top = text_grob(TiTxt, face = "bold", size = 14),
                  bottom = text_grob(InfoTxt)
  )
}




# function to simulate the correlations
# PopulationSize = size of the simulated sample. default N=20
# Nsim: number of simulations performed. Ideally 10000. For speed; default = 1000
# Rval: we are going to simulate two variables. Rval is the theoretical influence of these variables? Defaut = 0
# Problem can be: none, outlier, subpop(ulations) or skewed
# SDdist: case Problem = none -> useless
#         case Problem = outlier -> distance of outlier from population mean (unit is standard deviation)
#         case Problem = subpop -> distance between two subpopulation means (unit is standard deviation)
SimulateCorr <- function(PopSize=20,Nsim=1000,Rval=0,Problem="none",SDdist=3,Skewness=2,Solution="no"){
  # simulations of p-value distribution in the absence of correlations.
  Sigma <- matrix(c(1,Rval,Rval,1),2,2)
  
  dProblem <- switch(Problem, 
                     "none" = 0,
                     "outlier" = 1,
                     "subpop" = 2,
                     "skewed" = 3
  )
  dSolution <- switch(Solution,
                      "no" = 0,
                      "yes" = 1
  )
  POP <- switch(Problem, 
                "none" = c("Main"),
                "outlier" = c("Main","outlier"),
                "subpop" = c("Pop1","Pop2"),
                "skewed" = c("Main")
  )
  res <- data.frame(p=rep(NA,Nsim),R=rep(NA,Nsim),Sig=rep(NA,Nsim))
  Z <- rep(0,PopSize)
  for (i in c(1:Nsim)){
    if (dProblem <3){
      D <- mvrnorm(n = PopSize, rep(0, 2), Sigma, empirical = FALSE)
      if (dProblem == 1) {
        mMeans=apply(D,2,mean)
        mStd=apply(D,2,sd)
        D[PopSize,]<-mMeans+SDdist*mStd #rep(SDdist,2)
        Z <- c(rep(0,PopSize-1),1)
      }
      if (dProblem == 2) {
        SubN = round(PopSize/2)
        D[1:SubN,]<-D[1:SubN,]+matrix(rep(SDdist,len=SubN*2), nrow = SubN)
        Z <- c(rep(0,SubN),rep(1,PopSize-SubN))
        D1 <- D[,1]
        D2 <- D[,2]
        Add <- data.frame(D1,D2,Z)
      }
    }
    else { # in case of skewness, populations are supposed to be correlated
      X <- rsnorm(PopSize, mean = 0, sd = 1, xi = Skewness)
      Y <- rsnorm(PopSize, mean = 0, sd = 1, xi = 1)
      D = cbind(X,Y)
    }
    if (i==1){
      D1 <- D[,1]
      D2 <- D[,2]
      Zf = factor(Z)
      Scatter1 <- data.frame(D1,D2,Zf)
    }
    if (dSolution==0){
      R<-rcorr(D[,1],D[,2])
      res[i,"R"] <- R$r[1,2]
      res[i,"p"]<-R$P[1,2]
      res[i,"Sig"]<-res[i,"p"]<0.05
    } 
    else {
      if (dProblem <= 1){ # outlier
        R<-wincor(D[,1],D[,2],0.2)
        res[i,"R"] <- R$cor
        res[i,"p"]<-R$p.value
        res[i,"Sig"]<-res[i,"p"]<0.05
      }
      else { # two subpopulations
        Rrr <- lm(D2 ~ D1+Z, Add)
        SumRrr <- summary(Rrr, correlation = TRUE)
        res[i,"R"] <-SumRrr$coefficients[2,1] #Normalized coefficient
        res[i,"p"] <- SumRrr$coefficients[2,4] #p-value
        res[i,"Sig"]<-res[i,"p"]<0.05
      }
    }
  }
  Txt <- writeTXT(PopSize=PopSize,Nsim=Nsim,Rval=Rval,Problem=Problem,SDdist=SDdist,Skewness=Skewness,Solution=Solution)
  return(list(res,Scatter1,POP,Txt))
}

writeTXT <- function(PopSize=20,Nsim=1000,Rval=0,Problem="none",SDdist=3,Skewness=2,Solution="no"){
  Prob <- switch(Problem, 
                 "none" = NULL,
                 "outlier" = " with an outlier",
                 "subpop" = " with two subpopulations",
                 "skewed" = " with a skewed distribution")
  dProblem <- switch(Problem, 
                     "none" = 0,
                     "outlier" = 1,
                     "subpop" = 2,
                     "skewed" = 3)
  dSolution <- switch(Solution,
                      "no" = 0,
                      "yes" = 1)
  if (dSolution==1){
    if (dProblem<=1){
      TxtSol <- ".\nType: winsorized correlation"}
    else{
      TxtSol <- ".\nType: regression with subpopulation as factor"}
  }
  else{
    TxtSol <- ".\nType: Pearson correlation"}
  
  TitleTxt <- paste0("Simulations of correlation (R= ",Rval,
                     ")", Prob,TxtSol)
  ParamTxT <- paste0("# of simulations: ", Nsim,"; Population Size:",PopSize,
                     ";")
  return(list(TitleTxt,ParamTxT,Rval))
}