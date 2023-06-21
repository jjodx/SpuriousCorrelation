library(correlation) # for multilevel correlations
library(MASS) # for mvnorm

resMuCo <- data.frame(p=rep(NA,Nsim),R=rep(NA,Nsim),Sig=rep(NA,Nsim),NoDifButCorrel=rep(NA,Nsim))
D <- mvrnorm(n = PopSize, rep(0, 2), Sigma, empirical = FALSE)
# creating the two subgroups
SubN = round(PopSize/2)
D[1:SubN,]<-D[1:SubN,]+matrix(rep(SDdistSub,len=SubN*2), nrow = SubN)
Z <- c(rep(0,SubN),rep(1,PopSize-SubN))
Add <- data.frame(D1=D[,1],D2=D[,2],Z=factor(Z))
# computing multilevel correlation
Rmulti <- correlation(Add, multilevel = TRUE)
resMuCo[i,"R"] <-Rmulti$r #Normalized coefficient
resMuCo[i,"p"] <- Rmulti$p #p-value
resMuCo[i,"Sig"]<-resMuCo[i,"p"]<0.05