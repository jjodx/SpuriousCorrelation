library(MASS) # for mvnorm
library(Hmisc) # for rcorr
library('groundhog')
pkgs=c('foreach','doParallel')
groundhog.library(pkgs,'2022-07-01')
library(correlation)

# defining default parameters
PopSize=15 # size of the population/sample
Nsim=10000 # number of simulated populations/samples. 100 to go fast / should be 10000
Rval=0 # for simulating when there is no correlations between the measured parameters
Rval50=0.5 # for generating correlated measured parameters within a population
BinSize <- 0.01 # for histograms of p-values
SDdistSub <- 3
Sigma <- matrix(c(1,0,0,1),2,2)# covariance matrix



res <- data.frame(p=rep(NA,Nsim),R=rep(NA,Nsim),Sig=rep(NA,Nsim)) # pre-allocating space for the dataframe

starttime <- Sys.time ()
for (i in c(1:Nsim)){
  # drawing 15 elements for each random variable (x and y) with a covariance matrix equal to Sigma
  # the first column of D corresponds to the samples of the x parameter (e.g. the height) and the second to the samples of the y parameter (e.g. working memory capacity)
  D <- mvrnorm(n = PopSize, rep(0, 2), Sigma, empirical = FALSE) 
  SubN = round(PopSize/2)
  D[1:SubN,]<-D[1:SubN,]+matrix(rep(SDdistSub,len=SubN*2), nrow = SubN)
  Z <- c(rep(-1,SubN),rep(1,PopSize-SubN))
  Add <- data.frame(D1=D[,1],D2=D[,2],Z=factor(Z))
  Rmulti <- correlation::correlation(Add, multilevel = TRUE)
  res[i,"R"] <- Rmulti$r # storing the correlation value for future use
  res[i,"p"]<-Rmulti$p# storing the p-value for future use
  res[i,"Sig"]<-Rmulti$p<0.05 # tag p-values smaller than the type II error threshold (here 0.05)
}

stoptime = Sys.time () - starttime
print("for loop")
print(stoptime)

# inspired from http://datacolada.org/102
Ncores<- parallel::detectCores()

#get the core processors ready to go
registerDoParallel(cores=round(Ncores/2)) #I leave 1 free to do things while R runs
starttime <- Sys.time ()
#Actual loop
res = foreach(i = 1:Nsim, .combine='rbind') %dopar% {
  D <- MASS::mvrnorm(n = PopSize, rep(0, 2), Sigma, empirical = FALSE) 
  SubN = round(PopSize/2)
  D[1:SubN,]<-D[1:SubN,]+matrix(rep(SDdistSub,len=SubN*2), nrow = SubN)
  Z <- c(rep(-1,SubN),rep(1,PopSize-SubN))
  Add <- data.frame(D1=D[,1],D2=D[,2],Z=factor(Z))
  Rmulti <- correlation::correlation(Add, multilevel = TRUE)
  c(Rmulti$r,Rmulti$p,Rmulti$p<0.05)
}
colnames(res)<-c("R","P","Sig")


stoptime = Sys.time () - starttime
print("parallel loop")
print(stoptime)