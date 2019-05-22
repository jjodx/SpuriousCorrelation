library(Hmisc) # for rcorr
library(MASS) # for mvnorm
library(psycho) # for bayes_cor.test
set.seed(123)# to make sure everybody gets the same results
Nsim=100 # number of simulations
PopSize = 15 #number os elements in each sample
Sigma <- matrix(c(1,0,0,1),2,2)# covariance matrix
SDdist = 3
res <- data.frame(p=rep(NA,Nsim),R=rep(NA,Nsim),Sig=rep(NA,Nsim)) # pre-allocating space for the dataframe
for (i in c(1:Nsim)){
  # drawing 15 elements for each random variable (x and y) with a covariance matrix equal to Sigma
  # the first column of D corresponds to the samples of the x variable and the second to the samples of the y variable
  D <- mvrnorm(n = PopSize, rep(0, 2), Sigma, empirical = FALSE)
  mMeans<-apply(D,2,mean)
  mStd<-apply(D,2,sd)
  D[PopSize,]<-mMeans+SDdist*mStd
  Z <- c(rep(0,PopSize-1),1)
  # DframeBayes <- data.frame(D1=D[,1],D2=D[,2])
  bf <- bayes_cor.test(D[,1],D[,2])
  #bf <- correlationBF(x=D[,1],y=D[,2])
  R<-rcorr(D[,1],D[,2]) # correlation between the samples of x and y 
  res[i,"R"] <- R$r[1,2] # storing the correlation value for future use
  res[i,"p"]<-R$P[1,2]# storing the p-value for future use
  res[i,"Sig"]<-res[i,"p"]<0.05 # tag p-values smaller than the type II error threshold (here 0.05)
  res[i,"Bf"]<-bf$values$bf#bf@bayesFactor$bf #bf$values$bf[1,2]
}
Fcolor <- c("lightblue","orange")
ggplot(res, aes(x=Bf, fill=Sig)) + 
  geom_histogram(aes(x=Bf,fill=Sig),binwidth=.05,boundary = 0.05,color="black")+
  scale_fill_manual(values = Fcolor)
