library(mixtools)
setwd("D:/Data/FACS/ISE single reps x3 cbind 3 columns/")
folder = getwd()
dirs = list.dirs(path = folder, full.names = TRUE, recursive = TRUE)

for (j in 1:length(dirs)){
  setwd(dirs[j])
  print(dirs[j])
  folder = getwd()
  txt_list = list.files(path = folder, pattern ="*cbind.txt")

for (i in 1:length(txt_list)){
  file_name = txt_list[i]
  
  data1=read.table(file_name, header=F)   # reads the data into R; edit the filename 
  #data1[data1==0] <- NA
  
  lambda = c(0.9, 0.1)
  mu = c(8.6, 8.7)
  sigma = c(0.5, 0.6)
  
  nmix = 2
  
  #replicate 1
  temp = data1$V1
  temp[temp==0] <- NA
  x=na.omit(temp)
  nmixv1=normalmixEM(log(x), lambda, mu, sigma, k=nmix, maxit=10000, maxrestarts = 100, arbmean = TRUE, arbvar = TRUE)
  mu1 = nmixv1$mu[1]
  lambda1 = nmixv1$lambda[1]
  mu2 = nmixv1$mu[2] 
  
  mu1_weighted = mu1*lambda1
  mu2_weighted = mu2*(1-lambda1)
  mu_weighted_V1 = mu1_weighted + mu2_weighted
  
  #replicate 2
  temp = data1$V2
  temp[temp==0] <- NA
  x=na.omit(temp)
  nmixv1=normalmixEM(log(x), lambda, mu, sigma, k=nmix, maxit=10000, maxrestarts = 100, arbmean = TRUE, arbvar = TRUE)
  mu1 = nmixv1$mu[1]
  lambda1 = nmixv1$lambda[1]
  mu2 = nmixv1$mu[2] 
  
  mu1_weighted = mu1*lambda1
  mu2_weighted = mu2*(1-lambda1)
  mu_weighted_V2 = mu1_weighted + mu2_weighted
  
  #replicate 3
  temp = data1$V3
  temp[temp==0] <- NA
  x=na.omit(temp)
  nmixv1=normalmixEM(log(x), lambda, mu, sigma, k=nmix, maxit=10000, maxrestarts = 100, arbmean = TRUE, arbvar = TRUE)
  mu1 = nmixv1$mu[1]
  lambda1 = nmixv1$lambda[1]
  mu2 = nmixv1$mu[2] 
  
  mu1_weighted = mu1*lambda1
  mu2_weighted = mu2*(1-lambda1)
  mu_weighted_V3 = mu1_weighted + mu2_weighted
  
  #read data as vector
  dat1_v <- as.vector(log(data1[['V1']]))
  dat2_v <- as.vector(log(data1[['V2']]))
  dat3_v <- as.vector(log(data1[['V3']]))
  
  #remove NA
  dat1_v <- na.omit(dat1_v)
  dat2_v <- na.omit(dat2_v)
  dat3_v <- na.omit(dat3_v)
  
  
  #normalises for weighted
  dat1_v <- dat1_v-mu_weighted_V1
  dat2_v <- dat2_v-mu_weighted_V2
  dat3_v <- dat3_v-mu_weighted_V3
  
  #inverts log transformation
  dat1_v <- exp(dat1_v)
  dat2_v <- exp(dat2_v)
  dat3_v <- exp(dat3_v)
  
  data_final <- c(dat1_v, dat2_v, dat3_v)
  
  file_name = paste("mu_weighted_normalised_",txt_list[i], sep="")
  write.table(data_final, file_name, append = TRUE, row.names=FALSE, col.names=FALSE, sep=",")
}
}
  ######################################################

