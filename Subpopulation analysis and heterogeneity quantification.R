library(mixtools)
setwd("D:/Data/FACS/ISE single reps x3 combined 1 column mu_weighted subtracted renamed/")
folder = getwd()
dirs = list.dirs(path = folder, full.names = TRUE, recursive = TRUE)

lambda = c(0.9, 0.1)
mu = c(8.5, 9)
sigma = c(8, 8.5)

nmix2=function(data, nset){
  bic=matrix(0, 7, nset)
  for (k in 1:nset){
    temp=data[, k]
    temp[temp==0] <- NA
    x=na.omit(temp)
    nx=length(x)
    mu1=mean(log(x))
    sigma1=sd(log(x))
    mixc1=dnorm(log(x), mean=mu1, sd=sigma1)
    mixc1.loglik=sum(log(mixc1))
    bic[1,k]=-2*mixc1.loglik+(3*1-1)*log(nx)
    mixc2=normalmixEM(log(x), lambda, mu, sigma, k=2, maxit=10000, maxrestarts = 100, arbmean = TRUE, arbvar = TRUE)
    mixc3=normalmixEM(log(x), lambda, mu, sigma, k=3, maxit=10000, maxrestarts = 100, arbmean = TRUE, arbvar = TRUE)
    mixc4=normalmixEM(log(x), lambda, mu, sigma, k=4, maxit=10000, maxrestarts = 100, arbmean = TRUE, arbvar = TRUE)
    mixc5=normalmixEM(log(x), lambda, mu, sigma, k=5, maxit=10000, maxrestarts = 100, arbmean = TRUE, arbvar = TRUE)
    mixc6=normalmixEM(log(x), lambda, mu, sigma, k=6, maxit=10000, maxrestarts = 100, arbmean = TRUE, arbvar = TRUE)
    mixc7=normalmixEM(log(x), lambda, mu, sigma, k=7, maxit=10000, maxrestarts = 100, arbmean = TRUE, arbvar = TRUE)
    bic[2,k]=-2*mixc2$loglik+(3*2-1)*log(nx)
    bic[3,k]=-2*mixc3$loglik+(3*3-1)*log(nx)
    bic[4,k]=-2*mixc4$loglik+(3*4-1)*log(nx)
    bic[5,k]=-2*mixc5$loglik+(3*5-1)*log(nx)
    bic[6,k]=-2*mixc6$loglik+(3*6-1)*log(nx)
    bic[7,k]=-2*mixc7$loglik+(3*7-1)*log(nx)
  }
  return(bic)
}

nmixden=function(xmin, xmax, nmix, lambda, mu, sigma){
  xx=seq(from=xmin, to=xmax, length=500)
  xxd=vector(mode="numeric", length=500)
  for(i in 1:500){
    for(k in 1:nmix){
      xxd[i]=xxd[i]+lambda[k]*dnorm(xx[i], mean=mu[k], sd=sigma[k])
    }
  }
  results=data.frame(xx, xxd)
  return(results)
}


ncount3=function(data, nset){
  res=0
  for (k in 1:nset){
    temp=data[, k]
    x=na.omit(temp)
    res[k]=length(x)
  }
  res
}

dstat1=function(data, nmix, lam, mean, stdev){
  data = data[!is.na(data)]
  m0=mean(data)
  sd0=sd(data)
  lam=c(-1, lam)
  mean=c(m0, mean)
  stdev=c(sd0, stdev)
  dstat=0
  for(i in 1:(nmix+1)){
    for(j in 1:(nmix+1)){
      t1=(lam[i]*lam[j])/sqrt(2*pi*(stdev[i]^2 + stdev[j]^2))
      t2=exp(-0.5*(mean[i]-mean[j])^2/(stdev[i]^2+stdev[j]^2))
      dstat=dstat+t1*t2
    }}
  return(dstat)
}

boot.bc=function(theta, bsdata, ci.level=0.95){
  bsdata=sort(bsdata)
  tail.p=0.5*(1-ci.level)
  B=length(bsdata)
  #calculate upper end-point
  ct=0
  for(i in 1:B){
    if(bsdata[i] < theta) ct=ct+1
  }
  b=qnorm(ct/(B+1))
  #print(b)
  Q=(B+1)*pnorm(2*b+qnorm(1-tail.p))
  Q=as.integer(floor(Q))
  ci.u=bsdata[Q]
  Q=(B+1)*pnorm(2*b+qnorm(tail.p))
  Q=as.integer(floor(Q))
  if(Q==0) Q=Q+1
  ci.l=bsdata[Q]
  bci=c(ci.l, ci.u)
  return(bci)
}

boot.ise13=function(data, nmix, nboot){
  #ci's based on bootstrap percentile method
  data = data[!is.na(data)]
  bdat=NULL
  boot.res=NULL
  boot.lambda=NULL
  boot.mu=NULL
  boot.sigma=NULL
  x=log(data)
  n=length(x)
  res.dat=normalmixEM(x, k=nmix, epsilon=1.0e-8, maxit=100000)
  ise.dat=dstat1(x, nmix, res.dat$lambda, res.dat$mu, res.dat$sigma)
  for(j in 1:nboot){
    w = rmultinom(n, size = 1, prob = res.dat$lambda)
    bdat = sapply(1:n, function(i) rnorm(1, mean = res.dat$mu[(w[, 
                                                                 i] == 1)], sd = res.dat$sigma[w[, i] == 1]))
    #bdat=sample(x, size=n, replace=TRUE)
    resb=try(normalmixEM(bdat, k=nmix, lambda=res.dat$lambda, mu=res.dat$mu, sigma=res.dat$sigma, 
                         epsilon=1.0e-8, maxit=100000), silent=TRUE)
    if (class(resb) == "try-error" || resb$restarts != 0) {
      j = j - 1            }
    else {
      boot.res[j]=dstat1(bdat, nmix, resb$lambda, resb$mu, resb$sigma)
      boot.lambda=cbind(boot.lambda, resb$lambda)
      boot.mu=cbind(boot.mu, resb$mu)
      boot.sigma=cbind(boot.sigma, resb$sigma)
    }}
  ise.ci=NULL
  ise.ci=quantile(boot.res, probs=c(0.025, 0.975), na.rm = TRUE, type=6) #alternatively ise.ci=quantile(boot.res, probs=c(0.025, 0.975), na.rm = TRUE, type=6)
  lambda.ci=matrix(0, nrow=2, ncol=nmix)
  mu.ci=matrix(0, nrow=2, ncol=nmix)
  sigma.ci=matrix(0, nrow=2, ncol=nmix)
  for(k in 1:nmix){
    lambda.ci[,k]=quantile(boot.lambda[k,], probs=c(0.025, 0.975), type=6)
    mu.ci[,k]=quantile(boot.mu[k,], probs=c(0.025, 0.975), type=6)
    sigma.ci[,k]=quantile(boot.sigma[k,], probs=c(0.025, 0.975), type=6)
  }
  param.ci=cbind(lambda.ci, mu.ci, sigma.ci)
  row.names(param.ci)=c("ci_lower", "ci_upper")
  results=list(res.dat[2:4], ise.dat, ise.ci, boot.res, param.ci)
  return(results)
}

#################################

#################################################
for (j in 2:2){
  setwd(dirs[j])
  print(dirs[j])
  folder = getwd()
  txt_list = list.files(path = folder, pattern =".txt") #pattern ="^mu_weighted"
  
  header = c("file_name","nobs","sub_pops","lambda1","CI_lambda1_low","CI_lambda1_high","lamda1_CI_min","lambda1_CI_max","lambda2", "CI_lambda2_low","CI_lambda2_high","lambda2_CI_min","lambda2_CI_max","mu1","CI_mu1_low","CI_mu1_high","mu1_CI_min","mu1_CI_max","mu2","CI_mu2_low","CI_mu2_high","mu2_CI_min","mu2_CI_max","sigma1","CI_sigma1_low","CI_sigma1_high","sigma1_CI_min","sigma1_CI_max","sigma2","CI_sigma2_low","CI_sigma2_high","sigma2_CI_min","sigma2_CI_max","ISE","ISE_CI_low","ISE_CI_high","ISE_CI_min","ISE_CI_max")
  FF <- as.matrix(t(header))
  write.table(FF, "file_save.txt", append = TRUE, col.names=FALSE, sep=",")

for (i in 1:length(txt_list)){
  file_name = txt_list[i]
  
  data1=read.table(file_name, header=F)   # reads the data into R; edit the filename 
  data1[data1==0] <- NA
  # appropriately before using this command
  
  nsamps=1    # sets the number of samples in the dataset, which can be edited accordingly
  
  nobs3=ncount3(data1, nsamps)   # calculates the number of observations for each sample in the dataset
  
  nobs3    #  prints out the counts
  
  con.res=nmix2(data1, nsamps)   # fits the Normal mixture model to the log data for each sample using  1-7 Normals
  
  con.res    # prints an array of BIC values with the rows corresponding to the number of Normals and the
  # columns to each sample
  subpops = which.min(con.res)
  
  #con.sums=rowSums(con.res)   # adds up the BIC values accross rows
  
  #con.sums    # prints out the row sums; the best fitting mixture model (ie. number of Normal pdf's it contains) 
  #  for all the samples corresponds to the entry with the smallest value; 
  # assign  this value to "nmix" in the next command
  
  #could do nmix=subpops here
  nmix=subpops       # the value of nmix can be edited here
  
  
  
  ################################################################
  
  boot1=boot.ise13(data1, nmix=subpops, nboot=1000)
  
  
  #  there are 8 samples in the example code here;  the next set of commands fits the mixture
  #  model having nmix Normal pdf's to each of the 8 samples and stores the results in 
  #  in a seperate object for each sample;  these objects are
  #  called nmixv1 - nmixv8;
  #  edit the code accordingly for a different number of samples
  
  
  v1=na.omit(data1[,1])
  # removes the missing values from the data 
  #                          in the first column of the data object
  length(v1)  # returns the number of actual observations in the first sample
  nmixv1=normalmixEM(log(v1), lambda, mu, sigma, k=nmix, maxit=10000, maxrestarts = 100, arbmean = TRUE, arbvar = TRUE)  # fits the mixture model with 
  m0=mean(log(v1))
  sd0=sd(log(v1))
  #                             nmix pdf's and stores the results in object nmixv1
  #nmixv1$lambda  # prints the lambda estimates
  #nmixv1$mu         # prints the mu estimates
  #nmixv1$sigma     # prints the sigma estimates
  
  
  ###################################################
  
  
  #  the next set of commands calculate the ISE between the fitted single and mixture
  #  Normal models for each of the 8 samples using the function dstat;  
  #  edit the code accordingly for a different number of samples
  
  lambda1 = boot1[1][[1]]$lambda[1]  
  lambda1_CI_low = as.numeric(boot1[[5]])[1]
  lambda1_CI_high = as.numeric(boot1[[5]])[2]
  lambda1_CI_min = lambda1-lambda1_CI_low
  lambda1_CI_max = lambda1_CI_high-lambda1
  
  lambda2 = boot1[1][[1]]$lambda[2] 
  lambda2_CI_low = as.numeric(boot1[[5]])[3]
  lambda2_CI_high = as.numeric(boot1[[5]])[4]
  lambda2_CI_min = lambda2-lambda2_CI_low
  lambda2_CI_max = lambda2_CI_high-lambda2
  
  mu1 = boot1[1][[1]]$mu[1]
  mu1_CI_low = as.numeric(boot1[[5]])[5]
  mu1_CI_high = as.numeric(boot1[[5]])[6]
  mu1_CI_min = mu1-mu1_CI_low
  mu1_CI_max = mu1_CI_high-mu1
  
  mu2 = boot1[1][[1]]$mu[2]
  mu2_CI_low = as.numeric(boot1[[5]])[7]
  mu2_CI_high = as.numeric(boot1[[5]])[8]
  mu2_CI_min = mu2-mu2_CI_low
  mu2_CI_max = mu2_CI_high-mu2
  
  sigma1 = boot1[1][[1]]$sigma[1]
  sigma1_CI_low = as.numeric(boot1[[5]])[9]
  sigma1_CI_high = as.numeric(boot1[[5]])[10]
  sigma1_CI_min = sigma1-sigma1_CI_low
  sigma1_CI_max = sigma1_CI_high-sigma1
  
  sigma2 = boot1[1][[1]]$sigma[2]
  sigma2_CI_low = as.numeric(boot1[[5]])[11]
  sigma2_CI_high = as.numeric(boot1[[5]])[12]
  sigma2_CI_min = sigma2-sigma2_CI_low
  sigma2_CI_max = sigma2_CI_high-sigma2
  
  ISE = as.numeric(boot1[2])
  ISE_CI_low = as.numeric(boot1[[3]])[1]
  ISE_CI_high = as.numeric(boot1[[3]])[2]
  ISE_CI_min = ISE-ISE_CI_low
  ISE_CI_max = ISE_CI_high-ISE
  
  output<- data.frame(file_name,nobs3, subpops, lambda1, lambda1_CI_low, lambda1_CI_high, lambda1_CI_min, lambda1_CI_max, lambda2, lambda2_CI_low, lambda2_CI_high, lambda2_CI_min, lambda2_CI_max, mu1, mu1_CI_low, mu1_CI_high, mu1_CI_min, mu1_CI_max, mu2, mu2_CI_low, mu2_CI_high, mu2_CI_min, mu2_CI_max, sigma1, sigma1_CI_low, sigma1_CI_high, sigma1_CI_min, sigma1_CI_max, sigma2, sigma2_CI_low, sigma2_CI_high, sigma2_CI_min, sigma2_CI_max, ISE, ISE_CI_low, ISE_CI_high, ISE_CI_min, ISE_CI_max)
  data <- assign(file_name, output, .GlobalEnv)
  
  out_dir = paste(dirs[j],"/file_save.txt", sep="")
  write.table(data, out_dir, append = TRUE, col.names=FALSE, sep=",")
  ######################################################
  
  
  # the next sets of commands produce plots of the results
  
  
  #######################################################
  
  #  run this tiff command to store the subsequent plots to tiff files
  file_name = unlist(strsplit(file_name, split='.', fixed=TRUE))[1]
  tmp = paste(file_name, "tif", sep = ".")
  
  tiff(filename = tmp,
       width = 1024, height = 1024, units = "px", pointsize = 32,
       compression = "none",
       bg = "white", res = NA, family = "", restoreConsole = TRUE,
       type = "windows")
  
  
  #######################################################
  
  #  the argument xlim sets the range of x-values and ylim the range of y-values for the plot;
  #  argument main sets the title for the plot
  x_start = log(min(v1))-2
  x_end = log(max(v1))+2
  
  opar <- par(lwd=2)
  
  hist(log(v1), freq=F, xlim=c(x_start, x_end), ylim=c(0,1.2), main= tmp, xlab="log(fluorescence intensity)", ylab="density")
  par(opar)
  
  #  the first two argumets in the following call to the function nmixden are the min amd max x-values set for the plot;
  #  nmix specifies the number of mixtures;
  #  the last three arguments specify the parameter estimates for the Normal mixture model contained in the object 
  #  nmixv1 here and which were calculated previously;
  #  object plv1 is used to add the full mixture density to the hitogram plot 
  
  #plots theoretical normal against fitted mixture distribution
  xxc=seq(from=x_start, to=x_end, length.out=500)
  dxxc=nmixv1$lambda[1]*dnorm(xxc, mean=nmixv1$mu[1], sd=nmixv1$sigma[1])+nmixv1$lambda[2]*dnorm(xxc, mean=nmixv1$mu[2], sd=nmixv1$sigma[2])
  
  lines(xxc, dxxc, lwd=5) #mixture model plot
  
  dxxs=dnorm(xxc, mean=m0, sd=sd0)
  
  lines(xxc, dxxs, col="blue", lwd=5) #theoretical normal distribution
  
  # the following code adds the first Normal density, scaled by the relevant estimate of lambda to the plot
  
  xx1=seq(from=x_start, to=x_end, length=500)
  xxd1=nmixv1$lambda[1]*dnorm(xx1, nmixv1$mu[1], nmixv1$sigma[1])
  lines(xx1, xxd1, col="red", lwd=3, lty=1)
  
  #  the following code adds the second Normal density, scaled by the relevant estimate of lambda, to the plot
  
  xx2=seq(from=x_start, to=x_end, length=500)
  xxd2=nmixv1$lambda[2]*dnorm(xx2, nmixv1$mu[2], nmixv1$sigma[2])
  lines(xx2, xxd2, col="green", lwd=3, lty=1)
  
  dev.off()   # switches off the write to tiff file command if it has been used
}
}