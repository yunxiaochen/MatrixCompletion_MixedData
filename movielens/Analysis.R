set.seed(1)
load("movielens100K.Rdata")
source("Binomial.R")
#Split data into 80% training and 20% testing. 
library(InspectChangepoint)
library(parallel)
data = as.matrix(data)
obs.entry = which(Omega==1)
test.num = round(length(obs.entry)/5)

test.ind = sample(1:length(obs.entry), test.num)
test.entry = obs.entry[test.ind]
train.entry = obs.entry[-test.ind]

data.train = data
data.test = data
 
 
n=nrow(data)
p=ncol(data)
tot = rep(1,n) %*% t(rep(4,p))

Omega.train <- Omega 
Omega.test <- Omega

Omega.train[test.entry] = 0
Omega.test[train.entry] = 0

M0 = matrix(runif(n*p, -0.5, 0.5),n,p)

r.vec = seq(10, 90, by = 10)

 

analysis <- function(r){
  
  set.seed(1)
  print(r)
  
  result = rep(0, 8)
  rho = 4*r
  C2 <- 2*sqrt(r/p)
  C = 2*sqrt(r)
  
  res1 = NBE(rho, r, data.train, M0, Omega, tot)
  result[1] = lik(res1$M,data.test,Omega.test,tot)
  
  res2 = refi.nosp(res1$M, r, data.train, Omega,  C2,tot,C)
  result[2] = lik(res2$M,data.test,Omega.test,tot)
  
  res3.1 = refi.sp(res1$M, r, data.train, Omega,  C2,tot,C)
  res3.2 = refi.sp(res1$M, r, data.train, Omega,  C2,tot,C)
  res3.3 = refi.sp(res1$M, r, data.train, Omega,  C2,tot,C)
  res3.4 = refi.sp(res1$M, r, data.train, Omega,  C2,tot,C)
  res3.5 = refi.sp(res1$M, r, data.train, Omega,  C2,tot,C)
  
  result[3] = lik(res3.1$M,data.test,Omega.test,tot)
  
  res4 = (res3.1$M + res3.2$M + res3.3$M + res3.4$M+ res3.5$M)/5
  
  result[4] = lik(res4,data.test,Omega.test,tot)
  
  svd.res = svd(round(res4,8))
  Theta0 = proj((matrix(svd.res$u[,1:r], ncol=r) *sqrt(n)) %*% diag(sqrt(svd.res$d[1:r]/sqrt(n*p)),ncol=r,nrow=r),C)
  A0 = proj((matrix(svd.res$v[,1:r],ncol=r) *sqrt(p)) %*% diag(sqrt(svd.res$d[1:r]/sqrt(n*p)),nrow = r,ncol=r),C)

  res5 = CJMLE(Theta0, A0, r, data.train, Omega,  C,tot)
  result[5] = lik(res5$M,data.test,Omega.test,tot)
  
  res6 = refi.nosp(res5$M, r, data.train, Omega,  C2,tot,C)
  result[6] = lik(res6$M,data.test,Omega.test,tot)
  
  res7.1 = refi.sp(res5$M, r, data.train, Omega,  C2,tot,C)
  res7.2 = refi.sp(res5$M, r, data.train, Omega,  C2,tot,C)
  res7.3 = refi.sp(res5$M, r, data.train, Omega,  C2,tot,C)
  res7.4 = refi.sp(res5$M, r, data.train, Omega,  C2,tot,C)
  res7.5 = refi.sp(res5$M, r, data.train, Omega,  C2,tot,C)
  result[7] = lik(res7.1$M,data.test,Omega.test,tot)
  
  
  res8 = (res7.1$M + res7.2$M+ res7.3$M+ res7.4$M+ res7.5$M)/5
  
  result[8] = lik(res8,data.test,Omega.test,tot)
  filename = paste("rank", r, ".Rdata", sep="")
  save(result, file = filename)
}

r = mclapply(r.vec, analysis, mc.cores=9)

 
r.vec = seq(5, 90, by = 10)
C.vec = c(1,sqrt(2),2,3,4,5)
for(r in r.vec){
  for(C in C.vec){
    Theta0 = proj((matrix(svd.res$u[,1:r], ncol=r) *sqrt(n)) %*% diag(sqrt(svd.res$d[1:r]/sqrt(n*p)),ncol=r,nrow=r),C)
    A0 = proj((matrix(svd.res$v[,1:r],ncol=r) *sqrt(p)) %*% diag(sqrt(svd.res$d[1:r]/sqrt(n*p)),nrow = r,ncol=r),C)
    
    res5 = CJMLE(Theta0, A0, r, data.train, Omega,  C,tot)
    result[5] = lik(res5$M,data.test,Omega.test,tot)
    print(r)
    print(C)
    print(lik(res5$M,data.test,Omega.test,tot))
  }
} 


  