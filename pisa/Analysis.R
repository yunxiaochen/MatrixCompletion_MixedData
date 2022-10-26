set.seed(1)
load("OECD2018.Rdata")
source("Binomial.R")
#Split data into 80% training and 20% testing. 
library(InspectChangepoint)

data = as.matrix(data)
obs.entry = which(!is.na(data))
test.num = round(length(obs.entry)/5)

test.ind = sample(1:length(obs.entry), test.num)
test.entry = obs.entry[test.ind]
train.entry = obs.entry[-test.ind]

data.train = data
data.train[test.entry] = NA

data.test = data
data.test[train.entry] = NA

r.vec = 1:4

n=nrow(data)
p=ncol(data)
tot = rep(1,n) %*% t(variables$tot)

Omega <- !is.na(data.train)
data.train[is.na(data.train)] = 0

Omega.test <- !is.na(data.test)
data.test[is.na(data.test)] = 0

M0 = data.train
M0[data.train==0] = 0.2
M0[data.train!=0] = 0.8



analysis <- function(r){
  
  set.seed(1)
  print(r)
  
  result = rep(0, 8)
  C2 <- 2*sqrt(r/p)
  C = 2*sqrt(r)
  rho = C^2
  
  res1 = NBE(rho, r, data.train, M0, Omega, tot)
  result[1] = lik(res1$M,data.test,Omega.test,tot)
  
  res2 = refi.nosp(res1$M, r, data.train, Omega,  C2,tot,C)
  result[2] = lik(res2$M,data.test,Omega.test,tot)
  
  res3.1 = refi.sp.nbe(res1$M, r, data.train, Omega,  C2,tot,rho)
  res3.2 = refi.sp.nbe(res1$M, r, data.train, Omega,  C2,tot,rho)
  res3.3 = refi.sp.nbe(res1$M, r, data.train, Omega,  C2,tot,rho)
  res3.4 = refi.sp.nbe(res1$M, r, data.train, Omega,  C2,tot,rho)
  res3.5 = refi.sp.nbe(res1$M, r, data.train, Omega,  C2,tot,rho)
  
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
  
  res7.1 = refi.sp.jml(res5$Theta,res5$A, r, data.train, Omega,  C2,tot,C)
  res7.2 = refi.sp.jml(res5$Theta,res5$A, r, data.train, Omega,  C2,tot,C)
  res7.3 = refi.sp.jml(res5$Theta,res5$A, r, data.train, Omega,  C2,tot,C)
  res7.4 = refi.sp.jml(res5$Theta,res5$A, r, data.train, Omega,  C2,tot,C)
  res7.5 = refi.sp.jml(res5$Theta,res5$A, r, data.train, Omega,  C2,tot,C)
  
  result[7] = lik(res7.1$M,data.test,Omega.test,tot)

  res8 = (res7.1$M + res7.2$M+ res7.3$M+ res7.4$M+ res7.5$M)/5
  
  result[8] = lik(res8,data.test,Omega.test,tot)
  filename = paste("rank", r, ".Rdata", sep="")
  save(result, file = filename)
}
library(parallel)
r = mclapply(r.vec, analysis, mc.cores=4)
