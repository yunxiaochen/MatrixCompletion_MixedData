simulation <- function(i){
  
  #for(i in 1:100){
  set.seed(i)
  
  source("Mixed.R")
  max.norm <- matrix(0,3,8)
  f.norm <- matrix(0,3,8)
  
  
  
  #setting1
  n=400
  p=200
  pi = 0.2
  r = 5
  tot = 5
  rho = r
  
  Theta <- matrix(runif(n*r,-0.9,0.9),n,r)
  A <- matrix(runif(p*r,-0.9,0.9),p,r)
  
  C2 <- 2*sqrt(r/p)
  C = sqrt(r)
  
  temp = Theta %*% t(A)
  prob = 1/(1+exp(-temp))
  
  data <- cbind(matrix(rbinom(n*p,tot,prob), n,p/2), matrix(temp[,-(1:(p/2))] + rnorm(n*(p/2)), n,p/2))
  Omega <- matrix(rbinom(n*p,1,pi), n,p) 
  
  M0 = temp
  
  
  res1 = NBE(rho, r, data, M0, Omega, tot)
  res2 = refi.nosp(res1$M, r, data, Omega,  C2,tot)
  
  res3.1 = refi.sp.nbe(res1$M, r, data, Omega,  C2,tot,rho)
  res3.2 = refi.sp.nbe(res1$M, r, data, Omega,  C2,tot,rho)
  res3.3 = refi.sp.nbe(res1$M, r, data, Omega,  C2,tot,rho)
  res3.4 = refi.sp.nbe(res1$M, r, data, Omega,  C2,tot,rho)
  res3.5 = refi.sp.nbe(res1$M, r, data, Omega,  C2,tot,rho)
  
  res4 = (res3.1$M + res3.2$M + res3.3$M + res3.4$M+ res3.5$M)/5
  
  svd.res = svd(round(res2$M,8))
  Theta0 = proj((svd.res$u[,1:r] *sqrt(n)) %*% diag(sqrt(svd.res$d[1:r]/sqrt(n*p))),C)
  A0 = proj((svd.res$v[,1:r] *sqrt(p)) %*% diag(sqrt(svd.res$d[1:r]/sqrt(n*p))),C)
  
  res5 = CJMLE(Theta0, A0, r, data, Omega,  C,tot)
  
  res6 = refi.nosp(res5$M, r, data, Omega,  C2,tot)
  
  res7.1 = refi.sp.jml(res5$Theta,res5$A, r, data, Omega,  C2,tot,C)
  res7.2 = refi.sp.jml(res5$Theta,res5$A, r, data, Omega,  C2,tot,C)
  res7.3 = refi.sp.jml(res5$Theta,res5$A, r, data, Omega,  C2,tot,C)
  res7.4 = refi.sp.jml(res5$Theta,res5$A, r, data, Omega,  C2,tot,C)
  res7.5 = refi.sp.jml(res5$Theta,res5$A, r, data, Omega,  C2,tot,C)
  
  
  res8 = (res7.1$M + res7.2$M+ res7.3$M+ res7.4$M+ res7.5$M)/5
  
  
  max.norm[1,1] = M.accu(res1$M, temp)
  max.norm[1,2] = M.accu(res2$M, temp)
  max.norm[1,3] = M.accu(res3.1$M, temp)
  max.norm[1,4] = M.accu(res4, temp)
  max.norm[1,5] = M.accu(res5$M, temp)
  max.norm[1,6] = M.accu(res6$M, temp)
  max.norm[1,7] = M.accu(res7.1$M, temp)
  max.norm[1,8] = M.accu(res8, temp) 
  
  
  f.norm[1,1] = F.accu(res1$M, temp)
  f.norm[1,2] = F.accu(res2$M, temp)
  f.norm[1,3] = F.accu(res3.1$M, temp)
  f.norm[1,4] = F.accu(res4, temp)
  f.norm[1,5] = F.accu(res5$M, temp)
  f.norm[1,6] = F.accu(res6$M, temp)
  f.norm[1,7] = F.accu(res7.1$M, temp)
  f.norm[1,8] = F.accu(res8, temp) 
  
  #setting2
  n=800
  p=400 
  tot = 5
  rho = r
  
  Theta <- matrix(runif(n*r,-0.9,0.9),n,r)
  A <- matrix(runif(p*r,-0.9,0.9),p,r)
  
  C2 <- 2*sqrt(r/p)
  C = sqrt(r)
  
  temp = Theta %*% t(A)
  prob = 1/(1+exp(-temp))
  
  data <- cbind(matrix(rbinom(n*p,tot,prob), n,p/2), matrix(temp[,-(1:(p/2))] + rnorm(n*(p/2)), n,p/2))
  Omega <- matrix(rbinom(n*p,1,pi), n,p) 
  
  
  M0 = temp
  
  
  res1 = NBE(rho, r, data, M0, Omega, tot)
  res2 = refi.nosp(res1$M, r, data, Omega,  C2,tot)
  
  res3.1 = refi.sp.nbe(res1$M, r, data, Omega,  C2,tot,rho)
  res3.2 = refi.sp.nbe(res1$M, r, data, Omega,  C2,tot,rho)
  res3.3 = refi.sp.nbe(res1$M, r, data, Omega,  C2,tot,rho)
  res3.4 = refi.sp.nbe(res1$M, r, data, Omega,  C2,tot,rho)
  res3.5 = refi.sp.nbe(res1$M, r, data, Omega,  C2,tot,rho)
  
  res4 = (res3.1$M + res3.2$M + res3.3$M + res3.4$M+ res3.5$M)/5
  
  svd.res = svd(round(res2$M,8))
  Theta0 = proj((svd.res$u[,1:r] *sqrt(n)) %*% diag(sqrt(svd.res$d[1:r]/sqrt(n*p))),C)
  A0 = proj((svd.res$v[,1:r] *sqrt(p)) %*% diag(sqrt(svd.res$d[1:r]/sqrt(n*p))),C)
  
  res5 = CJMLE(Theta0, A0, r, data, Omega,  C,tot)
  
  res6 = refi.nosp(res5$M, r, data, Omega,  C2,tot)
  
  res7.1 = refi.sp.jml(res5$Theta,res5$A, r, data, Omega,  C2,tot,C)
  res7.2 = refi.sp.jml(res5$Theta,res5$A, r, data, Omega,  C2,tot,C)
  res7.3 = refi.sp.jml(res5$Theta,res5$A, r, data, Omega,  C2,tot,C)
  res7.4 = refi.sp.jml(res5$Theta,res5$A, r, data, Omega,  C2,tot,C)
  res7.5 = refi.sp.jml(res5$Theta,res5$A, r, data, Omega,  C2,tot,C)
  
  
  res8 = (res7.1$M + res7.2$M+ res7.3$M+ res7.4$M+ res7.5$M)/5
  
  
  max.norm[2,1] = M.accu(res1$M, temp)
  max.norm[2,2] = M.accu(res2$M, temp)
  max.norm[2,3] = M.accu(res3.1$M, temp)
  max.norm[2,4] = M.accu(res4, temp)
  max.norm[2,5] = M.accu(res5$M, temp)
  max.norm[2,6] = M.accu(res6$M, temp)
  max.norm[2,7] = M.accu(res7.1$M, temp)
  max.norm[2,8] = M.accu(res8, temp) 
  
  
  f.norm[2,1] = F.accu(res1$M, temp)
  f.norm[2,2] = F.accu(res2$M, temp)
  f.norm[2,3] = F.accu(res3.1$M, temp)
  f.norm[2,4] = F.accu(res4, temp)
  f.norm[2,5] = F.accu(res5$M, temp)
  f.norm[2,6] = F.accu(res6$M, temp)
  f.norm[2,7] = F.accu(res7.1$M, temp)
  f.norm[2,8] = F.accu(res8, temp)   
  
  #setting3
  n=1600
  p=800 
  tot = 5
  rho = r
  
  Theta <- matrix(runif(n*r,-0.9,0.9),n,r)
  A <- matrix(runif(p*r,-0.9,0.9),p,r)
  
  C2 <- 2*sqrt(r/p)
  C = sqrt(r)
  
  temp = Theta %*% t(A)
  prob = 1/(1+exp(-temp))
  
  data <- cbind(matrix(rbinom(n*p,tot,prob), n,p/2), matrix(temp[,-(1:(p/2))] + rnorm(n*(p/2)), n,p/2))
  Omega <- matrix(rbinom(n*p,1,pi), n,p) 
  
  
  M0 = temp
  
  
  res1 = NBE(rho, r, data, M0, Omega, tot)
  res2 = refi.nosp(res1$M, r, data, Omega,  C2,tot)
  
  res3.1 = refi.sp.nbe(res1$M, r, data, Omega,  C2,tot,rho)
  res3.2 = refi.sp.nbe(res1$M, r, data, Omega,  C2,tot,rho)
  res3.3 = refi.sp.nbe(res1$M, r, data, Omega,  C2,tot,rho)
  res3.4 = refi.sp.nbe(res1$M, r, data, Omega,  C2,tot,rho)
  res3.5 = refi.sp.nbe(res1$M, r, data, Omega,  C2,tot,rho)
  
  res4 = (res3.1$M + res3.2$M + res3.3$M + res3.4$M+ res3.5$M)/5
  
  svd.res = svd(round(res2$M,8))
  Theta0 = proj((svd.res$u[,1:r] *sqrt(n)) %*% diag(sqrt(svd.res$d[1:r]/sqrt(n*p))),C)
  A0 = proj((svd.res$v[,1:r] *sqrt(p)) %*% diag(sqrt(svd.res$d[1:r]/sqrt(n*p))),C)
  
  res5 = CJMLE(Theta0, A0, r, data, Omega,  C,tot)
  
  res6 = refi.nosp(res5$M, r, data, Omega,  C2,tot)
  
  res7.1 = refi.sp.jml(res5$Theta,res5$A, r, data, Omega,  C2,tot,C)
  res7.2 = refi.sp.jml(res5$Theta,res5$A, r, data, Omega,  C2,tot,C)
  res7.3 = refi.sp.jml(res5$Theta,res5$A, r, data, Omega,  C2,tot,C)
  res7.4 = refi.sp.jml(res5$Theta,res5$A, r, data, Omega,  C2,tot,C)
  res7.5 = refi.sp.jml(res5$Theta,res5$A, r, data, Omega,  C2,tot,C)
  
  
  res8 = (res7.1$M + res7.2$M+ res7.3$M+ res7.4$M+ res7.5$M)/5
  
  
  max.norm[3,1] = M.accu(res1$M, temp)
  max.norm[3,2] = M.accu(res2$M, temp)
  max.norm[3,3] = M.accu(res3.1$M, temp)
  max.norm[3,4] = M.accu(res4, temp)
  max.norm[3,5] = M.accu(res5$M, temp)
  max.norm[3,6] = M.accu(res6$M, temp)
  max.norm[3,7] = M.accu(res7.1$M, temp)
  max.norm[3,8] = M.accu(res8, temp) 
  
  
  f.norm[3,1] = F.accu(res1$M, temp)
  f.norm[3,2] = F.accu(res2$M, temp)
  f.norm[3,3] = F.accu(res3.1$M, temp)
  f.norm[3,4] = F.accu(res4, temp)
  f.norm[3,5] = F.accu(res5$M, temp)
  f.norm[3,6] = F.accu(res6$M, temp)
  f.norm[3,7] = F.accu(res7.1$M, temp)
  f.norm[3,8] = F.accu(res8, temp) 
  
  filename = paste("sim", i, ".Rdata", sep="")
  save(max.norm, f.norm, file = filename)
  
  #}
}

library(InspectChangepoint)
library(parallel)

r = mclapply(1:100, simulation, mc.cores=9)
