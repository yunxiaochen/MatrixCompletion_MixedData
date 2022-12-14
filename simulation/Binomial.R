#==========================================Other useful functions

rownorm <- function(X){
  sqrt(rowSums(X^2))
}

#==========================================Functions NBE for binary/ordinal data

proj.max <- function(M, rho){
  M[M>rho] = rho
  M[M< -rho] = -rho 
  M
}


proj.nuc <- function(M, const){
  
  sigtot = sum(svd(M)$d)
  
  if(sigtot <= const){
    M
  }else{
    PiS(M/const) * const
  }
}

lik <- function(M, data, Omega, tot){
  sum(Omega * data * M) - tot*sum(Omega * log((1+exp(M))))
}


linesearch <- function(obj0, M0, grad, data, Omega,tot, rho, r, n, p, step = 1, times =15){
  
  n = nrow(data)
  
  M1 = M0 + step * grad
  M1 = proj.nuc(M1, rho * sqrt(r*n*p))
  M1 = proj.max(M1, rho)
  obj1 = lik(M1, data, Omega,tot)
  
  z = 0
  
  while(obj1 < obj0 & z < times){
    step = step *0.5
    z = z+1
    M1 = M0 + step * grad
    M1 = proj.nuc(M1, rho * sqrt(r*n*p))
    M1 = proj.max(M1, rho)
    obj1 = lik(M1, data, Omega,tot)
  }
  if(z == times){
    M1=M0
  }
  #print(z)
  M1
}

grad <- function(M0, data, Omega,tot){
  prob0 = 1/(1+exp(-M0))
  grad = (data - tot*prob0)/(tot*prob0*(1-prob0))*Omega
  grad
}

NBE <- function(rho, r, data, M0, Omega, tot, step = 1){
  
  obj0 = lik(M0, data, Omega, tot)
  n = nrow(data)
  p = ncol(data)
  
  grad = grad(M0, data, Omega,tot)
  
  M1 = linesearch(obj0, M0, grad, data, Omega,tot, rho, r, n, p)
  obj1 = lik(M1, data, Omega, tot)
  
  while(obj1-obj0 > 5e-3){
    #print(obj1-obj0)
    obj0 = obj1
    M0 = M1
    grad = grad(M0, data, Omega,tot)
    M1 = linesearch(obj0, M0, grad, data, Omega,tot, rho, r, n, p)
    obj1 = lik(M1, data, Omega, tot)
  }
  list(M = M1, obj = obj1)
}

#============================================Functions for refinement for binary/ordinal data; Without splitting

grad.Theta <- function(data, Omega, tot, A0, Theta0){
  M0 = Theta0%*% t(A0)
  prob0 = 1/(1+exp(-M0))
  grad.M = (data - tot*prob0)*Omega
  grad.Theta = grad.M %*% A0
  grad.Theta
}

grad2.Theta <- function(data, Omega, tot, A0, Theta0){
  M0 = Theta0%*% t(A0)
  prob0 = 1/(1+exp(-M0))
  grad2.M = tot*prob0 *(1-prob0)*Omega
  grad2.Theta = grad2.M %*% (A0^2)
  grad2.Theta
}

grad.A <- function(data, Omega, tot, A0, Theta0){
  M0 = Theta0%*% t(A0)
  prob0 = 1/(1+exp(-M0))
  grad.M = (data - tot*prob0)*Omega
  grad.A = t(grad.M) %*% Theta0
  grad.A
}

grad2.A <- function(data, Omega, tot, A0, Theta0){
  M0 = Theta0%*% t(A0)
  prob0 = 1/(1+exp(-M0))
  grad2.M = tot*prob0 *(1-prob0)*Omega
  grad2.A = t(grad2.M) %*% (Theta0^2)
  grad2.A
}

proj <- function(Theta1, C){
  
  Theta1.norm = sqrt(rowSums(Theta1^2))  
  Theta1[Theta1.norm>C,] = Theta1[Theta1.norm>C,]/Theta1.norm[Theta1.norm>C]*C
  Theta1
}

lik.row <- function(M, data, Omega, tot){
  rowSums(Omega * data * M) - tot*rowSums(Omega * log((1+exp(M))))
}

lik.col <- function(M, data, Omega, tot){
  colSums(Omega * data * M) - tot*colSums(Omega * log((1+exp(M))))
}

search.Theta <- function(grad, grad2, data, Omega, tot, A0, Theta0, step, times){
  n = nrow(Theta0)
  lik.row0 = lik.row(Theta0%*% t(A0), data, Omega, tot)
  step.vec = step * rep(1,n) 
  Theta1 = Theta0 + step.vec * grad/grad2
  lik.row1 = lik.row(Theta1%*% t(A0), data, Omega, tot)
  
  z = 0
  
  while( min(lik.row1 - lik.row0) < 0 & z< times){
    step.vec[lik.row1 - lik.row0 < 0] =step.vec[lik.row1 - lik.row0 < 0] /2
    z = z+1
    Theta1 = Theta0 + step.vec * grad/grad2
    lik.row1 = lik.row(Theta1%*% t(A0), data, Omega, tot)
  }
  Theta1[lik.row1 - lik.row0 < 0,] = Theta0[lik.row1 - lik.row0 < 0,]
  Theta1
}

search.A <- function(grad, grad2, data, Omega, tot, A0, Theta0, step, times){
  p = nrow(A0)
  lik.col0 = lik.col(Theta0%*% t(A0), data, Omega, tot)
  step.vec = step * rep(1,p) 
  A1 = A0 + step.vec * grad/grad2
  lik.col1 = lik.col(Theta0%*% t(A1), data, Omega, tot)
  
  z = 0
  
  while( min(lik.col1 - lik.col0) < 0 & z< times){
    step.vec[lik.col1 - lik.col0 < 0] =step.vec[lik.col1 - lik.col0 < 0] /2
    z = z+1
    A1 = A0 + step.vec * grad/grad2
    lik.col1 = lik.col(Theta0%*% t(A1), data, Omega, tot)
  }
  A1[lik.col1 - lik.col0 < 0,] = A0[lik.col1 - lik.col0 < 0,]
  A1
}

search.Theta.C <- function(grad, grad2, data, Omega, tot, A0, Theta0, C, step, times){
  n = nrow(Theta0)
  lik.row0 = lik.row(Theta0%*% t(A0), data, Omega, tot)
  step.vec = step * rep(1,n) 
  Theta1 = proj(Theta0 + step.vec * grad/grad2,C)
  lik.row1 = lik.row(Theta1%*% t(A0), data, Omega, tot)
  
  z = 0
  
  while( min(lik.row1 - lik.row0) < 0 & z< times){
    step.vec[lik.row1 - lik.row0 < 0] =step.vec[lik.row1 - lik.row0 < 0] /2
    z = z+1
    Theta1 = proj(Theta0 + step.vec * grad/grad2,C)
    lik.row1 = lik.row(Theta1%*% t(A0), data, Omega, tot)
  }
  Theta1[lik.row1 - lik.row0 < 0,] = Theta0[lik.row1 - lik.row0 < 0,]
  Theta1
}

search.A.C <- function(grad, grad2, data, Omega, tot, A0, Theta0, C, step, times){
  p = nrow(A0)
  lik.col0 = lik.col(Theta0%*% t(A0), data, Omega, tot)
  step.vec = step * rep(1,p) 
  A1 = proj(A0 + step.vec * grad/grad2,C)
  lik.col1 = lik.col(Theta0%*% t(A1), data, Omega, tot)
  
  z = 0
  
  while( min(lik.col1 - lik.col0) < 0 & z< times){
    step.vec[lik.col1 - lik.col0 < 0] =step.vec[lik.col1 - lik.col0 < 0] /2
    z = z+1
    A1 = proj(A0 + step.vec * grad/grad2,C)
    lik.col1 = lik.col(Theta0%*% t(A1), data, Omega, tot)
  }
  A1[lik.col1 - lik.col0 < 0,] = A0[lik.col1 - lik.col0 < 0,]
  #print(lik(Theta0%*% t(A1), data, Omega, tot))
  
  A1
}


update.Theta <- function(data, Omega, tot, A0, Theta0, step = 1, times = 10){
  
  obj0 = lik(Theta0%*% t(A0), data, Omega, tot)
  grad = grad.Theta(data, Omega, tot, A0, Theta0)
  grad2 = grad2.Theta(data, Omega, tot, A0, Theta0)
  
  Theta1 =  search.Theta(grad, grad2, data, Omega, tot, A0, Theta0, step, times)
  obj1 = lik(Theta1%*% t(A0), data, Omega, tot)
  
  while(abs(obj1-obj0) > 1e-3){
    obj0 = obj1
    Theta0 = Theta1
    grad = grad.Theta(data, Omega, tot, A0, Theta0)
    grad2 = grad2.Theta(data, Omega, tot, A0, Theta0)
    Theta1 =  search.Theta(grad, grad2, data, Omega, tot, A0, Theta0, step, times)
    obj1 = lik(Theta1%*% t(A0), data, Omega, tot)
  }
  Theta1
}


update.A <- function(data, Omega, tot, A0, Theta0, step = 1, times = 10){
  obj0 = lik(Theta0%*% t(A0), data, Omega, tot)
  grad = grad.A(data, Omega, tot, A0, Theta0)
  grad2 = grad2.A(data, Omega, tot, A0, Theta0)
  A1 = search.A(grad, grad2, data, Omega, tot, A0, Theta0, step, times)
  obj1 = lik(Theta0%*% t(A1), data, Omega, tot)
  while(abs(obj1-obj0) > 1e-3){
    obj0 = obj1
    A0 = A1
    grad = grad.A(data, Omega, tot, A0, Theta0)
    grad2 = grad2.A(data, Omega, tot, A0, Theta0)
    A1 = search.A(grad, grad2, data, Omega, tot, A0, Theta0, step, times)
    obj1 = lik(Theta0%*% t(A1), data, Omega, tot)
  }
  A1
}


update.Theta.C <- function(data, Omega, tot, A0, Theta0, C, step = 1, times = 10){
  
  obj0 = lik(Theta0%*% t(A0), data, Omega, tot)
  
  grad = grad.Theta(data, Omega, tot, A0, Theta0)
  grad2 = grad2.Theta(data, Omega, tot, A0, Theta0)
  
  Theta1 =  search.Theta.C(grad, grad2, data, Omega, tot, A0, Theta0, C, step, times)
  obj1 = lik(Theta1%*% t(A0), data, Omega, tot)
  
  while(abs(obj1-obj0) > 1e-3){
    obj0 = obj1
    Theta0 = Theta1
    grad = grad.Theta(data, Omega, tot, A0, Theta0)
    grad2 = grad2.Theta(data, Omega, tot, A0, Theta0)
    Theta1 =  search.Theta.C(grad, grad2, data, Omega, tot, A0, Theta0,C, step, times)
    obj1 = lik(Theta1%*% t(A0), data, Omega, tot)
  }
  Theta1
}


update.A.C <- function(data, Omega, tot, A0, Theta0, C, step = 1, times = 10){
  obj0 = lik(Theta0%*% t(A0), data, Omega, tot)
  
  grad = grad.A(data, Omega, tot, A0, Theta0)
  grad2 = grad2.A(data, Omega, tot, A0, Theta0)
  A1 = search.A.C(grad, grad2, data, Omega, tot, A0, Theta0, C, step, times)
  obj1 = lik(Theta0%*% t(A1), data, Omega, tot)
  
  while(abs(obj1-obj0) > 1e-3){
    obj0 = obj1
    A0 = A1
    grad = grad.A(data, Omega, tot, A0, Theta0)
    grad2 = grad2.A(data, Omega, tot, A0, Theta0)
    A1 = search.A.C(grad, grad2, data, Omega, tot, A0, Theta0, C, step, times)
    obj1 = lik(Theta0%*% t(A1), data, Omega, tot)
  }
  A1
}

refi.nosp <- function(M, r, data, Omega,  C2, tot){
  
  svd.res = svd(round(M,8))
  A0 = proj(svd.res$v[,1:r], C2)
  Theta0 = svd.res$u[,1:r] %*% diag(svd.res$d[1:r] )
  Theta0 = update.Theta(data, Omega, tot, A0, Theta0)
  A0 = update.A(data, Omega, tot, A0, Theta0)
  list(M = Theta0%*% t(A0))
}

#===================================================Functions for CJMLE


CJMLE <- function(Theta0, A0, r, data, Omega,  C,tot){
  
  obj0 = lik(Theta0%*% t(A0), data, Omega, tot)
  Theta0 = update.Theta.C(data, Omega, tot, A0, Theta0,C)
  A0 = update.A.C(data, Omega, tot, A0, Theta0,C)
  
  obj1 = lik(Theta0%*% t(A0), data, Omega, tot)
  
  while(abs(obj1-obj0) > 1e-2){
    obj0 = obj1
    Theta0 = update.Theta.C(data, Omega, tot, A0, Theta0,C)
    A0 = update.A.C(data, Omega, tot, A0, Theta0,C)
    obj1 = lik(Theta0%*% t(A0), data, Omega, tot)
  #  print(obj1)
  }
  list(M = Theta0 %*% t(A0), Theta = Theta0, A = A0)  
}



#===================================================Functions for refinement for binary/ordinal data; With splitting


refi.sp.nbe <- function(M, r, data, Omega,  C2,tot,rho){
  
  n=nrow(data)
  split = rbinom(n, 1, 0.5)
  
  data1 = data[split==0,]
  data2 = data[split==1,]
  
  Omega1 = Omega[split==0,]
  Omega2 = Omega[split==1,]
  
  M1 = NBE(rho, r, data1, M[split==0,], Omega1, tot)$M
  M2 = NBE(rho, r, data2, M[split==1,], Omega2, tot)$M
  
  svd.res1 = svd(round(M1,8))
  A10 = proj(svd.res1$v[,1:r], C2)
  Theta10 = svd.res1$u[,1:r] %*% diag(svd.res1$d[1:r])
  
  
  svd.res2 = svd(round(M2,8))
  A20 = proj(svd.res2$v[,1:r], C2)
  Theta20 = svd.res2$u[,1:r] %*% diag(svd.res2$d[1:r])
  
  
  Theta2 = update.Theta(data2, Omega2, tot, A10, Theta20)
  A1 = update.A(data2, Omega2, tot, A10, Theta2)
  
  
  Theta1 = update.Theta(data1, Omega1, tot, A20, Theta10)
  A2 = update.A(data1, Omega1, tot, A20, Theta1)
  
  M[split==0,] = Theta1 %*% t(A2)
  M[split==1,] = Theta2 %*% t(A1)
  
  list(M = M)
}



refi.sp.jml <- function(Theta0, A0, r, data, Omega,  C2,tot,C){
  
  
  
  
  n=nrow(data)
  split = rbinom(n, 1, 0.5)
  
  data1 = data[split==0,]
  data2 = data[split==1,]
  
  Omega1 = Omega[split==0,]
  Omega2 = Omega[split==1,]
  
  M = Theta0 %*% t(A0)
  est1 = CJMLE(Theta0[split==0,], A0, r, data1, Omega1,  C,tot)
  est2 = CJMLE(Theta0[split==1,], A0, r, data2, Omega2,  C,tot)
  
#  svd.res1 = svd(round(M1,8))
#  A10 = proj(svd.res1$v[,1:r], C2)
#  Theta10 = svd.res1$u[,1:r] %*% diag(svd.res1$d[1:r])
  
  
#  svd.res2 = svd(round(M2,8))
#  A20 = proj(svd.res2$v[,1:r], C2)
#  Theta20 = svd.res2$u[,1:r] %*% diag(svd.res2$d[1:r])
  
  
  Theta2 = update.Theta(data2, Omega2, tot, est1$A, est2$Theta)
  A1 = update.A(data2, Omega2, tot, est1$A, Theta2)
  
  
  Theta1 = update.Theta(data1, Omega1, tot, est2$A, est1$Theta)
  A2 = update.A(data1, Omega1, tot, est2$A, Theta1)
  
  M[split==0,] = Theta1 %*% t(A2)
  M[split==1,] = Theta2 %*% t(A1)
  
  list(M = M)
}

#====================Evaluation

F.accu <- function(M.hat, M){
  sqrt(mean((M.hat-M)^2))
}

M.accu  <- function(M.hat, M){
  max(abs(M.hat-M)) 
}
