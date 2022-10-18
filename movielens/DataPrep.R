raw <- read.table(file = "udata.txt", head=FALSE)

n = 943
p = 1682

Omega <- matrix(0,n,p)
data <- matrix(0,n,p)


for(i in 1:nrow(raw)){
  Omega[raw[i,1],raw[i,2]] = 1
  data[raw[i,1],raw[i,2]] = raw[i,3]-1
}

save(data, Omega, file = "movielens100K.Rdata")
