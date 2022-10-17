library("haven")
set.seed(1)

data <- read_sas("cy07_msu_stu_cog.sas7bdat") 
#The dataset can be downloaded from https://webfs.oecd.org/pisa2018/SAS_STU_COG.zip

index <- which(data$OECD==1)
data <- data[index,]
variables <- read.csv("Items.csv") 
  
data = data[,c("CNT", variables$ItemID)]


index <- sample(1:nrow(data), 10000)
data <- data[index, ]

hist(rowSums(!is.na(data[,-1])))

count = rowSums(!is.na(data[,-1]))
data = data[count >= 20,]

ave = colMeans(data[,-1], na.rm=TRUE)
times = colSums(!is.na(data[,-1]))

data.cnt = data
data = data[,-1]
data = data[,ave>=0.1 & times >= 200]

variables=variables[ave>=0.1 & times >= 200,]
data.cnt = data.cnt[,c(1,which(ave>=0.1 & times >= 200))]


save(data.cnt, data, variables, file = "OECD2018.Rdata")
