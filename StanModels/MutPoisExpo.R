library(rstan)
library(bayesplot)
library(ggplot2)
library(rstanarm)
library(tidyverse)
library(reshape2)
library(gridExtra)

data = readRDS("lung_sbs_matrix.rds")
store = data[,1:10]
range(store)

model = "
data{
  int<lower=0> U;
  int<lower=0> I;
  int<lower=0> K;
  int<lower=0> y[U,I];
  real<lower=0> lambda0;
  real<lower=0> alpha0;
}

transformed data{
  vector<lower=0>[K] alpha0_vec;
  for (k in 1:K){
    alpha0_vec[k] = alpha0;
  }
}

parameters {
  simplex[K] theta[U];
  vector<lower=0>[K] beta[I];
}

model {
  for (u in 1:U)
    theta[u] ~ dirichlet(alpha0_vec);
  for (i in 1:I)
    beta[i] ~ exponential(lambda0);
    
  for(u in 1:U){
    for(i in 1:I){
      y[u,i] ~ poisson(theta[u]'*beta[i]);
    }
  }
}
"
storerows = rowSums(store)
storecols = colSums(store)

N = nrow(store)
S = ncol(store)
K = 5
y = store
a = 1000
b = 0.00001
alpha= 1/K

please = rstan::stan(model_code=model,
                     data = list(U=N, I=S,K=K,y=y,lambda0=0.1,alpha0=alpha),iter=6000)

vals2 = as.matrix(please)
dim(vals2)
displayval2=data.frame(names=rep(NA,ncol(vals2)),value=rep(NA,ncol(vals2)),madsd=rep(NA,ncol(vals2)),dim1=rep(NA,ncol(vals2)),dim2=rep(NA,ncol(vals2)),variable=rep(NA,ncol(vals2)))
for(i in 1:ncol(vals2)){
  displayval2$names[i] = colnames(vals2)[i]
  displayval2$value[i] = median(vals2[,i])
  displayval2$madsd[i] = mad(vals2[,i])
  displayvalues$dim1[i] = gsub("[^0-9.-]","",strsplit(colnames(vals2)[i],",")[[1]])[1]
  displayvalues$dim2[i] = gsub("[^0-9.-]","",strsplit(colnames(vals2)[i],",")[[1]])[2]
  displayvalues$variable[i] = gsub("[^a-zA-Z]","",colnames(vals2)[i])
}
displayvalues$dim1 = as.numeric(displayvalues$dim1)
displayvalues$dim2 = as.numeric(displayvalues$dim2)


thetamatrix = dcast(displayvalues %>% filter(variable=="theta"), dim1 ~ dim2)[,-1]
colnames(thetamatrix) = paste("sig",1:ncol(thetamatrix),sep="")
thetamatrix$names = rownames(store)
thetamatrix$group = thetamatrix$names %>% strsplit("_") %>% sapply("[",1)

p1 = ggplot(thetamatrix) + aes(x=names,y=sig1,fill=group) + geom_bar(stat="identity") + theme(axis.text.x = element_blank(),legend.position="none") + xlab("") + ylab("Signature 1")
p2 = ggplot(thetamatrix) + aes(x=names,y=sig2,fill=group) + geom_bar(stat="identity") + theme(axis.text.x = element_blank(),legend.position="none") + xlab("") + ylab("Signature 2")
p3 = ggplot(thetamatrix) + aes(x=names,y=sig3,fill=group) + geom_bar(stat="identity") + theme(axis.text.x = element_blank(),legend.position="none") + xlab("") + ylab("Signature 3")
p4 = ggplot(thetamatrix) + aes(x=names,y=sig4,fill=group) + geom_bar(stat="identity") + theme(axis.text.x = element_blank(),legend.position="none") + xlab("") + ylab("Signature 4")
p5 = ggplot(thetamatrix) + aes(x=names,y=sig5,fill=group) + geom_bar(stat="identity") + theme(axis.text.x = element_blank(),legend.position="none") + xlab("") + ylab("Signature 5")

grid.arrange(p1,p2,p3,p4,p5,ncol=1)


