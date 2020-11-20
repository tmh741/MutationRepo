library(rstan)
library(bayesplot)
library(ggplot2)
library(rstanarm)
library(tidyverse)
library(reshape2)
library(gridExtra)

data = readRDS("lung_sbs_matrix.rds")
store = as.matrix(t(data[,1]))
range(store)

model = "
data{
int<lower=0> N;
int<lower=0> S;
int<lower=0> K;
int<lower=0> y[N,S];
real<lower=0> a;
real<lower=0> b;
}

parameters{
  positive_ordered[K] theta[N];
  vector<lower=0>[K] beta[S];
}

model {
  for (n in 1:N)
    theta[n] ~ gamma(a,b);
  for (s in 1:(S))
    beta[s] ~ gamma(a,b);
    
  for (n in 1:N) {
    for(s in 1:S) {
      y[n,s] ~ poisson(theta[n]'*beta[s]);
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
a0= 1/K

please = rstan::stan(model_code=model,
                     data = list(N=N, S=S,K=K,y=y,a=a,b=b),iter=6000)
stan_rhat(please)

vals = as.matrix(please)
dim(vals)
displayvalues=data.frame(names=rep(NA,ncol(vals)),value=rep(NA,ncol(vals)),madsd=rep(NA,ncol(vals)),dim1=rep(NA,ncol(vals)),dim2=rep(NA,ncol(vals)),variable=rep(NA,ncol(vals)))
for(i in 1:ncol(vals)){
  displayvalues$names[i] = colnames(vals)[i]
  displayvalues$value[i] = median(vals[,i])
  displayvalues$madsd[i] = mad(vals[,i])
  displayvalues$dim1[i] = gsub("[^0-9.-]","",strsplit(colnames(vals)[i],",")[[1]])[1]
  displayvalues$dim2[i] = gsub("[^0-9.-]","",strsplit(colnames(vals)[i],",")[[1]])[2]
  displayvalues$variable[i] = gsub("[^a-zA-Z]","",colnames(vals)[i])
}
displayvalues$dim1 = as.numeric(displayvalues$dim1)
displayvalues$dim2 = as.numeric(displayvalues$dim2)


betamatrix = dcast(displayvalues %>% filter(variable=="beta"), dim1 ~ dim2)[,-1]
colnames(betamatrix) = paste("sig",1:ncol(betamatrix),sep="")
betamatrix$names = colnames(store)
betamatrix$group = betamatrix$names %>% strsplit("_") %>% sapply("[",1)

p1 = ggplot(betamatrix) + aes(x=names,y=sig1,fill=group) + geom_bar(stat="identity") + theme(axis.text.x = element_blank(),legend.position="none") + xlab("") + ylab("Signature 1")
p2 = ggplot(betamatrix) + aes(x=names,y=sig2,fill=group) + geom_bar(stat="identity") + theme(axis.text.x = element_blank(),legend.position="none") + xlab("") + ylab("Signature 2")
p3 = ggplot(betamatrix) + aes(x=names,y=sig3,fill=group) + geom_bar(stat="identity") + theme(axis.text.x = element_blank(),legend.position="none") + xlab("") + ylab("Signature 3")
p4 = ggplot(betamatrix) + aes(x=names,y=sig4,fill=group) + geom_bar(stat="identity") + theme(axis.text.x = element_blank(),legend.position="none") + xlab("") + ylab("Signature 4")
p5 = ggplot(betamatrix) + aes(x=names,y=sig5,fill=group) + geom_bar(stat="identity") + theme(axis.text.x = element_blank(),legend.position="none") + xlab("") + ylab("Signature 5")

grid.arrange(p1,p2,p3,p4,p5,ncol=1)
