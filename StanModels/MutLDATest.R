library(rstan)
library(bayesplot)
library(ggplot2)
library(rstanarm)
library(tidyverse)
library(reshape2)
library(gridExtra)

data = data.frame(t(readRDS("lung_sbs_matrix.rds")))[1:10,]
data$doc = rownames(data)
datalong = data %>% group_by(doc) %>% pivot_longer(cols=colnames(data)[1:96],names_to="Mutation") %>% filter(value!=0)
datalong$Mutation = as.numeric(as.factor(datalong$Mutation))
datalong$doc = as.numeric(as.factor(datalong$doc))


lda = "
data{
  int<lower=2> K;               // num topics
  int<lower=2> V;               // num words
  int<lower=1> M;               // num docs
  int<lower=1> N;               // total word instances
  int<lower=1,upper=V> w[N];    // word n
  int<lower=1,upper=M> doc[N];  // doc ID for word n
  vector<lower=0>[K] alpha;     // topic prior
  vector<lower=0>[V] beta;      // word prior
}
parameters {
  simplex[K] theta[M];   // topic dist for doc m
  simplex[V] phi[K];     // word dist for topic k
}

model {
  for (m in 1:M)
    theta[m] ~ dirichlet(alpha);  // prior
  for (k in 1:K)
    phi[k] ~ dirichlet(beta);     // prior
  for (n in 1:N) {
    real gamma[K];
    for (k in 1:K)
      gamma[k] = log(theta[doc[n], k]) + log(phi[k, w[n]]);
    target += log_sum_exp(gamma);  // likelihood;
  }
}
"

data = list(K = 5,
            V = length(unique(datalong$Mutation)),
            M = length(unique(datalong$doc)),
            N = sum(datalong$value),
            w = rep(datalong$Mutation,datalong$value),
            doc = rep(datalong$doc,datalong$value),
            alpha = rep(50/5,5),
            beta = rep(1,length(unique(datalong$Mutation))))

stan.model <- stan_model(model_code = lda)
stan.model.fit = sampling(stan.model,
                          data=data,
                          iter=1000,
                          chains=5,
                          seed=2019)


vals2 = as.matrix(stan.model.fit)
dim(vals2)
displayvalues=data.frame(names=rep(NA,ncol(vals2)),value=rep(NA,ncol(vals2)),madsd=rep(NA,ncol(vals2)),dim1=rep(NA,ncol(vals2)),dim2=rep(NA,ncol(vals2)),variable=rep(NA,ncol(vals2)))
for(i in 1:ncol(vals2)){
  displayvalues$names[i] = colnames(vals2)[i]
  displayvalues$value[i] = median(vals2[,i])
  displayvalues$madsd[i] = mad(vals2[,i])
  displayvalues$dim1[i] = gsub("[^0-9.-]","",strsplit(colnames(vals2)[i],",")[[1]])[1]
  displayvalues$dim2[i] = gsub("[^0-9.-]","",strsplit(colnames(vals2)[i],",")[[1]])[2]
  displayvalues$variable[i] = gsub("[^a-zA-Z]","",colnames(vals2)[i])
}
displayvalues$dim1 = as.numeric(displayvalues$dim1)
displayvalues$dim2 = as.numeric(displayvalues$dim2)

phimat = dcast(displayvalues %>% filter(variable=="phi"), dim1 ~ dim2)[,-1]
colnames(phimat) = paste("doc",1:ncol(phimat),sep="")
rownames(phimat) = paste("topic",1:nrow(phimat),sep="")
displayphi = data.frame(t(phimat),names=colnames(phimat))

p1 = ggplot(displayphi) + aes(x=names,y=topic1) + geom_bar(stat="identity") + theme(axis.text.x = element_blank(),legend.position="none") + xlab("") + ylab("Signature 1")
p2 = ggplot(displayphi) + aes(x=names,y=topic2) + geom_bar(stat="identity") + theme(axis.text.x = element_blank(),legend.position="none") + xlab("") + ylab("Signature 2")
p3 = ggplot(displayphi) + aes(x=names,y=topic3) + geom_bar(stat="identity") + theme(axis.text.x = element_blank(),legend.position="none") + xlab("") + ylab("Signature 3")
p4 = ggplot(displayphi) + aes(x=names,y=topic4) + geom_bar(stat="identity") + theme(axis.text.x = element_blank(),legend.position="none") + xlab("") + ylab("Signature 4")
p5 = ggplot(displayphi) + aes(x=names,y=topic5) + geom_bar(stat="identity") + theme(axis.text.x = element_blank(),legend.position="none") + xlab("") + ylab("Signature 5")



grid.arrange(p1,p2,p3,p4,p5,ncol=1)

