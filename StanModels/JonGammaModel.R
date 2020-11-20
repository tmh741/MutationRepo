library(rstan)

data = readRDS("lung_sbs_matrix.rds")
counts = as.matrix(data[,1])

jonathanmodel = "
data {
    int<lower=1> I;
    int<lower=1> J;
    int<lower=1> K;
    int<lower=0> X[I, J];
    real<lower=0> a;
    real<lower=0> alpha;
    real<lower=0, upper=1> lik_power;
}

transformed data {
    real<lower=0> rho[J]; 
    real<lower=0> total_count = 0;
    real<lower=0> scale;
    real<lower=0> a0;
    real<lower=0> b0;
    vector<lower=0>[I] alpha_array = rep_vector(alpha, I);
    for (j in 1:J)
      rho[j] = sum(X[,j]);
      
    for (i in 1:I)
        for (j in 1:J)
            total_count += X[i, j];
    scale = total_count / mean(rho) / J;
    a0 = 1./K;
    b0 = 1./scale;
}

parameters {
    matrix<lower=0>[K, J] theta;
    simplex[I] r[K];
    vector<lower=0>[K] mu;
}

model {
    real mutation_rate;

    for (k in 1:K) {
        mu[k] ~ gamma(a0, b0);
        r[k] ~ dirichlet(alpha_array);
    }
    for (j in 1:J)
        for (k in 1:K)
            theta[k, j] ~ gamma(a, a / mu[k]);

    for (i in 1:I)
        for (j in 1:J) {
            mutation_rate = 0;
            for (k in 1:K)
                mutation_rate += r[k,i]*theta[k,j];
            if (lik_power == 1)
                X[i, j] ~ poisson(rho[j] * mutation_rate);
            else {
                if (X[i, j] > 0)
                    target += lik_power * X[i, j] * log(mutation_rate);
                target += -lik_power * rho[j] * mutation_rate;
            }
        }
}
"

a=1000
alpha=1/K
a0=2
b0=4
lik_power=1
I = dim(counts)[1]
J = dim(counts)[2]

test1 = rstan::stan(model_code=jonathanmodel,
                      data = list(I=I, J=J, K=K, X=counts, a=a, alpha=alpha,lik_power=lik_power),
                    iter=10000)

stan_rhat(test1)

vals2= as.matrix(test1)
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

rmat = dcast(displayvalues %>% filter(variable=="r"), dim1 ~ dim2)[,-1]

colnames(rmat) = rownames(counts)
rownames(rmat) = paste("doc",1:5,sep="")
displayr = data.frame(t(rmat))
displayr$names = colnames(rmat)
for(i in 1:nrow(displayr)){
displayr$alex[i] = strsplit(displayr$names[i],"_",perl=T)[[1]][1]
}
p1 = ggplot(displayr) + aes(x=names,y=doc1,fill=as.character(alex)) + geom_bar(stat="identity") + theme(axis.text.x = element_blank(),legend.position="none") + xlab("") + ylab("Signature 1")
p2 = ggplot(displayr) + aes(x=names,y=doc2,fill=as.character(alex)) + geom_bar(stat="identity") + theme(axis.text.x = element_blank(),legend.position="none") + xlab("") + ylab("Signature 2")
p3 = ggplot(displayr) + aes(x=names,y=doc3,fill=as.character(alex)) + geom_bar(stat="identity") + theme(axis.text.x = element_blank(),legend.position="none") + xlab("") + ylab("Signature 3")
p4 = ggplot(displayr) + aes(x=names,y=doc4,fill=as.character(alex)) + geom_bar(stat="identity") + theme(axis.text.x = element_blank(),legend.position="none") + xlab("") + ylab("Signature 4")
p5 = ggplot(displayr) + aes(x=names,y=doc5,fill=as.character(alex)) + geom_bar(stat="identity") + theme(axis.text.x = element_blank(),legend.position="none") + xlab("") + ylab("Signature 5")

grid.arrange(p1,p2,p3,p4,p5,ncol=1)

