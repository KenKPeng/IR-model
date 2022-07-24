### Multivariate normal

library(MAsigma)

############ Model fitting #############

# Define the parameters and the data:
n = length(y) # how many days where the performance data is available
data = list("n" =n,"y"=y,"train_onday"=train_onday,"day"=day) # identify data

params = c("p0", "k","t","sigma","rho") # identify parameters

# Define the model:
cat(paste("model
            {
            y[1:n] ~ dmnorm(mu[1:n],s[,])
            for (i in 1:n) {
            for (j in 1:",dim(BC)[1],"){
            m[i,j] <- (k[1]*exp(-(day[i]-j+1)/t[1])-k[2]*exp(-(day[i]-j+1)/t[2]))*train_onday[j]
            }
            for (ii in 1:n){
            s1[i,ii] <- (1/tau)*(pow(rho,abs(day[i]-day[ii])))
            }
            mu[i] <- p0 + sum(m[i,1:day[i]])
            }
            s[1:n,1:n] <- inverse(s1[,])
            # setting the priors:
            p0 ~ dnorm(1029,0.0025)
            k[1:2] ~ dmnorm(kmu[1:2],ks[,])
            kmu[1] <- 0.036
            kmu[2] <- 0.05
            ks1[1,1] <- 0.006^2  
            ks1[1,2] <- 0.0000216
            ks1[2,1] <- 0.0000216
            ks1[2,2] <- 0.004^2
            ks[1:2,1:2] <- inverse(ks1[,]) ## jags take the inverse of the variance
            t[1:2] ~ dmnorm(tmu[1:2],ts[,])
            tmu[1] <- 50
            tmu[2] <- 17
            ts1[1,1] <- 38^2
            ts1[1,2] <- 182
            ts1[2,1] <- 182
            ts1[2,2] <- 12^2
            ts[1:2,1:2] <- inverse(ts1[,])
            tau ~ dgamma(1.0E-3, 1.0E-3) # jags doesn't allow the improper prior. This make p(sigma) approximately = 1/sigma
            sigma <- 1/sqrt(tau)
            rho ~ dbeta(10,1)
            }",sep = ""), file="modelBC2018_noconstrain.txt")


# Set the random starting values for each chain:
inits <- function () {list(p0 = rnorm(1,mean = 1029,sd = 20),
                           k = mvrnorm(1,c(0.036,0.05),matrix(c(0.012^2,0.00017,0.00017,0.016^2), nrow = 2)),
                           t = mvrnorm(1,c(50,40),matrix(c(5^2,13.5,13.5,3^2), nrow = 2)),
                           tau = 0.0025,
                           rho = rbeta(1, 10, 1))
  }

## Fixed initial values
## inits <- function () {list(p0 = 1029,
##                           k = c(0.036,0.05),
##                           t = c(38,19),
##                           sigma = 20,
##                           rho = 0.9)}



jags.m <- jags.model(file = "modelBC2018_noconstrain.txt", data=data, inits=inits, n.chains=2, n.adapt = 1000)



post_sample = coda.samples(jags.m, 
                           c("p0", "k","t","sigma","rho"), # choose which variables are needed to be monitored
                           n.iter = 5000) # number of iterations 
head(post_sample)
plot(post_sample)
gelman.diag(post_sample)
gelman.plot(post_sample)
##


## burn-in first 2000 iterations
summary(window(post_sample, start=1001,end = 5000))


## Set Constraints:
# Create a empty matrix to save the iterations
Bay_sample = matrix(, nrow = 0, ncol = 7)
colnames(Bay_sample) = c("k1","k2","p0","rho","sigma","t1","t2")
iter = 0 # iteration number


# Procesigma
while (iter <= 1000){
  tryCatch({
  # Sampling from the model:
    jags.m <- jags.model( file = "modelBC2018_noconstrain.txt", data=data, inits=inits, n.chains=1, n.adapt = 0)
    
  post_sample = coda.samples(jags.m, 
                             c("p0", "k","t","sigma","rho"), # choose which variables are needed to be monitored
                             n.iter = 1) # number of iterations 
  samp = window(post_sample)
  k1 = unlist(samp[,1])
  k2 = unlist(samp[,2])
  t1 = unlist(samp[,6])
  t2 = unlist(samp[,7])
  p0 = unlist(samp[,3])
  rho = unlist(samp[,4])
  sigma = unlist(samp[,5])
  
  if (k1<k2 & t1>t2 & k1>0 & t2>0){ # If satisfy constraints, iteration + 1 and save the result, and update the initial value . Otherwise, re-sample it.
    inits <- function () {list(p0 = p0,
                               k = c(k1,k2),
                               t = c(t1,t2),
                               sigma = sigma,
                               rho = rho)}
    Bay_sample = rbind(Bay_sample,c(k1,k2,p0,rho,sigma,t1,t2))
    iter = iter + 1
  }
  }, error=function(e){})
}

write.table(Bay_sample,file = "D:/projects/IR model coding/sampling_results/2018 season/FF_sample_TRIMPi_2018.csv")



######################### Results
Bay_sample = read.csv("D:/projects/IR model coding/sampling_results/2018 season/BC_sample_TRIMPi_2018_repara.csv",sep = " ")
summary(Bay_sample)

# A function to summarize the results:
IR_result = function(Bay_sample){
########## Point estimate
est = data.frame(k1 = c(mean(Bay_sample[,1]),sd(Bay_sample[,1])),
                 k2 = c(mean(Bay_sample[,2]),sd(Bay_sample[,2])),
                 p0 = c(mean(Bay_sample[,3]),sd(Bay_sample[,3])),
                 rho = c(mean(Bay_sample[,4]),sd(Bay_sample[,4])),
                 sigma = c(mean(Bay_sample[,5]),sd(Bay_sample[,5])),
                 t1 = c(mean(Bay_sample[,6]),sd(Bay_sample[,6])),
                 t2 = c(mean(Bay_sample[,7]),sd(Bay_sample[,7])),
                 row.names = c("mean","sd"))


######### Get interval estimate
quantiles = data.frame(k1 = c(quantile(Bay_sample[,1],probs = c(0.025,0.5,0.975))),
                       k2 = c(quantile(Bay_sample[,2],probs = c(0.025,0.5,0.975))),
                       p0 = c(quantile(Bay_sample[,3],probs = c(0.025,0.5,0.975))),
                       rho = c(quantile(Bay_sample[,4],probs = c(0.025,0.5,0.975))),
                       sigma = c(quantile(Bay_sample[,5],probs = c(0.025,0.5,0.975))),
                       t1 = c(quantile(Bay_sample[,6],probs = c(0.025,0.5,0.975))),
                       t2 = c(quantile(Bay_sample[,7],probs = c(0.025,0.5,0.975))),
                       row.names = c("2.5%","50%","97.5%"))
##### tn and tg
value1 = ((est$t1[1]*est$t2[1])/(est$t1[1]-est$t2[1]))*log((est$k2[1]/est$k1[1])*(est$t1[1]/est$t2[1])) 
peak1 = ((Bay_sample[,6]*Bay_sample[,7])/(Bay_sample[,6]-Bay_sample[,7]))*log((Bay_sample[,2]/Bay_sample[,1])*(Bay_sample[,6]/Bay_sample[,7]))  #### tg
tg = c(quantile(peak1,probs = c(0.025,0.5,0.975)),mean = mean(peak1),sd = sd(peak1),`model estimated value` = round(value1))

value2 = ((est$t1[1]*est$t2[1])/(est$t1[1]-est$t2[1]))*log((est$k2[1]/est$k1[1]))
rec1 = ((Bay_sample[,6]*Bay_sample[,7])/(Bay_sample[,6]-Bay_sample[,7]))*log((Bay_sample[,2]/Bay_sample[,1]))      ##tn
tn = c(quantile(rec1,probs = c(0.025,0.5,0.975)),mean = mean(rec1),sd = sd(rec1),`model estimated value` = round(value2))


est = list(`point estimate` = t(est),`Interval estimate` = t(quantiles), tg = tg, tn = tn)
print(est)
}

IR_result(Bay_sample)





## plot
ts.plot(Bay_sample[,7])
plot(density(Bay_sample[,3]))

sd((Bay_sample[,3]))







#### re para
############ Model fitting #############

# Define the parameters and the data:
n = length(y) # how many days where the performance data is available

data = list("n" =n,"y"=y,"train_onday"=train_onday,"day"=day) # identify data

params = c("p0", "k1","t","sigma","rho","k2") # identify parameters

# Define the model:
cat(paste("model
            {
            y[1:n] ~ dmnorm(mu[1:n],s[,])
            for (i in 1:n) {
            for (j in 1:",dim(BC)[1],"){
            m[i,j] <- (k1*exp(-(day[i]-j+1)/t[1])-k1*theta*exp(-(day[i]-j+1)/t[2]))*train_onday[j]
            }
            for (ii in 1:n){
            s1[i,ii] <- (1/tau)*(pow(rho,abs(day[i]-day[ii])))
            }
            mu[i] <- p0 + sum(m[i,1:day[i]])
            }
            s[1:n,1:n] <- inverse(s1[,])
            # setting the priors:
            p0 ~ dnorm(1000,0.0025)
            k1 ~ dunif(0,10)
            theta ~ dnorm(4.137,0.4)I(1,)
            t[1:2] ~ dmnorm(tmu[1:2],ts[,])
            tmu[1] <- 50
            tmu[2] <- 13
            ts1[1,1] <- 38^2
            ts1[1,2] <- 182
            ts1[2,1] <- 182
            ts1[2,2] <- 12^2
            ts[1:2,1:2] <- inverse(ts1[,])
            tau ~ dgamma(1.0E-3, 1.0E-3) # jags doesn't allow the improper prior. This make p(sigma) approximately = 1/sigma
            sigma <- 1/sqrt(tau)
            rho ~ dbeta(10,1)
            k2 <- theta*k1
            }",sep = ""), file="modelBC2018_constrain.txt")


# Set the random starting values for each chain:
inits <- function () {list(p0 = rnorm(1,mean = 1029,sd = 20),
                           k1 = 0.5,
                           theta = 4.127,
                           t = mvrnorm(1,c(50,13),matrix(c(5^2,13,13,3^2), nrow = 2)),
                           tau = 0.0025,
                           rho = rbeta(1, 10, 1))
}

## Fixed initial values
inits <- function () {list(p0 = 1029,
k1 = 0.5,
theta = 2.259,
t = c(50,17),
tau = 0.0025,
rho = 0.9)}



jags.m <- jags.model( file = "modelBC2018_constrain.txt", data=data, inits=inits, n.chains=2, n.adapt = 10)

post_sample = coda.samples(jags.m, 
                           c("p0", "k1","t","sigma","rho", "k2"), # choose which variables are needed to be monitored
                           n.iter = 5) # number of iterations 
head(post_sample)
plot(post_sample)
summary(window(post_sample, start=1001,end = 5000))
gelman.diag(post_sample)
gelman.plot(post_sample)
##

write.table(Bay_sample,file = "D:/projects/IR model coding/sampling_results/2018 season/FF_sample_TRIMPi_2018.csv")


Bay_sample = read.table("D:/projects/IR model coding/sampling_results/2018 season/FF_sample_TRIMPi_2018.csv")
summary(Bay_sample)

## with constraints, complex
# Create a empty matrix to save the iterations
Bay_sample = matrix(, nrow = 0, ncol = 7)
colnames(Bay_sample) = c("k1","k2","p0","rho","sigma","t1","t2")
iter = 0 # iteration number


while (iter <= 1000){
  tryCatch({
    # Sampling from the model:
    jags.m <- jags.model( file = "modelBC2018_constrain.txt", data=data, inits=inits, n.chains=1, n.adapt = 0)
    
    post_sample = coda.samples(jags.m, 
                               c("p0", "k1","t","sigma","rho", "k2"), # choose which variables are needed to be monitored
                               n.iter = 1) # number of iterations 
    samp = window(post_sample)
    k1 = unlist(samp[,1])
    k2 = unlist(samp[,2])
    t1 = unlist(samp[,6])
    t2 = unlist(samp[,7])
    p0 = unlist(samp[,3])
    rho = unlist(samp[,4])
    sigma = unlist(samp[,5])
    
    if (t1>t2 & t1>5 & t2 >3 & t1<60 & t2<60){ # If satisfy constraints, iteration + 1 and save the result, and update the initial value . Otherwise, re-sample it.
      inits <- function () {list(p0 = p0,
                                 k = c(k1,k2),
                                 t = c(t1,t2),
                                 sigma = sigma,
                                 rho = rho)}
      Bay_sample = rbind(Bay_sample,c(k1,k2,p0,rho,sigma,t1,t2))
      iter = iter + 1
    }
  }, error=function(e){})
}






######### How about no correlation with t1 t2


#### re para
############ Model fitting #############

# Define the parameters and the data:
n = length(y) # how many days where the performance data is available

data = list("n" =n,"y"=y,"train_onday"=train_onday,"day"=day) # identify data

params = c("p0", "k1","t1","t2","sigma","rho","k2") # identify parameters

# Define the model:
cat(paste("model
            {
            y[1:n] ~ dmnorm(mu[1:n],s[,])
            for (i in 1:n) {
            for (j in 1:",dim(BC)[1],"){
            m[i,j] <- (k1*exp(-(day[i]-j+1)/t1)-k1*theta*exp(-(day[i]-j+1)/t2))*train_onday[j]
            }
            for (ii in 1:n){
            s1[i,ii] <- (1/tau)*(pow(rho,abs(day[i]-day[ii])))
            }
            mu[i] <- p0 + sum(m[i,1:day[i]])
            }
            s[1:n,1:n] <- inverse(s1[,])
            # setting the priors:
            p0 ~ dnorm(1000,0.0025)
            k1 ~ dunif(0,10)
            theta ~ dnorm(4.137,0.4)I(1,)
            t1 ~ dnorm(50,0.16)I(5,60)
            t2 ~ dnorm(13,0.29)I(3,60)
            tau ~ dgamma(1.0E-3, 1.0E-3) # jags doesn't allow the improper prior. This make p(sigma) approximately = 1/sigma
            sigma <- 1/sqrt(tau)
            rho ~ dbeta(10,1)
            k2 <- theta*k1
            tn = ((t1*t2)/(t1-t2))*log(k2/k1)
            tg = ((t1*t2)/(t1-t2))*log((k2/k1)*(t1/t2))
            }",sep = ""), file="modelBC2018_constrain1.txt")


# Set the random starting values for each chain:
inits <- function () {list(p0 = rnorm(1,mean = 1029,sd = 20),
                           k1 = 0.5,
                           theta = 2.25,
                           t1 = rnorm(1,mean = 50, sd = 10),
                           t2 = rnorm(1,mean = 13, sd = 10),
                           tau = 0.0025,
                           rho = rbeta(1, 10, 1))
}

## Fixed initial values
inits <- function () {list(p0 = 1029,
                           k1 = 0.5,
                           theta = 2.259,
                           t1 = 50,
                           t2 = 17,
                           tau = 0.0025,
                           rho = 0.9)}



jags.m <- jags.model( file = "modelBC2018_constrain1.txt", data=data, inits=inits, n.chains=2, n.adapt = 1000)

post_sample = coda.samples(jags.m, 
                           c("p0", "k1","t1","sigma","rho", "k2","t2", "tn", "tg","theta"), # choose which variables are needed to be monitored
                           n.iter = 5000) # number of iterations 
summary(post_sample)
plot(post_sample)
summary(window(post_sample, start=8001,end = 11000))
gelman.diag(post_sample)
gelman.plot(post_sample)
##


samp = window(post_sample)
k1 = unlist(samp[,1])
k2 = unlist(samp[,2])
t1 = unlist(samp[,6])
t2 = unlist(samp[,7])
p0 = unlist(samp[,3])
rho = unlist(samp[,4])
sigma = unlist(samp[,5])
tg = unlist(samp[,8])
tn = unlist(samp[,10])
theta = unlist(samp[,9])
Bay_sample = as.data.frame(cbind(k1,k2,p0,rho,sigma,t1,t2,tn,tg,theta))

write.table(Bay_sample,file = "D:/projects/IR model coding/sampling_results/2018 season/FF_sample_TRIMPi_2018.csv")













## link two univariate normal prior

#### re para
############ Model fitting #############

# Define the parameters and the data:
n = length(y) # how many days where the performance data is available

data = list("n" =n,"y"=y,"train_onday"=train_onday,"day"=day) # identify data

params = c("p0", "k1","t1","sigma","rho","k2","t2") # identify parameters

# Define the model:
cat(paste("model
            {
            y[1:n] ~ dmnorm(mu[1:n],s[,])
            for (i in 1:n) {
            for (j in 1:",dim(BC)[1],"){
            m[i,j] <- (k1*exp(-(day[i]-j+1)/t1)-k1*theta*exp(-(day[i]-j+1)/t2))*train_onday[j]
            }
            for (ii in 1:n){
            s1[i,ii] <- (1/tau)*(pow(rho,abs(day[i]-day[ii])))
            }
            mu[i] <- p0 + sum(m[i,1:day[i]])
            }
            s[1:n,1:n] <- inverse(s1[,])
            # setting the priors:
            p0 ~ dnorm(1000,0.0025)
            k1 ~ dunif(0,10)
            theta ~ dnorm(4.137,0.4)I(1,)
            
            
            
            ts1[1] <- 38
            taut[1] <- 1/ts1[1]^2
            ts1[2] <- 12
            taut[2] <- 1/ts1[2]^2
            
            rho <- 0.9
            var.eta <- taut[2]/(1 - rho^2)
            tmu[1] <- 50
            tmu[2] <- 13
            
            t1 ~ dnorm(tmu[1], taut[1])I(5,60)
            t22 ~ dnorm(tmu[2], var.eta)I(3,60)
            
            t2 = t22 + rho * ts1[2]/ts1[1] * (t1 - tmu[1])
      

            tau ~ dgamma(1.0E-3, 1.0E-3) # jags doesn't allow the improper prior. This make p(sigma) approximately = 1/sigma
            sigma <- 1/sqrt(tau)
            rho ~ dbeta(10,1)
            k2 <- theta*k1
            }",sep = ""), file="modelBC2018_constrain.txt")


# Set the random starting values for each chain:
inits <- function () {list(p0 = rnorm(1,mean = 1029,sd = 20),
                           k1 = 0.5,
                           theta = 4.127,
                           t1 = rnorm(1,mean = 50,sd = 5),
                           t22 = rnorm(1,mean = 13,sd = 5),
                           tau = 0.0025,
                           rho = rbeta(1, 10, 1))
}

## Fixed initial values
#inits <- function () {list(p0 = 1029,
#                           k1 = 0.5,
#                           theta = 2.259,
#                           t = c(50,17),
#                           tau = 0.0025,
#                           rho = 0.9)}



jags.m <- jags.model( file = "modelBC2018_constrain.txt", data=data, inits=inits, n.chains=2, n.adapt = 1000)

post_sample = coda.samples(jags.m, 
                           c("p0", "k1","t1","sigma","rho", "k2","t2"), # choose which variables are needed to be monitored
                           n.iter = 5000) # number of iterations 
head(post_sample)
plot(post_sample)
summary(window(post_sample, start=1001,end = 5000))
gelman.diag(post_sample)
gelman.plot(post_sample)
##

write.table(Bay_sample,file = "D:/projects/IR model coding/sampling_results/2018 season/FF_sample_TRIMPi_2018.csv")


#Bay_sample = window(post_sample, start=1001,end = 5000)

## with constraints, complex
# Create a empty matrix to save the iterations
Bay_sample = matrix(, nrow = 0, ncol = 7)
colnames(Bay_sample) = c("k1","k2","p0","rho","sigma","t1","t2")
iter = 0 # iteration number


while (iter <= 1000){
  tryCatch({
    # Sampling from the model:
    jags.m <- jags.model( file = "modelBC2018_constrain.txt", data=data, inits=inits, n.chains=1, n.adapt = 0)
    
    post_sample = coda.samples(jags.m, 
                               c("p0", "k1","t","sigma","rho", "k2"), # choose which variables are needed to be monitored
                               n.iter = 1) # number of iterations 
    samp = window(post_sample)
    k1 = unlist(samp[,1])
    k2 = unlist(samp[,2])
    t1 = unlist(samp[,6])
    t2 = unlist(samp[,7])
    p0 = unlist(samp[,3])
    rho = unlist(samp[,4])
    sigma = unlist(samp[,5])
    
    if (t1>t2 & t1>5 & t2 >3 & t1<60 & t2<60){ # If satisfy constraints, iteration + 1 and save the result, and update the initial value . Otherwise, re-sample it.
      inits <- function () {list(p0 = p0,
                                 k = c(k1,k2),
                                 t = c(t1,t2),
                                 sigma = sigma,
                                 rho = rho)}
      Bay_sample = rbind(Bay_sample,c(k1,k2,p0,rho,sigma,t1,t2))
      iter = iter + 1
    }
  }, error=function(e){})
}
