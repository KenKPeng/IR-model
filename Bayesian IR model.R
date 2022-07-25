# load package
library(R2jags)

# Read the data
BC = read.csv("BC_data.csv", sep = " ")
head(BC)
y <- BC$performance[which(!is.na(BC$performance))]
y # y is all the performance scores
train_onday = BC$training_load
train_onday # train_onday is the daily training load
day = which(BC$performance!=0)-1
day # day is the day numbers where performance is available

# Bayesian IR model

#### re para
############ Model fitting #############

# Define the parameters and the data:
n = length(y) # Draw the number of days where the performance data is available
data = list("n" =n,"y"=y,"train_onday"=train_onday,"day"=day) # identify data

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


# Initialize the MCMC
jags.m <- jags.model( file = "modelBC2018_constrain1.txt", data=data, inits=inits, n.chains=2, n.adapt = 1000)

# Draw the sample from MCMC
post_sample = coda.samples(jags.m, 
                           c("p0", "k1","theta","t1","t2", "rho","sigma", "tn", "tg"), # choose which variables are needed to be monitored
                           n.iter = 5000) # number of iterations 
# check trace plot
plot(post_sample)
gelman.diag(post_sample)
gelman.plot(post_sample)

# Inference
summary(window(post_sample, start=3001,end = 6000)) # burn-in. Remove  first 1000 adaptation iterations and first 2000 regular iterations




















####################### Optional, extra step to implement the constrains ################################
####################### Only needed when specifying multivariate normal priors ##########################

# Create a empty matrix to save the iterations
Bay_sample = matrix(, nrow = 0, ncol = 7)
colnames(Bay_sample) = c("k1","theta","p0","rho","sigma","t1","t2")
iter = 0 # iteration number


while (iter <= 1000){ #iteration number
  tryCatch({
    # Sampling from the model:
    jags.m <- jags.model( file = "modelBC2018_constrain.txt", data=data, inits=inits, n.chains=1, n.adapt = 0)
    
    post_sample = coda.samples(jags.m, 
                               c("p0", "k1","t","sigma","rho", "theta"), # choose which variables are needed to be monitored
                               n.iter = 1) # number of iterations 
    samp = window(post_sample)
    k1 = unlist(samp[,1])
    theta = unlist(samp[,2])
    t1 = unlist(samp[,6])
    t2 = unlist(samp[,7])
    p0 = unlist(samp[,3])
    rho = unlist(samp[,4])
    sigma = unlist(samp[,5])
    
    if (t1>t2 & t1>5 & t2 >3 & t1<60 & t2<60){ # If satisfy constraints, iteration + 1 and save the result, and update the initial value . Otherwise, re-sample it.
      inits <- function () {list(p0 = p0,
                                 k = c(k1,theta),
                                 t = c(t1,t2),
                                 sigma = sigma,
                                 rho = rho)}
      Bay_sample = rbind(Bay_sample,c(k1,theta,p0,rho,sigma,t1,t2))
      iter = iter + 1
    }
  }, error=function(e){})
}
