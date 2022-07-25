############ Data preparation #############
# 3 variables are needed: Daily total training load, performance, Day. 
# Following is an example data:
BC = read.csv("BC_data.csv", sep = " ")
head(BC)
############ Data preparation #############


############ Defining the nonlinear least square function of IR model ############ 
fn <- function(theta) {
  k1 = theta[1]
  k2 = theta[2]
  t1 = theta[3]
  t2 = theta[4]
  p0= theta[5]
  train_load = c()
  train_load_data = BC$performance[which(!is.na(BC$performance))]
  for (i in c(which(!is.na(BC$performance))-1)){
    day = 0:i
    train_onday = BC[1:(i+1),"training_load"]
    train_load = c(train_load,p0+sum((k1*exp(-(i-day)/t1)-k2*exp(-(i-day)/t2))*train_onday)) # Here is the Banister IR model
  }
  sum((train_load - train_load_data)^2) # Here is the square of the difference of calculated performance and true performance, we are trying to minimize this value
}
# Confidence interval can be calcluated by bootstrap
############ Defining the nonlinear least square function of IR model ############ 



############ use PORT routines to minimize this function ############ 
est1 = nlminb(start = c(0.036,0.050,38,19,1050), # starting value
              fn)
est1









############ Maximum likelihood Estimation with normal assumption############ 
fn <- function(theta) {
  k1 = theta[1]
  k2 = theta[2]
  t1 = theta[3]
  t2 = theta[4]
  s = theta[5]  # Another parameter(variance) is needed here since I suppose the performance is normal distributed
  p0= theta[6]
  train_load = c()
  for (i in c(which(!is.na(BC$performance))-1)){
    day = 0:i
    train_onday = BC[1:(i+1),"training_load"]
    train_load = c(train_load,p0+sum((k1*exp(-(i-day)/t1)-k2*exp(-(i-day)/t2))*train_onday))# Banister IR model
  }
  train_load_data = BC$performance[which(!is.na(BC$performance))] # Draw avaliable performance 
  sum(0.5*((train_load_data - train_load)/s)^2 + 0.5*log(s)) # Negative log likelihood function. To maximize log Likelihood, we want to minimize this one.
}


############ use PORT routines to minimize the Negative likelihood function ############ 
mle_est = nlminb(start = c(0.036,0.050,38,19,5,50),fn)
mle_est
