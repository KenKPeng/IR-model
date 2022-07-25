library(MASS)
library(ggplot2)
library(matrixStats)
library(reshape2)
library(ggthemes)
library(ggplot2)
library(gridExtra)
library(hrbrthemes)

# MCMC result
Bay_sample = read.csv("C:/Users/Ken Peng/Desktop/IR model coding/sampling_results/2018 season/BC_sample_TRIMPi_2018_repara.csv",sep = " ")


y = matrix(, nrow = 365, ncol = 1000)

for (j in 1:1000){
  p0 = Bay_sample$p0[j+1]
  k1 = Bay_sample$k1[j+1]
  k2 =  Bay_sample$k2[j+1]
  t1 = Bay_sample$t1[j+1]
  t2 =   Bay_sample$t2[j+1]
  
  train_load = p0
  for (i in 1:(dim(BC)[1]-1)){
    day = 0:i
    train_onday = BC[1:(i+1),"training_load"]
    train_load = c(train_load,p0+sum((k1*exp(-(i-day)/t1)-k2*exp(-(i-day)/t2))*train_onday))
  }
  
  sigma = matrix(, nrow = 365, ncol = 365)
  pp = Bay_sample$pp[j+1]
  ss = Bay_sample$ss[j+1]
  for (a in 1:365){
    for (b in 1:365){
      sigma[a,b] = sqrt((ss^2)*(pp^abs(a-b)))
    }
  }
  y[,j] = mvrnorm(1,train_load,sigma)
}


pi = rowQuantiles(y, probs = c(0.025,0.975))

p0 = mean(Bay_sample$p0)
k1 = mean(Bay_sample$k1)
k2 = mean(Bay_sample$k2)
t1 = mean(Bay_sample$t1)
t2 = mean(Bay_sample$t2)
train_load = p0
for (i in 1:(dim(BC)[1]-1)){
  day = 0:(i-1)
  train_onday = BC[1:(i+1),"training_load"]
  train_load = c(train_load,p0+sum((k1*exp(-(i-day)/t1)-k2*exp(-(i-day)/t2))*train_onday))
}
BC$p  =train_load



p1 <- ggplot(BC, aes(x=Day, y=p)) +
  geom_line(color="#69b3a2",size = 1.2) + 
  xlab("") +
  theme_minimal()+
  theme(text = element_text(size=18), axis.text.x=element_text(angle=0, hjust=1)) + ylim(800,1150)+
  geom_point(data = BC, size = 3,
             mapping = aes(x = Day, y = performance, color = "red"))+ theme(legend.position = "none")+
  labs(subtitle = "Bivariate normal prior on t1, t2. Based on Hellard et al")+ylab("")+xlab("Day")+
  geom_ribbon(aes(ymin=pi[,1], ymax=pi[,2]), fill = "#69b3a2", alpha=0.4)
p1







