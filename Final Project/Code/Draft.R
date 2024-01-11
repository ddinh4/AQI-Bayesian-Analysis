library(ggplot2)
library(dplyr)
library(MCMCpack)
df1 = read.csv("data_date.csv")
unique_data <- df1 %>% distinct(Country, .keep_all = TRUE)
aqi_2 <- unique_data[,c("Status","AQI_Value")]
aqi <- aqi_2[,2]
aqi <-  aqi[aqi <= 500]
aqi_df <- data.frame(aqi)
#Comparing two models
hist(aqi)
#Normal model for data with unknown mu and sigma^2
mu = mean(aqi)
sigma2 = var(aqi)
likelihood <- function(aqi, mu, sigma2)
  prod(dnorm(aqi, mu, sqrt(sigma2)))
       
a = 5
x <- seq(0, 50, by = 0.1) 
mu0 <- dnorm(x, mean = 70, sd = sqrt((70^2)/5))
sigma02 <- dinvgamma(x, 0.1, 0.1)
likelihood_M1 <- integrate(likelihood(aqi,mu,sigma2)*mu0*sigma02, 0, Inf)


#MCMC for Gamma data with unknown mu and sigma
posterior.mu <- function(mu,sigma2,mu0,tau2,y)
  prod(dgamma(y, shape = (mu^2)/sigma2, rate = mu/sigma2))*dnorm(mu,mu0,sqrt(tau2))

posterior.sigma2 <- function(sigma2, mu, alpha, beta, y) {
  if (sigma2 < 0) return(-Inf) else
    return(prod(dgamma(y, shape = (mu^2)/sigma2, rate = mu/sigma2))*dinvgamma(sigma2,alpha,beta))
}

y = aqi

mu0 = 70
tau2 = (70^2)/5
alpha = 0.1
beta = 0.1

sigma2 <- 5 # initial value for sigma2
mu<-64 # initial value for mu 
n.sim<-10000 ; mu.mcmc<-NULL ; sigma2.mcmc <- NULL; accept.mu <- 0;accept.sigma2 <- 0
delta1 <- 1.5 ; ## normal proposal var for mu
delta2 <- 10 ## normal proposal var for sigma2

for ( ite in 1 : n.sim) {
  
  mu.star<-rnorm(1, mu, sqrt(delta1)) 
  
  r <- posterior.mu(mu.star,sigma2,mu0,tau2,y)/posterior.mu(mu,sigma2,mu0,tau2,y)
  
  if ((runif(1))< r) { 
    mu <- mu.star 
    accept.mu <- accept.mu+1
  }
  mu.mcmc<-c(mu.mcmc, mu)
  
  sigma2.star<-rnorm(1, sigma2, sqrt(delta2) )
  
  r <- posterior.sigma2(sigma2.star,mu,alpha,beta,y)/posterior.sigma2(sigma2,mu,alpha,beta,y)
  
  if ((runif(1))< r) { 
    sigma2 <- sigma2.star 
    accept.sigma2 <- accept.sigma2+1
  }
  sigma2.mcmc<-c(sigma2.mcmc, sigma2)
  
  
}

accept.mu/n.sim
accept.sigma2/n.sim

par(mfrow=c(2,1))
plot(mu.mcmc,type='l')
plot(sigma2.mcmc,type='l')

mean(mu.mcmc[1000:10000])
mean(sigma2.mcmc[1000:10000])
var(mu.mcmc[1000:10000])
var(sigma2.mcmc[1000:10000])
