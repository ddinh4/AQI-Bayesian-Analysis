library(MCMCpack)
df1 = read.csv("data_date.csv")
unique_data <- df1 %>% distinct(Country, .keep_all = TRUE)
aqi_2 <- unique_data[,c("Status","AQI_Value")]
aqi <- aqi_2[,2]
aqi <-  aqi[aqi <= 500]
aqi_df <- data.frame(aqi)
a = 5
log.posterior.mu <- function(mu,sigma2,mu0,tau2,y)
  #sum(dnorm(y, mu, sqrt(sigma2), log=TRUE)) + dnorm(mu,mu0, sqrt(tau2),log=TRUE)
  sum(dgamma(y, shape = (mu^2)/((mu^2)/a), rate = mu/((mu^2)/a), log=TRUE)) + dnorm(mu,mu0, sqrt(tau2),log=TRUE) 

log.posterior.sigma2 <- function(sigma2,mu,alpha,beta,y) {
  if (sigma2<0) return(-Inf) else
    return(sum(dgamma(y, shape = (mu^2)/((mu^2)/a), rate = mu/((mu^2)/a), log=TRUE)) + log(dinvgamma(sigma2,alpha,beta)))
}

y <- aqi #data

## prior hyperparameters
mu0 <- 70
tau2 <- (70^2)/5
alpha <- 0.1
beta <- 0.1


sigma2 <- var(aqi) # initial value for sigma2
mu<- mean(aqi) # initial value for mu 
n.sim<-10000 ; mu.mcmc<-NULL ; sigma2.mcmc <- NULL; accept.mu <- 0;accept.sigma2 <- 0
delta1 <- 40 ; ## normal proposal var for mu
delta2 <- 6900 ## normal proposal var for sigma2

for ( ite in 1 : n.sim) {
  
  mu.star<-rnorm(1, mu, sqrt(delta1))
  
  log.r <- log.posterior.mu(mu.star,sigma2,mu0,tau2,y) - log.posterior.mu(mu,sigma2,mu0,tau2,y)
  
  if (log(runif(1))< log.r) { 
    mu <- mu.star 
    accept.mu <- accept.mu+1
  }
  mu.mcmc<-c(mu.mcmc, mu)
  
  sigma2.star<-rnorm(1, sigma2, sqrt(delta2) )
  
  log.r <- log.posterior.sigma2(sigma2.star,mu,alpha,beta,y) - log.posterior.sigma2(sigma2,mu,alpha,beta,y)
  
  if (log(runif(1))< log.r) { 
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

#when change to gamma prior for both, the chains are less extreme.