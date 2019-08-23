#Function that'll be used to generate Fréchet samples using the inverse transform method
rfrechet <- function(n,lambda,alpha){
  U <- runif(n,0,1)
  x <- ((1/lambda) * log(1/U))^(-1/alpha)
  return(x)
}

#function to simultaneously update both parameters of the distribution
updateFrechet <- function(X,Alpha,b = 1){
  count <- 0 #keeping track of accepted values to calculate the acceptance rate
  lambda <- rgamma(1, length(X), sum(X^(-Alpha))) #This is the Gibbs step of the MCMC
  log_post <- function(Alpha){ #Log posterior distribution of alpha
    N <- length(X)
    post <- (N - 2)*log(Alpha) - Alpha*sum(log(X)) - N*log(sum(X^(-Alpha)))
    return(post)
  }
  
  #Metropolis-Hastings
  proposal <- rgamma(1, b*Alpha, b)
  accept_ratio <- log_post(proposal) + dgamma(Alpha, b*proposal, b, log = TRUE) - 
    dgamma(proposal, b*Alpha, b, log = TRUE) - log_post(Alpha)
  accept_prob <- min(1, exp(accept_ratio))
  
  if(runif(1) < accept_prob){
    accept <- proposal
    count <- 1
  } else {
    accept <- Alpha
  }
  return(list(accept = accept, count = count, lambda = lambda))
}

###########################################################################################

X <- rfrechet(n = 50, lambda = 0.5, alpha = 0.25) #True values are alpha = 0.25 and lambda = 0.5
hist(X)

Iter <- 10000 #Number of iterations for the MCMC

Lambda.out <- rep(NA,Iter) #Vector to store posterior samples of Lambda
Alpha.out <- rep(NA, Iter) #Vector to store posterior samples of Alpha
Count <- rep(0, Iter) #Vector to keep track of accepted values

#Chain starting values
Lambda.out[1] <- 1
Alpha.out[1] <- 2

#MCMC
for(i in 2:Iter){
  Updating <- updateFrechet(X, Alpha.out[i-1])
  Alpha.out[i] <- Updating[[1]]
  Count[i] <- Updating[[2]]
  Lambda.out[i] <- Updating[[3]]
  print(i)
}

#Burn-in phase
burn <- seq(500, Iter, by = 10)

#Ploting samples to see if we were able to recover true parameter values
hist(Alpha.out[burn])
plot(Alpha.out[burn], type = 'l')
abline(h = 0.25, col = 'red')
acf(Alpha.out[burn])
quantile(Alpha.out, probs = c(0.025, 0.5, 0.975)) #True value is well within 95% CI

hist(Lambda.out[burn])
plot(Lambda.out[burn], type = 'l')
abline(h = 0.5, col = 'red')
acf(Lambda.out[burn])
quantile(Lambda.out, probs = c(0.025, 0.5, 0.975)) #True value is well within 95% CI
