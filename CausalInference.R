#Causal inference 
library(Matching)
library(plyr)
#To reproduce results
set.seed(1)

#Definition 
n = 2500 #number of observations/subject for 1 replication 
B = 1000 #number of replications 
k = 5 #number of stratums
gamma_0 = 0.5 #true effect

#1 Initialize empty matrices  to store estimated lambda's from each of the 1000 replication 
#col 1 will represent estimates from propensity model 1 and so on. 
gamma_match = matrix(0, nrow = B, ncol = 4)
gamma_ipw = matrix(0, nrow = B, ncol = 4)
gamma_strat = matrix(0, nrow = B, ncol = 4)

#Define functions which we will later use in the loop 

#Matching 

match_function<-function(matched_values){
  to_match<-Match(Y = y, Tr = T, X = matched_values, estimand = "ATE")
  index<- rbind(to_match$index.treated, to_match$index.control)
  outcome<-glm(y[index]~T[index], family = poisson, data = data)
  get_coefficient<-summary(outcome)$coef[2]
  return(get_coefficient)
}


#Inverse Probability Weighting 
#Function below to obtain weights
get_weights<-function(values){
  weight = (1/values)*T + (1-T)/(1-values)
  return(weight)
}

#Use weights to fit poisson model to get estimands 
ipw_function<-function(weights_){
  fit<-glm(y~T, weights = weights_, family = poisson, data = data )
  coef = summary(fit)$coef[2]
}

#Stratification
#1. We define k = 5 stratums. The stratums are split according to the propensity scores 
# from each of propensity score model
#2. Then within each quintile, we fit the poisson model. Thereafter, we take the mean of estimate 
# from each quintile as the estimand for that replication 

stratification<-function(value){
  group<-cut(value, breaks = quantile(value, probs = seq(0,1, by = 0.2)), labels = 1:5, include.lowest = TRUE )
  lambda_est<-mean(as.numeric(dlply(.data = data, .(group),
                                    .fun = function(DF) {
                                      glm(y ~ T, data = DF, family = poisson)$coef[2]
                                    })))
  return(lambda_est)
}

for(j in 1:B){
  #simulate the standard normal random variables as per requested in the question 
  
  x.1<-rnorm(n, 0, 1)
  x.2<-rnorm(n, 0, 1)
  x.3<-rnorm(n, 0 ,1)
  
  #A is the probability for each subject of getting T = 1 
  A<-exp(0.5*x.1+0.75*x.3)/(1+exp(0.5*x.1+0.75*x.3))
  T<-rbinom(n, 1, A) #Assign treatment status 
  lambda.rate <- exp(0.5 + 4/(1+exp(-3*x.1)) + x.2 + gamma_0*T) #
  y<-rpois(n, lambda.rate) #count outcome response
  
  data<-data.frame(y,T,x.1,x.2,x.3)
  #Estimate the propensity score model 
  
  ps.1 <-glm(T~x.1, family = binomial(link = "logit"), data = data)
  ps.2 <-glm(T~x.1 + x.3, family = binomial(link = "logit"), data = data)
  ps.3 <-glm(T~x.1 + x.2, family = binomial(link = "logit"), data = data)
  ps.4 <-glm(T~ x.1 + x.2 + x.3, family = binomial(link = "logit"), data = data)
  
  ps.1values <-ps.1$fitted
  ps.2values <-ps.2$fitted
  ps.3values <-ps.3$fitted
  ps.4values <-ps.4$fitted

  gamma_match[j,1] = match_function(ps.1values)
  gamma_match[j,2] = match_function(ps.2values)
  gamma_match[j,3] = match_function(ps.3values)
  gamma_match[j,4] = match_function(ps.4values)

  #Inverse Probability Weighting
  
  weights.1<-get_weights(ps.1values)
  weights.2<-get_weights(ps.2values)
  weights.3<-get_weights(ps.3values)
  weights.4<-get_weights(ps.4values)
  
  gamma_ipw[j,1]<-ipw_function(weights.1)
  gamma_ipw[j,2]<-ipw_function(weights.2)
  gamma_ipw[j,3]<-ipw_function(weights.3)
  gamma_ipw[j,4]<-ipw_function(weights.4)

  #Stratification 
  gamma_strat[j,1] = stratification(ps.1values)
  gamma_strat[j,2] = stratification(ps.2values)
  gamma_strat[j,3] = stratification(ps.3values)
  gamma_strat[j,4] = stratification(ps.4values)
} 

# Finally, from each method of finding the causal effect from each propensity score model, 
# we like to compute the bias, sd and MSE 

# To simplify things, once again, the use of matrix implementation is very useful 
# Specifically : we create a matrix for each method of estimating causal effects 

bias_sd_mse<-function(x){
  bias = 1/B*sum(x-0.5)
  sd = sqrt((1/(B-1))*sum((x-mean(x))^2))
  mse = (1/B)*sum((x - 0.5)^2)
  return(c(bias,sd,mse))
}

match_effect<-matrix(0, nrow = 3, ncol = 4)
#in match_effect[1,1] represents the bias coming from the first propensity score model
# [2,1] represents the sd  coming from the first propensity score model 
# [3,1] represents mse coming from the first propensity score model using matching 
ipw_effect<-matrix(0, nrow = 3, ncol = 4)
strat_effect<-matrix(0, nrow = 3, ncol = 4)

#computes the bias, sd, mse for each column 

match_effect[,1]<-bias_sd_mse(gamma_match[,1])
match_effect[,2]<-bias_sd_mse(gamma_match[,2])
match_effect[,3]<-bias_sd_mse(gamma_match[,3])
match_effect[,4]<-bias_sd_mse(gamma_match[,4])

ipw_effect[,1]<-bias_sd_mse(gamma_ipw[,1])
ipw_effect[,2]<-bias_sd_mse(gamma_ipw[,2])
ipw_effect[,3]<-bias_sd_mse(gamma_ipw[,3])
ipw_effect[,4]<-bias_sd_mse(gamma_ipw[,4])

strat_effect[,1]<-bias_sd_mse(gamma_strat[,1])
strat_effect[,2]<-bias_sd_mse(gamma_strat[,2])
strat_effect[,3]<-bias_sd_mse(gamma_strat[,3])
strat_effect[,4]<-bias_sd_mse(gamma_strat[,4])



