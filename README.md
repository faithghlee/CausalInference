# Causality Assignment
In this repository, I would like to show, through the use of R, the importance of variable selection in propensity score models used to estimate causal effects. 

# A few concepts and definitions in Causal Inference 
a) Causality: Given X a covariate and an outcome Y, we are interested in finding a causal connection between X and Y. Simply put, we want to analyse the effect of the response (Y) when X is changed. 

b) Definition of confounder: Say we have another variable Z and Z is a cause of both X and Y (ie. Z causes X and Z causes Y). So to see the causal connection between X and Y, we need to account for confounder Z. 

c) In studying the causal effect of X on Y, we need to remove confounding effects. Naturally, we consider confounders (affect both X and Y) and also other variables that may be either connected to X or Y. However, how important is variable selection? 

d) Further information is included in PDF documentation 

# Aim / Task: 
Show through simulation the importance of variable selection in removing confounding effect
