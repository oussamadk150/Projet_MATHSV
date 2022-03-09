## ----librariries-----------------------------------------------------
library(tidyverse) # Pour la manipulation de données
library(MASS) 

## ----simulation de Y-------------------------------------------------
# On fixe les dimensions 
n = 100
p = 10
k_etoile = 2

# Pour avoir les mêmes données
set.seed(123)

# On choisit sigma
sigma = diag(10**-2, p)
# On simule les epsilon
epsilon = mvrnorm(n, rep(0,p), sigma)

#On fixe lambda
lambda = matrix(sample(-1:1, 
                       size = p*k_etoile, 
                       replace = T),
                nrow = p, ncol = k_etoile)

#On simule les eta
eta = mvrnorm(n, rep(0,k_etoile), diag(1,k_etoile))

q <- 3

X = mvrnorm(n, rep(0, q), diag(1, q))

beta <- matrix(sample(-1:1, 
              size = p * q, 
              replace = TRUE),
       nrow = q, ncol = p)

#On agrège pour obtenir Y
Y = X %*% beta + eta %*% t(lambda) + epsilon