

## --------------------------------------------------------------------
rm(list = ls()) #Nettoyage de l'environnement

## ----librariries-----------------------------------------------------
library(tidyverse) # Pour la manipulation de données
library(dplyr)
library(MASS)      # Pour mvrnorm
library(ggplot2)   # Pour tracer la distribution des paramètres
source("utils_formatting_functions.R")

## ----simulation de Y-------------------------------------------------
# On fixe les dimensions 
n = 1000
p = 10
k = 2
n_iterations = 100

# Pour avoir les mêmes données

set.seed(123)

# On choisit sigma
sigma = diag(10**-2, p)
# On simule les epsilon
epsilon = t(mvrnorm(n, rep(0,p), sigma))

#On fixe lambda
lambda = matrix(sample(-1:1, 
                       size = p*k, 
                       replace = T),
                nrow = p, ncol = k)

#On simule les eta
eta = mvrnorm(n, rep(0,k), diag(1,k))

#On agrège pour obtenir Y
Y = t(lambda %*% t(eta) + epsilon)


## ----hyperparamètres-----------------------------------------------------------
#Choix projet
a_1=2
a_2=3
#Choix de l'article
nu = 3
a_sigma = 1
b_sigma = 0.3



## ----step_1----------------------------------------------------------
# mise à jour de lambda
step1 <- function(Y,eta,sigma_2,phi,tau){
  lambda_post = matrix(0,p,k)
  for (j in 1:p) {
    Dj = diag(phi[j,]*tau)
    var_post = solve(solve(Dj) + t(eta) %*% eta / sigma_2[j])
    moy_post = var_post %*% t(eta) %*% Y[,j] / sigma_2[j]
    lambda_post[j,] = mvrnorm(1, moy_post, var_post)
  }
  return(lambda_post)
}


## ----step_2----------------------------------------------------------
# mise à jour de sigma_2
step2 <- function(Y,eta,lambda){
  sigma_2_post = rep(0,p)
  a_post = a_sigma + n/2
  for (j in 1:p) {
    b_post = b_sigma + 0.5*sum((Y[,j]-t(lambda[j,])%*%t(eta))**2)
    sigma_2_post[j] = 1/(rgamma(1,a_post,b_post))
  }
  return(sigma_2_post)
}


## ----step_3----------------------------------------------------------
# mise à jour de eta
step3 <- function(Y,lambda,sigma_2){
  sigma = diag(sigma_2)
  eta_post = matrix(0,n,k)
  for (i in 1:n) {
    var_post = solve(diag(1,k) + 
                       t(lambda)%*%solve(sigma)%*%lambda)
    moy_post = var_post %*% t(lambda) %*%solve(sigma)%*%Y[i,]
    eta_post[i,] = mvrnorm(1,moy_post,var_post)
  }
  return(eta_post)
}


## ----step_4----------------------------------------------------------
# mise à jour de eta
step4 <- function(lambda,tau){
  phi_post = matrix(0,p,k)
  for (j in 1:p) {
    phi_post[j,] = rgamma(1,(nu+1)/2, (nu+tau*(lambda[j,])**2)/2)
  }
  return(phi_post)
}


## ----step 5----------------------------------------------------------
# mise à jour de delta
step5 <- function(lambda,phi, delta){
  delta_post = delta # Initialisation à NA
  phi_lambda2 <- phi * lambda^2
  # cumprod donne le produit cumulé
  tau_m1 <- cumprod(delta_post) / delta[1] # Vecteur t^{(1)} de l'article
  omega <- colSums(phi_lambda2) # Vecteur des sommes de colonnes
  delta_post[1] = rgamma(1, a_1 +p*k/2,
                         1 + 0.5 * sum(tau_m1 * omega))
  for (h in 2:k){
    tau_mh <- (cumprod(delta_post) / delta[h])[-(1:(h-1))]
    delta_post[h] = rgamma(1, a_2 + p/2 * (k-h+1),
                           1 + 0.5 * sum(tau_mh * omega[-(1:(h - 1))]))
  }
  return(delta_post)
}


## ----Initialisation--------------------------------------------------
# Output
lambda_output = array(dim = c(p,k,n_iterations))
sigma_2_output = array(dim = c(p,n_iterations))
eta_output = array(dim = c(n,k,n_iterations))
phi_output = array(dim = c(p,k,n_iterations))
tau_output = array(dim = c(k,n_iterations))

#Initialisation
#On simule tau
delta_courant = c(rgamma(1,a_1,1),rgamma(k-1,a_2,1))

lambda_output[,,1] = rnorm(p*k,0,1) # ne sert à rien en principe
sigma_2_output[,1] = 1/(rgamma(p,a_sigma,b_sigma)) # tirage selon le prior
eta_output[,,1] = rnorm(n*k,0,1) # tirage selon le prior
phi_output[,,1] = matrix(rgamma(p*k,nu/2,nu/2),p,k)
tau_output[,1] = cumprod(delta_courant)


## --------------------------------------------------------------------
for (i in 2:n_iterations){
  lambda_output[,,i] = step1(Y,eta_output[,,i-1],
                             sigma_2_output[,i-1],
                             phi_output[,,i-1],tau_output[,i-1])
  sigma_2_output[,i] = step2(Y,eta_output[,,i-1],
                             lambda_output[,,i])
  eta_output[,,i] = step3(Y,lambda_output[,,i],
                          sigma_2_output[,i])
  phi_output[,,i] = step4(lambda_output[,,i],tau_output[,i-1])
  delta_courant = step5(lambda_output[,,i],phi_output[,,i],
                        delta_courant)
  
  tau_output[,i] = cumprod(delta_courant)
}

##-------SVD--------------
#On s'assure qu'aprés svd on retrouve les mémes que ceux simuler
A=svd(eta) # svd(eta)= UDV^T
D <- diag(A$d) # matrice diagonale de eta*t(eta)
B=A$u %*% D %*% t(lambda %*% A$v) # B=U*D*t(lambda *V)
Y_svd=B+t(epsilon) 


## ----plot lambda-----------------------------------------------------
# On calcule la moyenne a posteriori (estimateur)
apply(lambda_output, 
      MARGIN = c(1, 2), 
      mean) 
print(lambda)
lambda_df <- format_array(array_ = lambda_output, 
                          param_name = "Lambda") 
lambda_df  
# On affiche la chaine
ggplot(lambda_df) + 
  aes(x = iteration, y = Estimate) + 
  facet_wrap(~Parameter, 
             labeller = label_parsed,  
             nrow = 10) + 
  geom_line() 

# On affiche les densités
ggplot(lambda_df) + 
  aes(x = Estimate) + 
  facet_wrap(~Parameter, 
             labeller = label_parsed, 
             scales = "free_y", 
             nrow = 10) + 
  geom_density() 


## ----plot sigma_2----------------------------------------------------
# On calcule la moyenne a posteriori
mean_sigma_2 <- rowMeans(sigma_2_output)
mean_sigma_2

# Plots

# Premièrement, on formate en data.frame (df)
sigma_2_df <- format_matrix(matrix_ = sigma_2_output, 
                            param_name = "sigma", 
                            suffix_ = "^2") 
view(sigma_2_df)

# On affiche la chaine
ggplot(sigma_2_df) + 
  aes(x = iteration, y = Estimate) + 
  facet_wrap(~Parameter, 
             labeller = label_parsed, 
             nrow = 2) + 
  geom_line() 

# On affiche les densités
ggplot(sigma_2_df) + 
  aes(x = Estimate) + 
  facet_wrap(~Parameter, 
             labeller = label_parsed,  
             scales = "free_y", 
             nrow = 2) + 
  geom_density() 


## ----plot phi--------------------------------------------------------
# On calcule la moyenne a posteriori
apply(phi_output, 
      MARGIN = c(1, 2), 
      mean) 

eta_df <- format_array(array_ = phi_output, 
                       param_name = "phi") 
eta_df  
# On affiche la chaine
ggplot(eta_df) + 
  aes(x = iteration, y = Estimate) + 
  facet_wrap(~Parameter, 
             labeller = label_parsed,  
             nrow = 10) + 
  geom_line() 

# On affiche les densités
ggplot(eta_df) + 
  aes(x = Estimate) + 
  facet_wrap(~Parameter, 
             labeller = label_parsed, 
             scales = "free_y", 
             nrow = 10) + 
  geom_density() 


## ----plot tau--------------------------------------------------------
# On calcule la moyenne a posteriori
mean_tau <- rowMeans(tau_output)
mean_tau

# Plots

# Premièrement, on formate en data.frame (df)
tau_df <- format_matrix(matrix_ = tau_output, 
                        param_name = "tau", 
                        suffix_ = "^2") 
view(tau_df)

# On affiche la chaine
ggplot(tau_df) +  
  aes(x = iteration, y = Estimate) + 
  facet_wrap(~Parameter, 
             labeller = label_parsed, 
             nrow = 2) + 
  geom_line() 

# On affiche les densités
ggplot(tau_df) + 
  aes(x = Estimate) + 
  facet_wrap(~Parameter, 
             labeller = label_parsed,  
             scales = "free_y", 
             nrow = 2) + 
  geom_density() 

## ----plot eta------------------------------------------------
apply(eta_output[1:3,,], c(1,2), mean)
print(eta[1:3,])
