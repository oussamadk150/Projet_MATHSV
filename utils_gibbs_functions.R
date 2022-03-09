## ----librariries-----------------------------------------------------
library(tidyverse) # Pour la manipulation de donn?es
library(MASS)      # Pour mvrnorm


## ----step_1----------------------------------------------------------
# mise ? jour de lambda
step1 <- function(Y,eta,sigma_2,phi,tau,k){
  p <- ncol(Y)
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
# mise ? jour de sigma_2
step2 <- function(Y,eta,lambda,a_sigma,b_sigma){
  p <- ncol(Y)
  n <- nrow(Y)
  sigma_2_post = rep(0,p)
  a_post = a_sigma + n/2
  for (j in 1:p) {
    b_post = b_sigma + 0.5*sum((Y[,j]-t(lambda[j,])%*%t(eta))**2)
    sigma_2_post[j] = 1/(rgamma(1,a_post,b_post))
  }
  return(sigma_2_post)
}


## ----step_3----------------------------------------------------------
# mise ? jour de eta
step3 <- function(Y,lambda,sigma_2,k){
  n <- nrow(Y)
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
# mise ? jour de eta
step4 <- function(lambda,tau,k,nu){
  p <- ncol(Y)
  phi_post = matrix(0,p,k)
  for (j in 1:p) {
    phi_post[j,] = rgamma(1,(nu+1)/2, (nu+tau*(lambda[j,])**2)/2)
  }
  return(phi_post)
}


## ----step 5----------------------------------------------------------
# mise ? jour de delta
step5 <- function(lambda,phi, delta,k,a_1,a_2){
  p <- ncol(Y)
  delta_post = delta # Initialisation ? NA
  phi_lambda2 <- phi * lambda^2
  # cumprod donne le produit cumul?
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
## ----step 6----------------------------------------------------------
# mise ? jour de beta
step6 <- function(Y,X,lambda,eta,sigma_2,sigma_2_beta){
  p <- ncol(Y)
  q <- ncol(X)
  beta_post = matrix(0, nrow = p, ncol = q)
  precision_prior_beta <- diag(1 / sigma_2_beta, nrow = q)
  for (j in 1:p) {
    var_post = solve(precision_prior_beta + t(X) %*% X / sigma_2[j])
    moy_post = var_post %*% t(X) %*% (Y[,j] - eta %*% lambda[j,]) / sigma_2[j]
    beta_post[j,] = mvrnorm(1, moy_post, var_post)
  }
  return(t(beta_post))
}

gibbs_sampler <- function(Y, X, k, n_iterations,
                          a_1, a_2, sigma_2_beta, 
                          nu, a_sigma, b_sigma){
  ## ----Initialisation--------------------------------------------------
  # Output
  p <- ncol(Y)
  n <- nrow(Y)
  q <- ncol(X)
  lambda_output = array(dim = c(p, k, n_iterations))
  sigma_2_output = array(dim = c(p,n_iterations))
  eta_output = array(dim = c(n,k,n_iterations))
  phi_output = array(dim = c(p,k,n_iterations))
  tau_output = array(dim = c(k,n_iterations))
  beta_output = array(dim = c(q,p,n_iterations))
  
  #Initialisation
  #On simule tau
  delta_courant = c(rgamma(1,a_1,1),rgamma(k-1,a_2,1))
  
  lambda_output[,,1] = rnorm(p*k,0,1) # ne sert ? rien en principe
  sigma_2_output[,1] = 1/(rgamma(p,a_sigma,b_sigma)) # tirage selon le prior
  eta_output[,,1] = rnorm(n*k,0,1) # tirage selon le prior
  phi_output[,,1] = matrix(rgamma(p*k,nu/2,nu/2),p,k)
  tau_output[,1] = cumprod(delta_courant)
  beta_output[,,1] = rnorm(q*p,0,sigma_2_beta)
  
  ## --------------------------------------------------------------------
  for (i in 2:n_iterations){
    lambda_output[,,i] = step1(Y-X%*%beta_output[,,i-1],
                               eta_output[,,i-1],
                               sigma_2_output[,i-1],
                               phi_output[,,i-1],tau_output[,i-1],k)
    sigma_2_output[,i] = step2(Y-X%*%beta_output[,,i-1],
                               eta_output[,,i-1],
                               lambda_output[,,i],a_sigma,b_sigma)
    eta_output[,,i] = step3(Y-X%*%beta_output[,,i-1],
                            lambda_output[,,i],
                            sigma_2_output[,i],k)
    phi_output[,,i] = step4(lambda_output[,,i],tau_output[,i-1],k,nu)
    delta_courant = step5(lambda_output[,,i],phi_output[,,i],
                          delta_courant,k,a_1,a_2)
    
    tau_output[,i] = cumprod(delta_courant)
    beta_output[,,i] = step6(Y,X, lambda_output[,,i],
                             eta_output[,,i], sigma_2_output[,i],sigma_2_beta)
  }
  # Return a named list
  return(list(lambda = lambda_output,
              sigma2 = sigma_2_output,
              eta = eta_output,
              phi = phi_output,
              tau = tau_output,
              beta = beta_output))
}

# Test (section à supprimer ou déplacer) ----------------------------------

source("utils_generation_donnees_simulees.R") # Generation de x et y

all_outputs = gibbs_sampler(Y, X, k = 2, n_iterations=1000,
                            a_1=2, a_2 = 3, sigma_2_beta = 1, 
                            nu = 3, a_sigma=1, b_sigma = 0.3)


# Pour extraire un élément, on utilise ma_liste$nom_element
lambda_output <- all_outputs$lambda

# Calcul de moyenne a posteriori

apply(lambda_output, 
      MARGIN = c(1, 2), 
      mean) 

# Formattage pour ggplot

source("utils_formatting_functions.R") # Chargement de fonctions de formattage
lambda_df <- format_array(lambda_output, "Lambda")
