

## --------------------------------------------------------------------
rm(list = ls()) #Nettoyage de l'environnement

## ----librariries-----------------------------------------------------
library(tidyverse) # Pour la manipulation de données
library(MASS)      # Pour mvrnorm
library(mixtools) # Pour dmvnorm (densite d'une loi normale multivariee)
source("utils_formatting_functions.R")
source("utils_generation_donnees_simulees.R")
source("utils_gibbs_functions.R")

evaluation_k_BIC <- function(Y, K){
  n <- nrow(Y)
  p <- ncol(Y)
  BIC = matrix(NA, nrow = K, ncol = 2)

  for (k in 2:(K)){
    all_outputs = gibbs_sampler(Y, X, k, n_iterations=100,
                                a_1=2, a_2 = 3, sigma_2_beta = 1, 
                                nu = 3, a_sigma=1, b_sigma = 0.3)
    
    lambda_output <- all_outputs$lambda
    sigma_2_output <- all_outputs$sigma2
    
    omega = apply(lambda_output, 
                  MARGIN = c(1, 2), 
                  mean)%*%t(apply(lambda_output, 
                                  MARGIN = c(1, 2), 
                                  mean)) + diag(rowMeans(sigma_2_output))
    
    L = sum(apply(Y,1,function(y){
      mixtools::logdmvnorm(y,rep(0,p),omega)
    }))
    BIC[k,] = c(k,L - ((k*p + p)/2)*log(n))
  }
  colnames(BIC) <- c("iterations","BIC")
  as.data.frame(BIC)
}
valeurs_BIC <- evaluation_k_BIC(Y, 10)


ggplot(valeurs_BIC, aes(x = iterations)) + 
  geom_point(aes(y = BIC, colour = "BIC")) +
  labs(x = "Valeur de k", y = "Valeur BIC", title = "Valeur du BIC en fonction de k",
       colour = "Mon nom")
