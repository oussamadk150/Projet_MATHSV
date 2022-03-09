

## --------------------------------------------------------------------
rm(list = ls()) #Nettoyage de l'environnement

## ----librariries-----------------------------------------------------
library(tidyverse) # Pour la manipulation de données
library(MASS)      # Pour mvrnorm
library(mixtools) # Pour dmvnorm (densite d'une loi normale multivariee)
source("utils_formatting_functions.R")
source("utils_generation_donnees_simulees.R")
source("utils_gibbs_functions.R")

## ----fonction qui calcul le BIC du modèle en fonction de k------------------------------------------------
evaluation_k_BIC <- function(Y, K){
  n <- nrow(Y)
  p <- ncol(Y)
  BIC = matrix(NA, nrow = K, ncol = 2)
  
  for (k in 2:(K)){
    lois = gibbs_sampler(Y, X, 2, n_iterations=1000,
                     a_1=2, a_2=3, sigma_2_beta=1, 
                     nu=3, a_sigma=1, b_sigma=0.3)
    lambda_output = lois[1]
    sigma_2_output = lois[2]
    omega = apply(lambda_output, 
                  MARGIN = c(1, 2), 
                  mean)%*%t(apply(lambda_output, 
                                  MARGIN = c(1, 2), 
                                  mean)) + diag(rowMeans(sigma_2_output))
    #Warning
    #On ne prend pas encore compte beta dans le calcul de la vraissemblance
    L = sum(apply(Y,1,function(y){
      mixtools::logdmvnorm(y,rep(0,p),omega)
    }))
    BIC[k,] = c(k,L - ((k*p + p)/2)*log(n))
  }
  colnames(BIC) <- c("iterations","BIC")
  as.data.frame(BIC)
}
valeurs_BIC <- evaluation_k_BIC(Y, 10)

## ----On trace la valeur du BIC en fonction de k-----------------------------------------------------
ggplot(valeurs_BIC, aes(x = iterations)) + 
  geom_point(aes(y = BIC, colour = "BIC")) +
  labs(x = "Valeur de k", y = "Valeur BIC", title = "Valeur du BIC en fonction de k",
       colour = "Mon nom")
