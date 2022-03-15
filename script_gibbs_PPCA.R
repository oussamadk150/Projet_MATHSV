

## --------------------------------------------------------------------
rm(list = ls()) #Nettoyage de l'environnement

## ----librariries-----------------------------------------------------
library(tidyverse) # Pour la manipulation de données
library(parallel)
source("utils_generation_donnees_simulees.R") # Generation de x et y
source("utils_gibbs_functions.R") # Pour les fonctions de gibbs
source("utils_formatting_functions.R") # Pour le formattage des résultas

# Test on differents vaslues of k
tested_ks <- 2:6
all_outputs = mclapply(tested_ks,
                       function(k_){
                         gibbs_sampler(Y, X, k = k_, n_iterations=1000,
                                       a_1 = 2, a_2 = 3, sigma_2_beta = 1, 
                                       nu = 3, a_sigma=1, b_sigma = 0.3)},
                       mc.cores = detectCores())

# Si mclapply ne marche pas sur windows, faites 
# Ce sera plus long
# all_outputs = lapply(tested_ks,
#                        function(k_){
#                          gibbs_sampler(Y, X, k = k_, n_iterations=1000,
#                                        a_1 = 2, a_2 = 3, sigma_2_beta = 1, 
#                                        nu = 3, a_sigma=1, b_sigma = 0.3)})

# BIC ---------------------------------------------------------------------

BICs <- tibble(K = tested_ks,
               BIC = map_dbl(all_outputs, function(x) x$BIC))
ggplot(BICs) +
  aes(x = K, y = BIC) + 
  geom_point()

# posterior distributions of taus

# We focus in the case k = 6

(all_outputs[[which(tested_ks == 6)]]$tau) %>% 
  format_matrix("tau") %>% 
  ggplot(aes(x = Parameter, y = Estimate)) +
  geom_boxplot()

# Check results on a specified model

mon_k <- 6

output <- all_outputs[[which(tested_ks == mon_k)]]

# Estimation of beta
beta_df <- format_array(output$beta, "beta")
ggplot(beta_df) +
  aes(x = Estimate) +
  geom_density(aes(color = "Estimation")) +
  facet_wrap(~Parameter, nrow = nrow(beta),
             labeller = label_parsed) +
  geom_vline(data = data.frame(truth = as.numeric(t(beta)),
                               Parameter = levels(beta_df$Parameter)),
             aes(xintercept = truth, color = "Truth")) +
  labs(y = "Estimated density", colour = "") + 
  scale_color_manual(values = c("black", "red"))

# Estimation of sigma2

# On checke la convergence sur les chaines

sigma2_df <- format_matrix(output$sigma2, "sigma", "^2")
ggplot(sigma2_df) +
  aes(x = iteration, y = Estimate) +
  geom_line() +
  facet_wrap(~Parameter, nrow = 2, labeller = label_parsed) + 
  scale_y_continuous(trans = "log10")

# Estimation of Omega

Lambda_posterior <- output$lambda
omega_post <- apply(Lambda_posterior, 
                    MARGIN = c(1, 2), 
                    mean)%*%t(apply(Lambda_posterior, 
                                    MARGIN = c(1, 2), 
                                    mean)) + 
  diag(rowMeans(output$sigma2))

omega_true <- round(lambda%*% t(lambda) + sigma, 2)

par(mfrow = c(1, 2))
corrplot::corrplot(cov2cor(omega_post))
corrplot::corrplot(cov2cor(omega_true))
par(mfrow = c(1, 1))
