# Fonction principale -----------------------------------------------------

eta_lambda_normalise<-function(eta,lambda)
{
  # Initialement, eta = UDV'
  # Donc eta Lambda' = UDV'Lambda' = UD (Lambda V)'
  # donc eta_svd = UD et Lambda_svd = Lambda V
  resultat_svd = svd(eta) # svd(eta)= UDV^T
  D <- diag(resultat_svd$d) # matrice diagonale de eta*t(eta)
  V <- resultat_svd$v
  # Ici, la version normalisée est obtenue en multipliant à droite par V
  # On capitalise sur le fait que V'V = I
  eta_svd <- eta %*% V # = UDV'V = UD
  lambda_svd <- lambda %*% V
  return(list(eta = eta_svd,
              lambda = lambda_svd))
}


# Test (section à virer ou remplacer) -------------------------------------

source("utils_generation_donnees_simulees.R") # Pour obtenir eta et lambda

exemple <- eta_lambda_normalise(eta, lambda)
eta_norm <- exemple$eta
# On vérifie que les 2 colonnes sont orthogonales
sum(eta_norm[, 1] * eta_norm[, 2]) # 0 à la précision machine près

