################################################################################
#                      TP : INTERVALLES DE CONFIANCE 
################################################################################

rm(list=ls())
set.seed(2026)

alpha <- 0.05
nrep <- 1000

vec_n <- seq(20, 100, by=20)

############# 1. LOI DE BERNOULLI B(p) ########################################

vec_p_bern <- seq(0, 0.5, by=0.1)

f_bernoulli <- function(n, p){
  bornes_IC_1 <- matrix(0, nrep, 2)
  bornes_IC_2 <- matrix(0, nrep, 2)
  
  for(i in 1:nrep){
    X <- rbinom(n, size=1, prob=p)
    p_tilde <- mean(X)
    
    # IC 1 : Conservateur (majoration de la variance par 1/4)
    bornes_IC_1[i,] <- p_tilde + c(-1, 1) * (1 / (2 * sqrt(n))) * qnorm(1 - alpha/2)
    
    # IC 2 : Asymptotique (Slutsky, variance estimable)
    bornes_IC_2[i,] <- p_tilde + c(-1, 1) * sqrt(p_tilde * (1 - p_tilde) / n) * qnorm(1 - alpha/2)
  }
  
  return(c(
    mean((p <= bornes_IC_1[,2]) * (p >= bornes_IC_1[,1])),
    mean((p <= bornes_IC_2[,2]) * (p >= bornes_IC_2[,1]))
  ))
}

tableau_bern_1 <- matrix(0, length(vec_n), length(vec_p_bern))
tableau_bern_2 <- matrix(0, length(vec_n), length(vec_p_bern))

for(i in 1:length(vec_n)){
  for(j in 1:length(vec_p_bern)){
    res <- f_bernoulli(vec_n[i], vec_p_bern[j])
    tableau_bern_1[i,j] <- res[1]
    tableau_bern_2[i,j] <- res[2]
  }
}

tableau_bern_1 <- as.data.frame(tableau_bern_1)
tableau_bern_2 <- as.data.frame(tableau_bern_2)
colnames(tableau_bern_1) <- colnames(tableau_bern_2) <- paste0("p=", vec_p_bern)
rownames(tableau_bern_1) <- rownames(tableau_bern_2) <- paste0("n=", vec_n)

cat("--- BERNOULLI : IC Conservateur ---\n")
print(tableau_bern_1)
cat("\n--- BERNOULLI : IC Slutsky ---\n")
print(tableau_bern_2)


############# 2. LOI DE POISSON P(lambda) ######################################

vec_lambda <- c(0.5, 1, 2, 5, 10)

f_poisson <- function(n, lambda){
  bornes_IC <- matrix(0, nrep, 2)
  
  for(i in 1:nrep){
    X <- rpois(n, lambda=lambda)
    lambda_hat <- mean(X)
    
    bornes_IC[i,] <- lambda_hat + c(-1, 1) * qnorm(1 - alpha/2) * sqrt(lambda_hat / n)
  }
  
  return(mean((lambda <= bornes_IC[,2]) * (lambda >= bornes_IC[,1])))
}

tableau_poisson <- matrix(0, length(vec_n), length(vec_lambda))

for(i in 1:length(vec_n)){
  for(j in 1:length(vec_lambda)){
    tableau_poisson[i,j] <- f_poisson(vec_n[i], vec_lambda[j])
  }
}

tableau_poisson <- as.data.frame(tableau_poisson)
colnames(tableau_poisson) <- paste0("lambda=", vec_lambda)
rownames(tableau_poisson) <- paste0("n=", vec_n)

cat("\n--- POISSON : Niveau de confiance empirique (IC Slutsky) ---\n")
print(tableau_poisson)


############# 3. LOI EXPONENTIELLE E(lambda) ###################################


vec_lambda_exp <- c(0.5, 1, 2, 3)

f_exponentielle <- function(n, lambda){
  bornes_IC_exact_asympt <- matrix(0, nrep, 2)
  bornes_IC_approx <- matrix(0, nrep, 2)
  
  for(i in 1:nrep){
    X <- rexp(n, rate=lambda)
    lambda_hat <- 1 / mean(X)
    q <- qnorm(1 - alpha/2)
    
    # Version exacte issue de la Delta-Méthode
    bornes_IC_exact_asympt[i,] <- c(lambda_hat / (1 + q / sqrt(n)), lambda_hat / (1 - q / sqrt(n)))
    
    # Version approximée symétrique
    bornes_IC_approx[i,] <- lambda_hat + c(-1, 1) * q * lambda_hat / sqrt(n)
  }
  
  return(c(
    mean((lambda <= bornes_IC_exact_asympt[,2]) * (lambda >= bornes_IC_exact_asympt[,1])),
    mean((lambda <= bornes_IC_approx[,2]) * (lambda >= bornes_IC_approx[,1]))
  ))
}

tableau_exp_1 <- matrix(0, length(vec_n), length(vec_lambda_exp))
tableau_exp_2 <- matrix(0, length(vec_n), length(vec_lambda_exp))

for(i in 1:length(vec_n)){
  for(j in 1:length(vec_lambda_exp)){
    res <- f_exponentielle(vec_n[i], vec_lambda_exp[j])
    tableau_exp_1[i,j] <- res[1]
    tableau_exp_2[i,j] <- res[2]
  }
}

tableau_exp_1 <- as.data.frame(tableau_exp_1)
tableau_exp_2 <- as.data.frame(tableau_exp_2)
colnames(tableau_exp_1) <- colnames(tableau_exp_2) <- paste0("lambda=", vec_lambda_exp)
rownames(tableau_exp_1) <- rownames(tableau_exp_2) <- paste0("n=", vec_n)

cat("\n--- EXPONENTIELLE : Delta-Méthode (Forme exacte) ---\n")
print(tableau_exp_1)
cat("\n--- EXPONENTIELLE : Delta-Méthode (Forme approximée) ---\n")
print(tableau_exp_2)


############# 4. LOI UNIFORME U[0, theta] #####################################


vec_theta <- c(1, 2, 5, 10)

f_uniforme_theta <- function(n, theta){
  bornes_IC <- matrix(0, nrep, 2)
  
  # Quantiles de la loi Beta(n, 1)
  a_quant <- qbeta(alpha / 2, shape1=n, shape2=1)
  b_quant <- qbeta(1 - alpha / 2, shape1=n, shape2=1)
  
  for(i in 1:nrep){
    X <- runif(n, min=0, max=theta)
    theta_hat <- max(X)
    
    # Intervalle exact non asymptotique
    bornes_IC[i,] <- c(theta_hat / b_quant, theta_hat / a_quant)
  }
  
  return(mean((theta <= bornes_IC[,2]) * (theta >= bornes_IC[,1])))
}

tableau_unif_theta <- matrix(0, length(vec_n), length(vec_theta))

for(i in 1:length(vec_n)){
  for(j in 1:length(vec_theta)){
    tableau_unif_theta[i,j] <- f_uniforme_theta(vec_n[i], vec_theta[j])
  }
}

tableau_unif_theta <- as.data.frame(tableau_unif_theta)
colnames(tableau_unif_theta) <- paste0("theta=", vec_theta)
rownames(tableau_unif_theta) <- paste0("n=", vec_n)

cat("\n--- UNIFORME U[0, theta] : Niveau de confiance exact (Variable pivotale Beta) ---\n")
print(tableau_unif_theta)


############# 5. LOI UNIFORME U[a, b] (Borne inf a inconnue, b connu) ##########
vec_a <- c(0, 1, 2)
b_fixe <- 10

f_uniforme_ab <- function(n, a_val){
  bornes_IC <- matrix(0, nrep, 2)
  
  for(i in 1:nrep){
    X <- runif(n, min=a_val, max=b_fixe)
    a_hat <- min(X)
    
    # IC bilatéral avec pivot exponentiel
    bornes_IC[i,] <- c(
      a_hat + log(alpha / 2) * (b_fixe - a_hat) / n,
      a_hat + log(1 - alpha / 2) * (b_fixe - a_hat) / n
    )
  }
  
  return(mean((a_val <= bornes_IC[,2]) * (a_val >= bornes_IC[,1])))
}

tableau_unif_ab <- matrix(0, length(vec_n), length(vec_a))

for(i in 1:length(vec_n)){
  for(j in 1:length(vec_a)){
    tableau_unif_ab[i,j] <- f_uniforme_ab(vec_n[i], vec_a[j])
  }
}

tableau_unif_ab <- as.data.frame(tableau_unif_ab)
colnames(tableau_unif_ab) <- paste0("a=", vec_a)
rownames(tableau_unif_ab) <- paste0("n=", vec_n)

cat("\n--- UNIFORME U[a, b] : Estimation de a (Pivot Exp(1)) ---\n")
print(tableau_unif_ab)


############# 6. LOI NORMALE N(mu, sigma^2) : ESTIMATION DE MU #################

vec_mu <- c(-2, 0, 2, 5)
sigma_vrai <- 2

# ------------------------------------------------------------------------------
# 6.1 Cas 1 : Variance sigma^2 CONNUE
# ------------------------------------------------------------------------------

f_normale_mu_connu <- function(n, mu){
  bornes_IC <- matrix(0, nrep, 2)
  q <- qnorm(1 - alpha/2)
  
  for(i in 1:nrep){
    X <- rnorm(n, mean=mu, sd=sigma_vrai)
    mu_hat <- mean(X)
    
    bornes_IC[i,] <- mu_hat + c(-1, 1) * q * sigma_vrai / sqrt(n)
  }
  
  return(mean((mu <= bornes_IC[,2]) * (mu >= bornes_IC[,1])))
}

tableau_norm_mu_connu <- matrix(0, length(vec_n), length(vec_mu))

for(i in 1:length(vec_n)){
  for(j in 1:length(vec_mu)){
    tableau_norm_mu_connu[i,j] <- f_normale_mu_connu(vec_n[i], vec_mu[j])
  }
}

tableau_norm_mu_connu <- as.data.frame(tableau_norm_mu_connu)
colnames(tableau_norm_mu_connu) <- paste0("mu=", vec_mu)
rownames(tableau_norm_mu_connu) <- paste0("n=", vec_n)

cat("\n--- NORMALE N(mu, sigma^2) : Estimation de mu (sigma^2 connu) ---\n")
print(tableau_norm_mu_connu)


# ------------------------------------------------------------------------------
# 6.2 Cas 2 : Variance sigma^2 INCONNUE (Théorème de Cochran & Student)
# ------------------------------------------------------------------------------

f_normale_mu_inconnu <- function(n, mu){
  bornes_IC <- matrix(0, nrep, 2)
  q_student <- qt(1 - alpha/2, df = n - 1)
  
  for(i in 1:nrep){
    X <- rnorm(n, mean=mu, sd=sigma_vrai)
    mu_hat <- mean(X)
    S_n <- sd(X) # Ecart-type empirique corrigé
    
    bornes_IC[i,] <- mu_hat + c(-1, 1) * q_student * S_n / sqrt(n)
  }
  
  return(mean((mu <= bornes_IC[,2]) * (mu >= bornes_IC[,1])))
}

tableau_norm_mu_inconnu <- matrix(0, length(vec_n), length(vec_mu))

for(i in 1:length(vec_n)){
  for(j in 1:length(vec_mu)){
    tableau_norm_mu_inconnu[i,j] <- f_normale_mu_inconnu(vec_n[i], vec_mu[j])
  }
}

tableau_norm_mu_inconnu <- as.data.frame(tableau_norm_mu_inconnu)
colnames(tableau_norm_mu_inconnu) <- paste0("mu=", vec_mu)
rownames(tableau_norm_mu_inconnu) <- paste0("n=", vec_n)

cat("\n--- NORMALE N(mu, sigma^2) : Estimation de mu (sigma^2 inconnu, Student T(n-1)) ---\n")
print(tableau_norm_mu_inconnu)


################################################################################
############# 7. LOI NORMALE N(mu, sigma^2) : ESTIMATION DE SIGMA^2 ############

vec_sigma2 <- c(0.5, 1, 4, 9)
mu_fixe <- 0

f_normale_sigma2 <- function(n, sigma2_val){
  bornes_IC <- matrix(0, nrep, 2)
  q_low <- qchisq(alpha / 2, df = n - 1)
  q_high <- qchisq(1 - alpha / 2, df = n - 1)
  
  for(i in 1:nrep){
    X <- rnorm(n, mean=mu_fixe, sd=sqrt(sigma2_val))
    S_n2 <- var(X) # Variance empirique corrigée
    
    bornes_IC[i,] <- c( (n - 1) * S_n2 / q_high, (n - 1) * S_n2 / q_low )
  }
  
  return(mean((sigma2_val <= bornes_IC[,2]) * (sigma2_val >= bornes_IC[,1])))
}

tableau_norm_sigma2 <- matrix(0, length(vec_n), length(vec_sigma2))

for(i in 1:length(vec_n)){
  for(j in 1:length(vec_sigma2)){
    tableau_norm_sigma2[i,j] <- f_normale_sigma2(vec_n[i], vec_sigma2[j])
  }
}

tableau_norm_sigma2 <- as.data.frame(tableau_norm_sigma2)
colnames(tableau_norm_sigma2) <- paste0("sigma2=", vec_sigma2)
rownames(tableau_norm_sigma2) <- paste0("n=", vec_n)

cat("\n--- NORMALE N(mu, sigma^2) : Estimation de sigma^2 (Loi Chi2(n-1)) ---\n")
print(tableau_norm_sigma2)