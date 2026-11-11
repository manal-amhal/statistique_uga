###############################################################################
#               TP1
###############################################################################

rm(list=ls())

############## Exercice 1 #####################################################

n=1000

theta=5

set.seed(2026)

X=rpois(n,theta)

theta_MV=mean(X) #max de vraisemblance
theta_mom_1=mean(X) #moment ordre 1
theta_mom_2=var(X) #moment ordre 2

c(theta_MV,theta_mom_1,theta_mom_2) #bilan des 3

M=1000

vec_theta_MV=rep(0,M)
vec_theta_mom_1=rep(0,M)
vec_theta_mom_2=rep(0,M)

i=1
while(i<=M){
  
  X=rpois(n,theta)
  
  vec_theta_MV[i]=mean(X)
  vec_theta_mom_1[i]=mean(X)
  vec_theta_mom_2[i]=var(X)
  
  i=i+1
}

#comparaison des 3
quartz()
boxplot(cbind(vec_theta_MV,vec_theta_mom_1,vec_theta_mom_2),
        names=c(expression(hat(theta)[MV]),expression(hat(theta)[1]),expression(hat(theta)[2])))
abline(a=theta,b=0,col="red",lwd=2)


epsilon=0.1
c(mean(abs(vec_theta_MV-theta)>epsilon),
  mean(abs(vec_theta_mom_1-theta)>epsilon),
  mean(abs(vec_theta_mom_2-theta)>epsilon))#proportion de valeurs éloignées

hist(sqrt(n)*(vec_theta_MV-theta),breaks=25,freq=F,
     xlab=expression(hat(theta)[MV]),main="Histogramme") #Histogramme
points(y=dnorm(seq(-7,7,length=500),sd=sqrt(theta)),x=seq(-7,7,length=500),
       type="l",col="red",lwd=2) #loi limite


########### Exercice 2 #####################################################

n=100

theta=1

set.seed(2026)

X=runif(n,0,theta)

theta_MV=max(X) #max de vraisemblance
theta_mom_1=2*mean(X) #moment ordre 1

thetaxe=seq(0,2,length=200)
f=function(theta){ #vraisemblance
  (1/theta^n)*prod((X>=0)*(X<=theta))
}

quartz()
plot(y=mapply(f,thetaxe),x=thetaxe,type="l",
     lwd=2,main="Vraisemblance",xlab=expression(theta),
     ylab="Vraisemblance")

abline(v=max(X),lty=5,lwd=2,col="red")

i=1
while(i<=9){
  
  X=runif(n,0,theta)
  
  points(y=mapply(f,thetaxe),x=thetaxe,type="l",lwd=2)
  
  i=i+1
}

c(theta_MV,theta_mom_1) #bilan des 2

M=1000

vec_theta_MV=rep(0,M)
vec_theta_mom_1=rep(0,M)

i=1
while(i<=M){
  
  X=runif(n,0,theta)
  
  vec_theta_MV[i]=max(X)
  vec_theta_mom_1[i]=2*mean(X)
  
  i=i+1
}

#comparaison des 2
quartz()
boxplot(cbind(vec_theta_MV,vec_theta_mom_1),
        names=c(expression(hat(theta)[MV]),expression(hat(theta)[1])))
abline(a=theta,b=0,col="red",lwd=2)

epsilon=0.1
c(mean(abs(vec_theta_MV-theta)>epsilon),
  mean(abs(vec_theta_mom_1-theta)>epsilon))#proportion de valeurs éloignées

quartz()
hist(sqrt(n)*(vec_theta_mom_1-theta),breaks=25,freq=F,
     xlab=expression(hat(theta)[MV]),main="Histogramme") #Histogramme
points(y=dnorm(seq(-2.5,2.5,length=500),sd=sqrt(theta^2/3)),x=seq(-2.5,2.5,length=500),
       type="l",col="red",lwd=2) #loi limite


###############################################################################
#           Initialisation commune pour la suite
###############################################################################

rm(list = ls())
set.seed(2026)

M <- 1000        # nombre de répétitions
epsilon <- 0.1   # seuil pour la proportion d'estimations éloignées


###############################################################################
#           TP Appliquer à d'autres lois
###############################################################################


############## Loi de Poisson  ################################################


n <- 1000
theta <- 5

vec_theta_MV  <- numeric(M)
vec_theta_MM1 <- numeric(M)
vec_theta_MM2 <- numeric(M)

for(i in 1:M) {
  X <- rpois(n, theta)
  vec_theta_MV[i]  <- mean(X)  # EMV = mean(X)
  vec_theta_MM1[i] <- mean(X)  # MM1 = E[X] = mean(X)
  vec_theta_MM2[i] <- var(X)   # MM2 = Var(X) = var(X)
}

# Boxplot comparative
quartz()
boxplot(cbind(vec_theta_MV, vec_theta_MM1, vec_theta_MM2),
        names = c(expression(hat(theta)[MV]), expression(hat(theta)[MM1]), expression(hat(theta)[MM2])),
        main = "Poisson P(5)")
abline(h = theta, col = "red", lwd = 2)

# Déviations epsilon
c(Poisson_MV  = mean(abs(vec_theta_MV - theta) > epsilon),
  Poisson_MM1 = mean(abs(vec_theta_MM1 - theta) > epsilon),
  Poisson_MM2 = mean(abs(vec_theta_MM2 - theta) > epsilon))

# Histogramme asymptotique TCL pour EMV
quartz()
hist(sqrt(n) * (vec_theta_MV - theta), breaks = 25, freq = FALSE,
     xlab = expression(sqrt(n)*(hat(theta)[MV] - theta)), main = "Loi limite - Poisson")
curve(dnorm(x, mean = 0, sd = sqrt(theta)), add = TRUE, col = "red", lwd = 2)


############## Loi de Bernoulli  ##############################################
n <- 1000
p <- 0.4

vec_p_MV <- numeric(M)
vec_p_MM <- numeric(M)

for(i in 1:M) {
  X <- rbinom(n, size = 1, prob = p)
  vec_p_MV[i] <- mean(X)  # EMV = mean(X)
  vec_p_MM[i] <- mean(X)  # MM = E[X] = mean(X)
}

quartz()
boxplot(cbind(vec_p_MV, vec_p_MM),
        names = c(expression(hat(p)[MV]), expression(hat(p)[MM])),
        main = "Bernoulli B(0.4)")
abline(h = p, col = "red", lwd = 2)

c(Bernoulli_MV = mean(abs(vec_p_MV - p) > epsilon),
  Bernoulli_MM = mean(abs(vec_p_MM - p) > epsilon))

quartz()
hist(sqrt(n) * (vec_p_MV - p), breaks = 25, freq = FALSE,
     xlab = expression(sqrt(n)*(hat(p)[MV] - p)), main = "Loi limite - Bernoulli")
curve(dnorm(x, mean = 0, sd = sqrt(p * (1 - p))), add = TRUE, col = "red", lwd = 2)


############## Loi Binomiale (k connue)  ######################################
n <- 1000
K <- 20
p <- 0.3

vec_p_MV <- numeric(M)
vec_p_MM <- numeric(M)

for(i in 1:M) {
  X <- rbinom(n, size = K, prob = p)
  vec_p_MV[i] <- mean(X) / K  # EMV = mean(X) / K
  vec_p_MM[i] <- mean(X) / K  # MM = E[X] / K
}

quartz()
boxplot(cbind(vec_p_MV, vec_p_MM),
        names = c(expression(hat(p)[MV]), expression(hat(p)[MM])),
        main = "Binomiale B(20, 0.3)")
abline(h = p, col = "red", lwd = 2)

c(Binomiale_MV = mean(abs(vec_p_MV - p) > epsilon),
  Binomiale_MM = mean(abs(vec_p_MM - p) > epsilon))

quartz()
hist(sqrt(n) * (vec_p_MV - p), breaks = 25, freq = FALSE,
     xlab = expression(sqrt(n)*(hat(p)[MV] - p)), main = "Loi limite - Binomiale")
curve(dnorm(x, mean = 0, sd = sqrt(p * (1 - p) / K)), add = TRUE, col = "red", lwd = 2)


############## Loi Géométrique ################################################
n <- 1000
p <- 0.25

vec_p_MV <- numeric(M)
vec_p_MM <- numeric(M)

for(i in 1:M) {
  X <- rgeom(n, prob = p) + 1  # rgeom dans R compte le nombre d'echecs, d'ou le +1
  vec_p_MV[i] <- 1 / mean(X)  # EMV = 1 / mean(X)
  vec_p_MM[i] <- 1 / mean(X)  # MM = 1 / E[X]
}

quartz()
boxplot(cbind(vec_p_MV, vec_p_MM),
        names = c(expression(hat(p)[MV]), expression(hat(p)[MM])),
        main = "Géométrique Geog(0.25)")
abline(h = p, col = "red", lwd = 2)

c(Geometrique_MV = mean(abs(vec_p_MV - p) > epsilon),
  Geometrique_MM = mean(abs(vec_p_MM - p) > epsilon))

# Delta-méthode pour la variance limite : p^2 * (1-p)
quartz()
hist(sqrt(n) * (vec_p_MV - p), breaks = 25, freq = FALSE,
     xlab = expression(sqrt(n)*(hat(p)[MV] - p)), main = "Loi limite - Géométrique")
curve(dnorm(x, mean = 0, sd = sqrt((p^2) * (1 - p))), add = TRUE, col = "red", lwd = 2)


############## Loi Uniforme discrète ##########################################
n <- 100
N <- 50

vec_N_MV <- numeric(M)
vec_N_MM <- numeric(M)

for(i in 1:M) {
  X <- sample(1:N, size = n, replace = TRUE)
  vec_N_MV[i] <- max(X)                # EMV = max(X)
  vec_N_MM[i] <- 2 * mean(X) - 1       # MM1 : E[X] = (N+1)/2 => N = 2*mean(X) - 1
}

quartz()
boxplot(cbind(vec_N_MV, vec_N_MM),
        names = c(expression(hat(N)[MV]), expression(hat(N)[MM])),
        main = "Uniforme Discrète U({1,...,50})")
abline(h = N, col = "red", lwd = 2)

c(UnifDiscr_MV = mean(abs(vec_N_MV - N) > epsilon),
  UnifDiscr_MM = mean(abs(vec_N_MM - N) > epsilon))

# Histogramme pour l'estimateur des moments (TCL)
quartz()
hist(sqrt(n) * (vec_N_MM - N), breaks = 25, freq = FALSE,
     xlab = expression(sqrt(n)*(hat(N)[MM] - N)), main = "Loi limite MM - Uniforme Discrète")
curve(dnorm(x, mean = 0, sd = sqrt((N^2 - 1) / 3)), add = TRUE, col = "red", lwd = 2)


############## Loi Uniforme continue ##########################################
n <- 100
theta <- 1

vec_theta_MV <- numeric(M)
vec_theta_MM <- numeric(M)

for(i in 1:M) {
  X <- runif(n, min = 0, max = theta)
  vec_theta_MV[i] <- max(X)         # EMV = max(X)
  vec_theta_MM[i] <- 2 * mean(X)    # MM1 = 2 * mean(X)
}

quartz()
boxplot(cbind(vec_theta_MV, vec_theta_MM),
        names = c(expression(hat(theta)[MV]), expression(hat(theta)[MM])),
        main = "Uniforme Continue U([0, 1])")
abline(h = theta, col = "red", lwd = 2)

c(UnifCont_MV = mean(abs(vec_theta_MV - theta) > epsilon),
  UnifCont_MM = mean(abs(vec_theta_MM - theta) > epsilon))

# Histogramme TCL pour l'estimateur des moments (vitesse sqrt(n))
quartz()
hist(sqrt(n) * (vec_theta_MM - theta), breaks = 25, freq = FALSE,
     xlab = expression(sqrt(n)*(hat(theta)[MM] - theta)), main = "Loi limite MM - Uniforme Continue")
curve(dnorm(x, mean = 0, sd = sqrt(theta^2 / 3)), add = TRUE, col = "red", lwd = 2)


############## Loi Exponentielle  #############################################
n <- 1000
lambda <- 2

vec_lambda_MV <- numeric(M)
vec_lambda_MM <- numeric(M)

for(i in 1:M) {
  X <- rexp(n, rate = lambda)
  vec_lambda_MV[i] <- 1 / mean(X)   # EMV = 1 / mean(X)
  vec_lambda_MM[i] <- 1 / mean(X)   # MM = 1 / E[X]
}

quartz()
boxplot(cbind(vec_lambda_MV, vec_lambda_MM),
        names = c(expression(hat(lambda)[MV]), expression(hat(lambda)[MM])),
        main = "Exponentielle E(2)")
abline(h = lambda, col = "red", lwd = 2)

c(Exponentielle_MV = mean(abs(vec_lambda_MV - lambda) > epsilon),
  Exponentielle_MM = mean(abs(vec_lambda_MM - lambda) > epsilon))

# Delta-méthode pour la loi limite de EMV (variance = lambda^2)
quartz()
hist(sqrt(n) * (vec_lambda_MV - lambda), breaks = 25, freq = FALSE,
     xlab = expression(sqrt(n)*(hat(lambda)[MV] - lambda)), main = "Loi limite - Exponentielle")
curve(dnorm(x, mean = 0, sd = lambda), add = TRUE, col = "red", lwd = 2)


############## Loi Normale ################################################
n <- 1000
mu <- 3
sigma2 <- 4
sigma <- sqrt(sigma2)

vec_mu_MV     <- numeric(M)
vec_sigma2_MV <- numeric(M)
vec_sigma2_MM <- numeric(M)

for(i in 1:M) {
  X <- rnorm(n, mean = mu, sd = sigma)
  vec_mu_MV[i]     <- mean(X)                      # EMV mu = mean(X)
  vec_sigma2_MV[i] <- mean((X - mean(X))^2)        # EMV sigma^2 (biaisé)
  vec_sigma2_MM[i] <- var(X)                       # MM / sans biais (n-1)
}

# Boxplot pour mu
quartz()
boxplot(vec_mu_MV, names = expression(hat(mu)[MV]), main = "Normale: mu = 3")
abline(h = mu, col = "red", lwd = 2)

# Boxplot pour sigma^2
quartz()
boxplot(cbind(vec_sigma2_MV, vec_sigma2_MM),
        names = c(expression(hat(sigma)[MV]^2), expression(hat(sigma)[MM]^2)),
        main = "Normale: sigma^2 = 4")
abline(h = sigma2, col = "red", lwd = 2)

c(Normale_mu_MV   = mean(abs(vec_mu_MV - mu) > epsilon),
  Normale_sig2_MV = mean(abs(vec_sigma2_MV - sigma2) > epsilon),
  Normale_sig2_MM = mean(abs(vec_sigma2_MM - sigma2) > epsilon))

# Loi limite de mu (TCL direct)
quartz()
hist(sqrt(n) * (vec_mu_MV - mu), breaks = 25, freq = FALSE,
     xlab = expression(sqrt(n)*(hat(mu)[MV] - mu)), main = "Loi limite - Normale (mu)")
curve(dnorm(x, mean = 0, sd = sigma), add = TRUE, col = "red", lwd = 2)


###############################################################################
#           TP1 Compléments : lois présentes dans les notes manuscrites
#           (même démarche : EMV vs moments, boxplot, déviations, loi limite)
###############################################################################


############## Loi Binomiale négative (r connu)  ##############################
n <- 1000
r <- 4
p <- 0.4

vec_p_MV <- numeric(M)
vec_p_MM <- numeric(M)

for(i in 1:M) {
  X <- rnbinom(n, size = r, prob = p)   # R compte le nombre d'echecs avant le r-ieme succes
  vec_p_MV[i] <- r / (r + mean(X))      # EMV = r / (r + mean(X))
  vec_p_MM[i] <- r / (r + mean(X))      # MM : E[X] = r(1-p)/p => p = r / (E[X] + r)
}

quartz()
boxplot(cbind(vec_p_MV, vec_p_MM),
        names = c(expression(hat(p)[MV]), expression(hat(p)[MM])),
        main = "Binomiale négative BN(4, 0.4)")
abline(h = p, col = "red", lwd = 2)

c(BinNeg_MV = mean(abs(vec_p_MV - p) > epsilon),
  BinNeg_MM = mean(abs(vec_p_MM - p) > epsilon))

# Delta-méthode pour la variance limite : p^2 * (1-p) / r
quartz()
hist(sqrt(n) * (vec_p_MV - p), breaks = 25, freq = FALSE,
     xlab = expression(sqrt(n)*(hat(p)[MV] - p)), main = "Loi limite - Binomiale négative")
curve(dnorm(x, mean = 0, sd = sqrt((p^2) * (1 - p) / r)), add = TRUE, col = "red", lwd = 2)


############## Loi Uniforme continue U([a, 0]) avec a < 0 #####################
n <- 100
a <- -5

vec_a_MV <- numeric(M)
vec_a_MM <- numeric(M)

for(i in 1:M) {
  X <- runif(n, min = a, max = 0)
  vec_a_MV[i] <- min(X)          # EMV = min(X)
  vec_a_MM[i] <- 2 * mean(X)     # MM1 : E[X] = a/2 => a = 2 * mean(X)
}

quartz()
boxplot(cbind(vec_a_MV, vec_a_MM),
        names = c(expression(hat(a)[MV]), expression(hat(a)[MM])),
        main = "Uniforme Continue U([-5, 0])")
abline(h = a, col = "red", lwd = 2)

c(UnifCont_a0_MV = mean(abs(vec_a_MV - a) > epsilon),
  UnifCont_a0_MM = mean(abs(vec_a_MM - a) > epsilon))

# Histogramme TCL pour l'estimateur des moments (vitesse sqrt(n))
quartz()
hist(sqrt(n) * (vec_a_MM - a), breaks = 25, freq = FALSE,
     xlab = expression(sqrt(n)*(hat(a)[MM] - a)), main = "Loi limite MM - Uniforme U([a, 0])")
curve(dnorm(x, mean = 0, sd = abs(a) / sqrt(3)), add = TRUE, col = "red", lwd = 2)


############## Loi Uniforme continue U([a, b]) ################################
n <- 100
a <- 2
b <- 7

vec_a_MV <- numeric(M)
vec_b_MV <- numeric(M)
vec_a_MM <- numeric(M)
vec_b_MM <- numeric(M)

for(i in 1:M) {
  X <- runif(n, min = a, max = b)
  vec_a_MV[i] <- min(X)                          # EMV a = min(X)
  vec_b_MV[i] <- max(X)                          # EMV b = max(X)
  vec_a_MM[i] <- mean(X) - sqrt(3 * var(X))      # MM : a = E[X] - sqrt(3 V(X))
  vec_b_MM[i] <- mean(X) + sqrt(3 * var(X))      # MM : b = E[X] + sqrt(3 V(X))
}

# Boxplot pour a
quartz()
boxplot(cbind(vec_a_MV, vec_a_MM),
        names = c(expression(hat(a)[MV]), expression(hat(a)[MM])),
        main = "Uniforme Continue U([2, 7]) : estimation de a")
abline(h = a, col = "red", lwd = 2)

# Boxplot pour b
quartz()
boxplot(cbind(vec_b_MV, vec_b_MM),
        names = c(expression(hat(b)[MV]), expression(hat(b)[MM])),
        main = "Uniforme Continue U([2, 7]) : estimation de b")
abline(h = b, col = "red", lwd = 2)

c(UnifAB_a_MV = mean(abs(vec_a_MV - a) > epsilon),
  UnifAB_a_MM = mean(abs(vec_a_MM - a) > epsilon),
  UnifAB_b_MV = mean(abs(vec_b_MV - b) > epsilon),
  UnifAB_b_MM = mean(abs(vec_b_MM - b) > epsilon))

# Histogrammes TCL pour les estimateurs des moments (variance limite : 2 (b-a)^2 / 15)
quartz()
par(mfrow = c(1, 2))
hist(sqrt(n) * (vec_a_MM - a), breaks = 25, freq = FALSE,
     xlab = expression(sqrt(n)*(hat(a)[MM] - a)), main = "Loi limite MM - a")
curve(dnorm(x, mean = 0, sd = (b - a) * sqrt(2 / 15)), add = TRUE, col = "red", lwd = 2)
hist(sqrt(n) * (vec_b_MM - b), breaks = 25, freq = FALSE,
     xlab = expression(sqrt(n)*(hat(b)[MM] - b)), main = "Loi limite MM - b")
curve(dnorm(x, mean = 0, sd = (b - a) * sqrt(2 / 15)), add = TRUE, col = "red", lwd = 2)
par(mfrow = c(1, 1))


############## Loi Gamma (a = forme, b = taux)  ###############################
n <- 1000
a <- 3
b <- 2

# Matrice de covariance asymptotique de (mean(X), S^2) : [[mu2, mu3], [mu3, mu4 - mu2^2]]
Sigma_asympt <- function(dens, lower, upper) {
  m   <- integrate(function(x) x * dens(x), lower, upper)$value
  mu2 <- integrate(function(x) (x - m)^2 * dens(x), lower, upper)$value
  mu3 <- integrate(function(x) (x - m)^3 * dens(x), lower, upper)$value
  mu4 <- integrate(function(x) (x - m)^4 * dens(x), lower, upper)$value
  list(m = m, v = mu2, Sigma = matrix(c(mu2, mu3, mu3, mu4 - mu2^2), 2, 2))
}

# Ecart-type limite par delta-méthode : sqrt(grad' Sigma grad)
sd_limite <- function(grad, Sigma) as.numeric(sqrt(t(grad) %*% Sigma %*% grad))

vec_b_MM1 <- numeric(M)   # cas 1 : b avec a connu
vec_a_MM2 <- numeric(M)   # cas 2 : a (a et b inconnus)
vec_b_MM2 <- numeric(M)   # cas 2 : b (a et b inconnus)
vec_a_MM3 <- numeric(M)   # cas 3 : a avec b connu

for(i in 1:M) {
  X <- rgamma(n, shape = a, rate = b)
  vec_b_MM1[i] <- a / mean(X)              # cas 1 : E[X] = a/b => b = a / E[X]
  vec_a_MM2[i] <- mean(X)^2 / var(X)       # cas 2 : a = E[X]^2 / V(X)
  vec_b_MM2[i] <- mean(X) / var(X)         # cas 2 : b = E[X] / V(X)
  vec_a_MM3[i] <- b * mean(X)              # cas 3 : a = b E[X]
}

# Boxplot pour b (cas 1 et 2)
quartz()
boxplot(cbind(vec_b_MM1, vec_b_MM2),
        names = c(expression(hat(b)[MM1]), expression(hat(b)[MM2])),
        main = "Gamma(3, 2) : estimation de b")
abline(h = b, col = "red", lwd = 2)

# Boxplot pour a (cas 2 et 3)
quartz()
boxplot(cbind(vec_a_MM2, vec_a_MM3),
        names = c(expression(hat(a)[MM2]), expression(hat(a)[MM3])),
        main = "Gamma(3, 2) : estimation de a")
abline(h = a, col = "red", lwd = 2)

c(Gamma_b_MM1 = mean(abs(vec_b_MM1 - b) > epsilon),
  Gamma_b_MM2 = mean(abs(vec_b_MM2 - b) > epsilon),
  Gamma_a_MM2 = mean(abs(vec_a_MM2 - a) > epsilon),
  Gamma_a_MM3 = mean(abs(vec_a_MM3 - a) > epsilon))

# Variances limites (cas 2 par delta-méthode, avec gradients de a = m^2/v et b = m/v)
res <- Sigma_asympt(function(x) dgamma(x, shape = a, rate = b), 0, Inf)
m0 <- res$m
v0 <- res$v
sd_a2 <- sd_limite(c(2 * m0 / v0, -m0^2 / v0^2), res$Sigma)
sd_b2 <- sd_limite(c(1 / v0, -m0 / v0^2), res$Sigma)

quartz()
par(mfrow = c(2, 2))
hist(sqrt(n) * (vec_b_MM1 - b), breaks = 25, freq = FALSE,
     xlab = expression(sqrt(n)*(hat(b)[MM1] - b)), main = "Loi limite - b (cas 1)")
curve(dnorm(x, mean = 0, sd = b / sqrt(a)), add = TRUE, col = "red", lwd = 2)
hist(sqrt(n) * (vec_b_MM2 - b), breaks = 25, freq = FALSE,
     xlab = expression(sqrt(n)*(hat(b)[MM2] - b)), main = "Loi limite - b (cas 2)")
curve(dnorm(x, mean = 0, sd = sd_b2), add = TRUE, col = "red", lwd = 2)
hist(sqrt(n) * (vec_a_MM2 - a), breaks = 25, freq = FALSE,
     xlab = expression(sqrt(n)*(hat(a)[MM2] - a)), main = "Loi limite - a (cas 2)")
curve(dnorm(x, mean = 0, sd = sd_a2), add = TRUE, col = "red", lwd = 2)
hist(sqrt(n) * (vec_a_MM3 - a), breaks = 25, freq = FALSE,
     xlab = expression(sqrt(n)*(hat(a)[MM3] - a)), main = "Loi limite - a (cas 3)")
curve(dnorm(x, mean = 0, sd = sqrt(a)), add = TRUE, col = "red", lwd = 2)
par(mfrow = c(1, 1))


############## Loi Bêta  ######################################################
n <- 1000
a <- 2
b <- 5

vec_b_MM1 <- numeric(M)   # cas 1 : b avec a connu
vec_a_MM2 <- numeric(M)   # cas 2 : a avec b connu
vec_a_MM3 <- numeric(M)   # cas 3 : a (a et b inconnus)
vec_b_MM3 <- numeric(M)   # cas 3 : b (a et b inconnus)

for(i in 1:M) {
  X <- rbeta(n, shape1 = a, shape2 = b)
  s_hat <- mean(X) * (1 - mean(X)) / var(X) - 1     # s = a + b
  vec_b_MM1[i] <- a * (1 / mean(X) - 1)             # cas 1 : b = a (1/E[X] - 1)
  vec_a_MM2[i] <- b * mean(X) / (1 - mean(X))       # cas 2 : a = b E[X] / (1 - E[X])
  vec_a_MM3[i] <- mean(X) * s_hat                   # cas 3 : a = E[X] s
  vec_b_MM3[i] <- (1 - mean(X)) * s_hat             # cas 3 : b = (1 - E[X]) s
}

# Boxplot pour b (cas 1 et 3)
quartz()
boxplot(cbind(vec_b_MM1, vec_b_MM3),
        names = c(expression(hat(b)[MM1]), expression(hat(b)[MM3])),
        main = "Bêta(2, 5) : estimation de b")
abline(h = b, col = "red", lwd = 2)

# Boxplot pour a (cas 2 et 3)
quartz()
boxplot(cbind(vec_a_MM2, vec_a_MM3),
        names = c(expression(hat(a)[MM2]), expression(hat(a)[MM3])),
        main = "Bêta(2, 5) : estimation de a")
abline(h = a, col = "red", lwd = 2)

c(Beta_b_MM1 = mean(abs(vec_b_MM1 - b) > epsilon),
  Beta_b_MM3 = mean(abs(vec_b_MM3 - b) > epsilon),
  Beta_a_MM2 = mean(abs(vec_a_MM2 - a) > epsilon),
  Beta_a_MM3 = mean(abs(vec_a_MM3 - a) > epsilon))

# Variances limites : cas 1 et 2 par delta-méthode simple, cas 3 par delta-méthode sur (mean, S^2)
sd_b1 <- sqrt(b * (a + b)^2 / (a * (a + b + 1)))
sd_a2 <- sqrt(a * (a + b)^2 / (b * (a + b + 1)))

res <- Sigma_asympt(function(x) dbeta(x, shape1 = a, shape2 = b), 0, 1)
m0 <- res$m
v0 <- res$v
sd_a3 <- sd_limite(c((2 * m0 - 3 * m0^2) / v0 - 1, -m0^2 * (1 - m0) / v0^2), res$Sigma)
sd_b3 <- sd_limite(c((1 - m0) * (1 - 3 * m0) / v0 + 1, -m0 * (1 - m0)^2 / v0^2), res$Sigma)

quartz()
par(mfrow = c(2, 2))
hist(sqrt(n) * (vec_b_MM1 - b), breaks = 25, freq = FALSE,
     xlab = expression(sqrt(n)*(hat(b)[MM1] - b)), main = "Loi limite - b (cas 1)")
curve(dnorm(x, mean = 0, sd = sd_b1), add = TRUE, col = "red", lwd = 2)
hist(sqrt(n) * (vec_a_MM2 - a), breaks = 25, freq = FALSE,
     xlab = expression(sqrt(n)*(hat(a)[MM2] - a)), main = "Loi limite - a (cas 2)")
curve(dnorm(x, mean = 0, sd = sd_a2), add = TRUE, col = "red", lwd = 2)
hist(sqrt(n) * (vec_a_MM3 - a), breaks = 25, freq = FALSE,
     xlab = expression(sqrt(n)*(hat(a)[MM3] - a)), main = "Loi limite - a (cas 3)")
curve(dnorm(x, mean = 0, sd = sd_a3), add = TRUE, col = "red", lwd = 2)
hist(sqrt(n) * (vec_b_MM3 - b), breaks = 25, freq = FALSE,
     xlab = expression(sqrt(n)*(hat(b)[MM3] - b)), main = "Loi limite - b (cas 3)")
curve(dnorm(x, mean = 0, sd = sd_b3), add = TRUE, col = "red", lwd = 2)
par(mfrow = c(1, 1))
