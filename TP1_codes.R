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

###########Exercice 2 #############

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


rm(list = ls())
set.seed(2026)

M <- 1000
epsilon <- 0.1

rm(list = ls())
set.seed(2026)

M <- 1000
epsilon <- 0.1

#rm(list = ls())
set.seed(2026)

M <- 1000
epsilon <- 0.1



###############################################################################
#           TP1 Appliquer à d'autres lois
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
