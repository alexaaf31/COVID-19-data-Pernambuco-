# Skewness
library(plot3Drgl)
library(moments)

rm(list=ls())


#PDF function:
f = function(x,par){
  lambda = par[1]
  beta   = par[2]
  alpha  = par[3]
  phi    = par[4]
  g <- (alpha/phi) * (x/phi)^(alpha-1) * exp(-(x/phi)^alpha)
  G <- 1 - exp(-(x/phi)^alpha)
  h = g*(1-G)^G*(G/(1-G) - log(1-G))
  H = 1 - (1 - G)^G
  pdf = lambda*beta*h*(H)^(beta-1)*(1-H)^(-(beta+1))*exp(-lambda*(H/(1-H))^beta)
  return(pdf)
}


lambda = seq(1,5,.1)
beta   = seq(1.5,5,.1)
x1 <- seq(0.01, 10, le = 2*10^4)
U1=f(x1,c(lambda,1.2,4.5,beta))
min(U1[U1 > 0])
which.min(U1)

G_skew<- function(lambda,beta) {
  x= x1[-c(which.min(U1):(2*10^4))]
  U=f(x,c(lambda,beta, 1.2,4.5))
  Skewness = skewness(U) 
  return(Skewness)}
G_skew<-Vectorize(G_skew, c('lambda','beta'))
Skewness<-outer(lambda,beta,G_skew)
persp3D(lambda,beta, Skewness,
        main=expression(paste(MKFW (lambda,beta,1.2,4.5))),
        xlab="lambda",ylab="beta",zlab="Skewness",
        phi = 20,
        col = "darkgray",
        border = "black")

# Kurtosis

# Clear workspace and load required packages
rm(list=ls())
library(plot3Drgl)
library(moments)

#PDF function:
f = function(x,par){
  lambda = par[1]
  beta   = par[2]
  alpha  = par[3]
  phi    = par[4]
  g <- (alpha/phi) * (x/phi)^(alpha-1) * exp(-(x/phi)^alpha)
  G <- 1 - exp(-(x/phi)^alpha)
  #h = g*(1-G)^G*(G/(1-G) - log(1-G))
  #H = 1 - (1 - G)^G
  pdf = lambda*beta*g*(1-G)^(-beta*G)*(1-(1-G)^G)^(beta-1)*(G/(1-G) - log(1-G))*exp(-lambda*((1-G)^(-G)-1)^beta)
  return(pdf)
}

# Define parameter ranges
lambda_range <- seq(1, 5, by = 0.1)
beta_range   <- seq(1.5, 5, by = 0.1)

# Generate x values for PDF evaluation
x <- seq(0.01, 10, length.out = 3 * 10^4)

# Evaluate PDF and find minimum positive value
U <- f(x, c(lambda_range[1], beta_range[1], 0.5,1.5))
min_U <- min(U[U > 0])

# Find index of minimum positive value
min_index <- which(U == min_U)

# Define function to calculate kurtosis
G_kurt <- function(lambda,beta) {
  x <- x[-seq(min_index, 3*10^4)]
  U <- f(x, c(lambda, beta, 0.5, 1.5))
  Kurtosis <- kurtosis(U)
  return(Kurtosis)
}

# Vectorize kurtosis function and calculate kurtosis values for all parameter combinations
G_kurt <- Vectorize(G_kurt, c('lambda', 'beta'))
Kurtosis <- outer(lambda_range, beta_range, G_kurt)

# Generate 3D plot of kurtosis values
persp3D(lambda_range, beta_range, Kurtosis,
        main = expression(paste(MKFW(lambda, beta, 0.5, 1.5))),
        xlab = "lambda", ylab = "beta", zlab = "Kurtosis",
        phi = 20, col = "darkgray", border = "black")


