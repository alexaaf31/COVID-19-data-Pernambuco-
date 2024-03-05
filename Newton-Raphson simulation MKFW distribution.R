## script to generate random samples (inverion method) using newton method

#parameters
lambda <- 1.2
beta   <- 0.5
alpha  <- 0.8
phi    <- 0.3

my_pdf <- function(x) {
  g <- (alpha/phi) * (x/phi)^(alpha-1) * exp(-(x/phi)^alpha)
  G <- 1 - exp(-(x/phi)^alpha)
  H <- 1 - (1 - G)^G
  h <- g*(1 - G)^G*(G/(1-G) - log(1-G))
  MKNFW <- lambda * beta*h*(H)^(beta-1)*(1 - H)^(-(beta + 1)) * exp(-lambda*(H/(1-H))^beta)
  return(MKNFW)
}

my_cdf <- function(x){
  G <- 1 - exp(-(x/phi)^alpha)
  cdfISW <- 1 -exp(-lambda * ((1-G)^(-G)-1)^beta)
  return(cdfISW)
}

# Apply the Newton-Raphson method
my_inv_cdf <- function(p) {
  # Initialize variables
  x0  <- 0
  x1  <- 0.1
  eps <- 1e-8
  
  # Apply the Newton-Raphson method
  while (abs(x1 - x0) > eps) {
    x0 <- x1
    # Check for invalid input values
    if (x0 < 0 || is.na(my_cdf(x0)) || is.na(my_pdf(x0))) {
      next
    }
    
    x1 <- x0 - (my_cdf(x0) - p) / my_pdf(x0)
    
    # Check if x1 is negative and set it to a small positive value if it is
    if (x1 < 0) {
      x1 <- eps
    }
  }
  return(x1)
}

# Generate a sample of size n from the distribution
n <- 1000
y <- replicate(n, my_inv_cdf(runif(1)))
y
# Plot the histogram of the sample
hist(y)

################################################################################
###############################graphically evaluate#############################
################################################################################
plot(ecdf(y),ylab = "Fn(y)", xlab = "y",main = "")
curve(my_cdf(x), add = TRUE, col = 2,lwd = 2) 
legend(0.2,0.2,lwd = c(1),lty = c(1),
       c("empirical cdf","ISW cdf"),
       col = c("black","red"),
       bty = "n")
################################################################################
#####################################pdf########################################
################################################################################
hist(y, probability =  TRUE, border = "black",col = "white",main = "",)
curve(my_pdf(x), add = TRUE, col = 2, lwd = 2)
legend("topright",lwd = c(1),lty = c(1),
       c("histogram","ISW pdf"),
       col = c("black","red"),
       bty ="n")
################################################################################

MCS = function(nrep=1000, nobs, semente=2011, lambda = 1.2, beta = 0.5, alpha = 0.8,
               phi = 0.3)
{ 
  tempo.inicio = Sys.time()
  logLikk = function(param){
    lambda = param[1]; beta = param[2]; alpha = param[3]; phi = param[4]; nobs = length(y)
    g <- (alpha/phi) * (y/phi)^(alpha-1) * exp(-(y/phi)^alpha)
    G <- 1 - exp(-(y/phi)^alpha)
    loglik <- log(lambda) + log(beta) + log(g) -(beta*G)*log(1-G) + (beta-1)*log(1 - (1 - G)^G) + 
      log(G/(1-G) - log(1 - G)) -lambda * ((1-G)^(-G)-1)^beta
    minus.logL <- sum(-loglik, na.rm = F)
    return(minus.logL)
  }
  emvlambda <- rep(0, nrep)
  emvbeta   <- rep(0, nrep)
  emvalpha  <- rep(0,nrep)
  emvphi    <- rep(0,nrep)
  set.seed(semente)
  contadorFalhas <- 0
  # laco de Monte Carlo
  i = 0
  
  while(i < nrep){
    
    y <- replicate(nobs, my_inv_cdf(runif(1)))
    
    # Optimize log-likelihood function
    
    ir <- suppressWarnings(optim(c(lambda,beta,alpha,phi),logLikk, method="BFGS"))
    if(ir$convergence == 0){
      i <- i + 1
      emvlambda[i]<- ir$par[1]
      emvbeta[i]  <- ir$par[2]
      emvalpha[i] <- ir$par[3]
      emvphi[i]   <- ir$par[4]
      
    }
    else{
      contadorFalhas <- contadorFalhas + 1
    }
  } 
  # fim do laco de Monte Carlo
  lambdamedio= mean(emvlambda)
  betamedio  = mean(emvbeta)
  alphamedio = mean(emvalpha)
  phimedio   = mean(emvphi)
  lambdavies = lambdamedio - lambda
  betavies   = betamedio   - beta
  alphavies  = alphamedio  - alpha
  phivies    = phimedio    - phi
  tempo.fim  = Sys.time()
  tempo.exec = tempo.fim - tempo.inicio
  lambdaEQM  = var(emvlambda) + (lambdavies)^2
  betaEQM    = var(emvbeta)   + (betavies)^2
  alphaEQM   = var(emvalpha)  + (alphavies)^2
  phiEQM     = var(emvphi)    + (phivies)^2
  resultado  = data.frame(nobs=nobs, nrep=nrep, semente=semente,lambda=lambda, beta=beta,
                          alpha=alpha,phi=phi,lambdamedio = lambdamedio, betamedio = betamedio, 
                          alphamedio=alphamedio,phimedio=phimedio,lambdavies = lambdavies,betavies=betavies, 
                          alphavies=alphavies,phivies=phivies, lambdaEQM= lambdaEQM,betaEQM=betaEQM, 
                          alphaEQM = alphaEQM,phiEQM=phiEQM, falhas=contadorFalhas, 
                          tempo=tempo.exec)
  return(resultado)
}
MCS(nobs=50)
MCS(nobs=100)
MCS(nobs=200) 
MCS(nobs=500)
