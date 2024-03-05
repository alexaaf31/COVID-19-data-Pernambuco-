# Set parameters
lambda <- 0.5
beta   <- 0.4
sigma  <- 4.0
b0     <- 1.5
b1     <- 2.5

# Generate a random number between 0 and 1
x1 <- runif(1)

# Calculate mu
mu <- b0 + b1 * x1

# Define the MKNFW distribution PDF
f <- function(x) {
  z <- (x - mu) / sigma
  g <- (1 / sigma) * exp(z - exp(z))
  G <- 1 - exp(-exp(z))
  H <- 1 - (1 - G)^G
  h <- g * (1 - G)^G * (G / (1 - G) - log(1 - G))
  MKNFW <- lambda * beta * h * (H)^(beta - 1) * (1 - H)^(-(beta + 1)) * exp(-lambda * (H / (1 - H))^beta)
  return(MKNFW)
}

# Define the g(x) PDF
g <- function(x) {
  z <- (x - mu) / sigma
  pdfW <- (1 / sigma) * exp(z - exp(z))
  return(pdfW)
}

# Define the g(x) quantile function
rLWeibull <- function(n) {
  u <- runif(n)
  Q <- sigma * (log(-log(1 - u))) + mu
  return(Q)
}

# Set the random seed
set.seed(1986)

# Find M
t <- rLWeibull(1000)
M <- max(f(t) / g(t))

# Define the rLMKNFW function to generate random samples from the MKNFW distribution
rLMKNFW <- function(f, g, M, n) {
  naccepts <- 0
  y <- numeric(n)
  while (naccepts < n) {
    x <- rLWeibull(1)
    u <- runif(1)
    w <- f(x) / (M * g(x))
    if (u <= w) {
      y[naccepts + 1] <- x
      naccepts <- naccepts + 1
    }
  }
  return(y)
}

# Generate 1000 random samples from the MKNFW distribution using the rLMKNFW function
y <- rLMKNFW(f, g, M, 1000)

y=d=pc2=NULL
u = 11.5
nobs = 50

# Generate random values from the rLMKNFW distribution with parameters f, g, M, and nobs
t = rLMKNFW(f, g, M, nobs) 

# Generate uniform random values between 0 and u
c = runif(nobs, 0, u) 

# Initialize y and d vectors with zeros using the numeric() function
y = numeric(nobs) 
d = numeric(nobs) 

# Generate y and d values based on the comparison of t and c
for (i in 1:nobs) {
  if (t[i] <= c[i]) { # if t is less than or equal to c
    y[i] = t[i] # set y to t
    d[i] = 1 # set d to 1
  } else { # if t is greater than c
    y[i] = c[i] # set y to c
    d[i] = 0 # set d to 0
  }
}

# Combine y and d vectors into a matrix using cbind() function
Dc = cbind(y, d) 

# Calculate the proportion of d values equal to 0
pc2 = mean(d == 0) 

# Return the Dc matrix and pc2 proportion as a list
list(Dc = Dc, pc2 = pc2)

################################################################################
#############################Monte Carlo Simulation (MCS)#######################
################################################################################
MCS = function(nrep=1000, nobs, semente=1986, lambda = 0.5, beta = 0.4, sigma = 4.0,
               b0 = 1.5, b1 = 2.5)
{ 
  D = function(u, f, g, M){
    
    # Generate random values from the rLMKNFW distribution with parameters f, g, M, and nobs
    t = rLMKNFW(f, g, M, nobs) 
    
    # Generate uniform random values between 0 and u
    c = runif(nobs, 0, u) 
    
    # Initialize y and d vectors with zeros using the numeric() function
    y = numeric(nobs) 
    d = numeric(nobs) 
    
    # Generate y and d values based on the comparison of t and c
    for (i in 1:nobs) {
      if (t[i] <= c[i]) { # if t is less than or equal to c
        y[i] = t[i] # set y to t
        d[i] = 1 # set d to 1
      } else { # if t is greater than c
        y[i] = c[i] # set y to c
        d[i] = 0 # set d to 0
      }
    }
    
    # Combine y and d vectors into a matrix using cbind() function
    Dc = cbind(y, d) 
    
    # Calculate the proportion of d values equal to 0
    pc2 = mean(d == 0) 
    
    # Return the Dc matrix and pc2 proportion as a list
    return(list(Dc = Dc, pc2 = pc2)) 
  } 
  
  tempo.inicio = Sys.time()
  logLikk = function(param, D){
    # Extract parameters from param vector
    lambda = param[1]
    beta   = param[2]
    sigma  = param[3]
    b0     = param[4]
    b1     = param[5]
    
    nobs = length(y) # Calculate number of observations
    
    # Check if any parameter is smaller than 1e-20, return large value if true
    if(any(param < 1e-20)) return(.Machine$double.xmax^.5)
    
    # Calculate mu, z, g, and G using extracted parameters
    mu = b0 + b1*x1
    z = (y-mu)/sigma
    g <- (1/sigma)*exp(z - exp(z))
    G <- 1 - exp(-exp(z))
    
    # Calculate the MKNFW density function f and the survival function s
    f = lambda*beta*g*(1-G)^(-beta*G)*(1 - (1 - G)^G)^(beta-1)*(G/(1-G) - log(1-G))*
      exp(-lambda*((1 - G)^(-G) - 1)^beta)
    s = exp(-lambda*((1 - G)^(-G) - 1)^beta)
    
    # Calculate the log-likelihood for each observation using censored data
    lv <- censur*(log(f))+(1-censur)*(log(s))
    
    # Return the negative sum of the log-likelihoods
    return(sum(-lv))
  }
  
  emvlambda <- emvbeta <- emvsigma <- emvb0 <- emvb1 <- numeric(nrep)
  
  set.seed(semente)
  contadorFalhas = 0
  # Start of the Monte Carlo Loop
  i = 0
  
  while(i < nrep){
    
    Dc=D(u, f, g, M)# dados completos
    Dados=data.frame(Dc[1]) #transformando em data.frame
    names(Dados)=c("temp","censur")
    Dobs=Dados[,]#dados observaveis sem variaveis latentes
    y = Dobs$temp
    censur = Dobs$censur
    
    # Optimize log-likelihood function
    ir =  suppressWarnings(optim(c(lambda,beta,sigma,b0,b1),logLikk, method="Nelder-Mead"))
    if(ir$convergence == 0){
      i = i + 1
      emvlambda[i]= ir$par[1]
      emvbeta[i]  = ir$par[2]
      emvsigma[i] = ir$par[3]
      emvb0[i]    = ir$par[4]
      emvb1[i]    = ir$par[5]
      
    }
    else{
      contadorFalhas = contadorFalhas + 1
    }
  } 
  # End of the Monte Carlo Loop
  lambdamedio= mean(emvlambda)
  betamedio  = mean(emvbeta)
  sigmamedio = mean(emvsigma)
  b0medio    = mean(emvb0)
  b1medio    = mean(emvb1)
  lambdavies = lambdamedio - lambda
  betavies   = betamedio   - beta
  sigmavies  = sigmamedio  - sigma
  b0vies     = b0medio     - b0
  b1vies     = b1medio     - b1
  tempo.fim  = Sys.time()
  tempo.exec = tempo.fim - tempo.inicio
  lambdaEQM  = var(emvlambda) + (lambdavies)^2
  betaEQM    = var(emvbeta)   + (betavies)^2
  sigmaEQM   = var(emvsigma)  + (sigmavies)^2
  b0EQM      = var(emvb0)     + (b0vies)^2
  b1EQM      = var(emvb1)     + (b1vies)^2
  resultado  = data.frame(nobs=nobs, nrep=nrep, semente=semente,lambda=lambda, beta=beta,
                          sigma=sigma,b0=b0,b1=b1,lambdamedio = lambdamedio, betamedio = betamedio, 
                          sigmamedio=sigmamedio,b0medio=b0medio,b1medio=b1medio,lambdavies = lambdavies,betavies=betavies, 
                          sigmavies=sigmavies,b0vies=b0vies,b1vies=b1vies,lambdaEQM= lambdaEQM,betaEQM=betaEQM, 
                          sigmaEQM = sigmaEQM,b0EQM=b0EQM,b1EQM=b1EQM,falhas=contadorFalhas, 
                          tempo=tempo.exec)
  return(resultado)
}
MCS(nobs=50)
MCS(nobs=100)
MCS(nobs=200) 
MCS(nobs=500)

