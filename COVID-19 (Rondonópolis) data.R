
######################################################################################################################
######################################################################################################################
######################################################################################################################

# Function to calculate information criteria (AIC, BIC, CAIC, HQIC)
criterios <- function(value, num_params, num_data_points) {
  # value - the value from optimization
  # num_params - the number of parameters
  # num_data_points - the number of data points
  
  # Calculate log-likelihood
  l <- 2 * value
  
  # Calculate AIC
  AIC <- l + 2 * num_params
  
  # Calculate BIC
  BIC <- l + num_params * log(num_data_points)
  
  # Calculate CAIC
  CAIC <- AIC + (2 * (num_params + 2) * (num_params + 3)) / (num_data_points - num_params - 3)
  
  # Calculate HQIC
  HQIC <- l + 2 * log(log(num_data_points)) * num_params
  
  # Combine results into a matrix
  result <- cbind(AIC, CAIC, BIC, HQIC)
  
  # Return results
  return(result)
}



# Function to calculate confidence intervals and p-values of parameters
IC <- function(parametros, hessiana, n) {
  # Calculate variance-covariance matrix
  Var <- solve(hessiana, tol = 1e-15)
  
  # Calculate standard errors of parameters
  SE <- sqrt(diag(Var))
  
  # Calculate t-values of parameters
  tvalue <- parametros / SE
  
  # Calculate p-values of parameters
  pvalue <- 2 * pt(abs(tvalue), n, lower.tail = FALSE)
  
  # Calculate confidence intervals of parameters
  LI <- parametros - qt(0.975, n) * SE
  LS <- parametros + qt(0.975, n) * SE
  
  # Combine results into a matrix
  resul <- cbind(parametros, SE, tvalue, pvalue, LI, LS)
  colnames(resul) <- c("Parameter", "SE", "t-value", "p-value", "Lower Bound", "Upper Bound")
  
  # Return results
  return(resul)
}

setwd("D:\\Syncmega\\Sync\\Doutorado UFPE\\Tese\\tese_escrita_inglês\\tese_NFMK-G\\aplicação\\data")

require("csv")
require("xlsx")
library(utils)
require(survival)
require(AdequacyModel)
require(zipfR)

dados = read.csv("RONCOVID.csv",1)
t = dados$time
hist(t,border = "black",col = "white", xlab = "age ( in years)",main = "")
y=log(t)
censur = dados$censur
mean(dados$censur==0)
descriptive(t)
hist(dados$idade,border = "black",col = "white", xlab = "age (in years)",main = "")
#########################kaplan-Meier from variable acomp######################
ekm=survfit(Surv(t)~dados$Cardiovascular)
plot(ekm,col = c("black","darkred","darkgreen"),xlim = c(0,50), xlab = " Survival time (days)", 
     ylab = "Survival probability ",lwd=c(2,2), lty=c(1:2))
legend("topright", c("0","1"), lty = 1:2, lwd = c(2,2), col = c("black","darkred","darkgreen"),bty="n")

#######################################################################################
##################################LMKNFW###############################################
#######################################################################################
MKNFW<-function(par){
  lambda= par[1]
  beta  = par[2]
  sigma = par[3]
  b0    = par[4]
  b1    = par[5]
  b2    = par[6]
  if(any(lambda< 1e-20)) return(.Machine$double.xmax^.5)
  if(any(beta  < 1e-20)) return(.Machine$double.xmax^.5)
  if(any(sigma < 1e-20)) return(.Machine$double.xmax^.5)
  mu=b0 + b1*dados$idade + b2*dados$Cardiovascular
  z = (y-mu)/(sigma)
  g = (1/sigma)*exp(z)*exp(-exp(z))
  G = 1 - exp(-exp(z))
  f = lambda*beta*g*(1-G)^(-beta*G)*(1 - (1 - G)^G)^(beta-1)*(G/(1-G) - log(1-G))*
    exp(-lambda*((1 - G)^(-G) - 1)^beta)
  s = exp(-lambda*((1 - G)^(-G) - 1)^beta)
  
  # Calculate the log-likelihood for each observation using censored data
  lv <- censur*(log(f))+(1-censur)*(log(s))
  
  # Return the negative sum of the log-likelihoods
  return(sum(-lv))
}
#valorin<-c(0.66569563,  1.92998915,  2.18748308,  4.34484597, -0.01692673,-0.26912549)
valorin<-c(2.40369648,  0.22824899,  0.24063441,  5.09918629, -0.01635904, -0.30738264)
LMKNFW<-optim(valorin,MKNFW,method="BFGS",hessian=T);
LMKNFW
n=nrow(dados)
criterios(LMKNFW$value,length(LMKNFW$par),n)
parametros=LMKNFW$par
hessiana=LMKNFW$hessian
IC(parametros,hessiana,n)

#######################################################################################
##################################LKWW#################################################
#######################################################################################
KWW<-function(par){
  a     = par[1]
  b     = par[2]
  sigma = par[3]
  b0    = par[4]
  b1    = par[5]
  b2    = par[6]
  if(any(a     < 1e-20)) return(.Machine$double.xmax^.5)
  if(any(b     < 1e-20)) return(.Machine$double.xmax^.5)
  if(any(sigma < 1e-20)) return(.Machine$double.xmax^.5)
  mu=b0 + b1*dados$idade + b2*dados$Cardiovascular
  z = (y-mu)/(sigma)
  g = (1/sigma)*exp(z)*exp(-exp(z))
  G = 1 - exp(-exp(z))
  f = a*b*g*G^(a-1)*(1 - G^a)^(b-1)
  s = (1 - G^a)^b
  lv<-censur*(log(f))+(1-censur)*(log(s))
  sum(-lv)
}

#valorin<-c(2.17631661,  0.16619760,  0.68421240,  3.35664936, -0.01710062, -0.24344652)
valorin<-c(1.3310590,  0.5509565,  0.6324386,  4.1835414, -0.0166783, -0.3029760)
LKWW<-optim(valorin,KWW,method="BFGS",hessian=T);
LKWW
n=nrow(dados)
criterios(LKWW$value,length(LKWW$par),n)
parametros1=LKWW$par
hessiana1=LKWW$hessian
IC(parametros1,hessiana1,n)

#######################################################################################
####################################LWW################################################
#######################################################################################
WW<-function(par){
  alpha = 1
  lambda= par[1]
  beta  = par[2]
  sigma = par[3]
  b0    = par[4]
  b1    = par[5]
  b2    = par[6]
  if(any(alpha < 1e-20)) return(.Machine$double.xmax^.5)
  if(any(lambda< 1e-20)) return(.Machine$double.xmax^.5)
  if(any(beta  < 1e-20)) return(.Machine$double.xmax^.5)
  if(any(sigma < 1e-20)) return(.Machine$double.xmax^.5)
  mu=b0 + b1*dados$idade+ b2*dados$Cardiovascular
  z = (y-mu)/(sigma)
  g = (1/sigma)*exp(z)*exp(-exp(z))
  G = 1 - exp(-exp(z))
  f = alpha*lambda*beta*g*G^(alpha*beta - 1)*(1-G)^(-(alpha*beta + 1))*
    exp(-lambda*(G/(1-G))^(alpha*beta))
  s = exp(-lambda*(G/(1-G))^(alpha*beta))
  lv<-censur*(log(f))+(1-censur)*(log(s))
  sum(-lv)
}
valorin<-c(0.03047366,  4.02213124,  3.44316399,  3.93440796, -0.01602248, -0.26759115)
LWW<-optim(valorin,WW,method="BFGS",hessian=T);
LWW
n=nrow(dados)
criterios(LWW$value,length(LWW$par),n)
parametros2=LWW$par
hessiana2=LWW$hessian
IC(parametros2,hessiana2,n)

#######################################################################################
####################################LBW################################################
#######################################################################################
BW<-function(par){
  a     = par[1]
  b     = par[2]
  sigma = par[3]
  b0    = par[4]
  b1    = par[5]
  b2    = par[6]
  if(any(a < 1e-20)) return(.Machine$double.xmax^.5)
  if(any(b< 1e-20)) return(.Machine$double.xmax^.5)
  if(any(sigma < 1e-20)) return(.Machine$double.xmax^.5)
  mu=b0 + b1*dados$idade + b2*dados$Cardiovascular
  z = (y-mu)/(sigma)
  g = (1/sigma)*exp(z)*exp(-exp(z))
  G = 1 - exp(-exp(z))
  f = (1/(beta(a,b)))*g*G^(a-1)*(1 - G)^(b-1)
  s = 1 - pbeta(G,a,b)
  lv<-censur*(log(f))+(1-censur)*(log(s))
  sum(-lv)
}
valorin<-c(1.40132291,  1.15283214,  0.67167665,  4.57461427, -0.01680173, -0.29600621)
LBW<-optim(valorin,BW,method="BFGS",hessian=T);
LBW
n=nrow(dados)
criterios(LBW$value,length(LBW$par),n)
parametros3=LBW$par
hessiana3=LBW$hessian
IC(parametros3,hessiana3,n)
################################################################################



############################residual quantile and QQ-plot#######################

FMKNFW<-function(par,y){
  lambda= par[1]
  beta  = par[2]
  sigma = par[3]
  b0    = par[4]
  b1    = par[5]
  b2    = par[6]
  if(any(lambda< 1e-20)) return(.Machine$double.xmax^.5)
  if(any(beta  < 1e-20)) return(.Machine$double.xmax^.5)
  if(any(sigma < 1e-20)) return(.Machine$double.xmax^.5)
  mu=b0 + b1*dados$idade + b2*dados$Cardiovascular
  z = (y-mu)/(sigma)
  G = 1 - exp(-exp(z))
  cdf = 1 - exp(-lambda*((1 - G)^(-G) - 1)^beta)
  return(cdf)
}

#########################residual index without package#########################
plot(qnorm(FMKNFW(LMKNFW$par,y)),
     pch = 16,
     ylab = "Quantile residuals",
     xlab = "Index",ylim = c(-4,4))

# Add horizontal lines at y = -3 and y = 3
abline(h = c(-3, 0, 3), col = "black", lty = 2)

# Add horizontal line at y = 0
abline(h = 0, col = "black",lty = 2)
##########################QQ-plot normal distribution###########################
require(car)

qqPlot(qnorm(FMKNFW(LMKNFW$par, y)),ylab = "Sample Quantiles", 
       xlab = "Theoretical Quantiles", envelope = F, grid = F, 
       pch =16,col.lines=carPalette()[1],id=F)

###########################second form##########################################

# Set seed for reproducibility
set.seed(77)

# Generate random normal data
q <- rnorm(1000)

# Pre-calculate qnorm_y
qnorm_y <- qnorm(FMKNFW(LMKNFW$par, y))

# Generate QQ plot
qqplot(q, qnorm_y,
       ylab = "Sample Quantiles",
       xlab = "Theoretical Quantiles",
       pch = 16,
       cex.axis = 1.2)

# Add QQ line to the plot
qqline(qnorm_y, lwd = 3)

#teste da razao de verosValuelhancas generalizado

pdf_MKNFW<-function(par,y){
  lambda= par[1]
  beta  = par[2]
  sigma = par[3]
  b0    = par[4]
  b1    = par[5]
  b2    = par[6]
  if(any(lambda< 1e-20)) return(.Machine$double.xmax^.5)
  if(any(beta  < 1e-20)) return(.Machine$double.xmax^.5)
  if(any(sigma < 1e-20)) return(.Machine$double.xmax^.5)
  mu=b0 + b1*dados$idade + b2*dados$Cardiovascular
  z = (y-mu)/(sigma)
  g = (1/sigma)*exp(z)*exp(-exp(z))
  G = 1 - exp(-exp(z))
  f = lambda*beta*g*(1-G)^(-beta*G)*(1 - (1 - G)^G)^(beta-1)*(G/(1-G) - log(1-G))*
    exp(-lambda*((1 - G)^(-G) - 1)^beta)
  return(f)
}

pdf_KWW<-function(par,y){
  a     = par[1]
  b     = par[2]
  sigma = par[3]
  b0    = par[4]
  b1    = par[5]
  b2    = par[6]
  if(any(a     < 1e-20)) return(.Machine$double.xmax^.5)
  if(any(b     < 1e-20)) return(.Machine$double.xmax^.5)
  if(any(sigma < 1e-20)) return(.Machine$double.xmax^.5)
  mu=b0 + b1*dados$idade + b2*dados$Cardiovascular
  z = (y-mu)/(sigma)
  g = (1/sigma)*exp(z)*exp(-exp(z))
  G = 1 - exp(-exp(z))
  f = a*b*g*G^(a-1)*(1 - G^a)^(b-1)
  return(f)
}

pdf_WW<-function(par,y){
  alpha = 1
  lambda= par[1]
  beta  = par[2]
  sigma = par[3]
  b0    = par[4]
  b1    = par[5]
  b2    = par[6]
  if(any(alpha < 1e-20)) return(.Machine$double.xmax^.5)
  if(any(lambda< 1e-20)) return(.Machine$double.xmax^.5)
  if(any(beta  < 1e-20)) return(.Machine$double.xmax^.5)
  if(any(sigma < 1e-20)) return(.Machine$double.xmax^.5)
  mu=b0 + b1*dados$idade + b2*dados$Cardiovascular
  z = (y-mu)/(sigma)
  g = (1/sigma)*exp(z)*exp(-exp(z))
  G = 1 - exp(-exp(z))
  f = alpha*lambda*beta*g*G^(alpha*beta - 1)*(1-G)^(-(alpha*beta + 1))*
    exp(-lambda*(G/(1-G))^(alpha*beta))
  return(f)
}

pdf_BW<-function(par,y){
  a     = par[1]
  b     = par[2]
  sigma = par[3]
  b0    = par[4]
  b1    = par[5]
  b2    = par[6]
  if(any(a < 1e-20)) return(.Machine$double.xmax^.5)
  if(any(b< 1e-20)) return(.Machine$double.xmax^.5)
  if(any(sigma < 1e-20)) return(.Machine$double.xmax^.5)
  mu=b0 + b1*dados$idade + b2*dados$Cardiovascular
  z = (y-mu)/(sigma)
  g = (1/sigma)*exp(z)*exp(-exp(z))
  G = 1 - exp(-exp(z))
  f = (1/(beta(a,b)))*g*G^(a-1)*(1 - G)^(b-1)
  return(f)
}


#MKNFW vs KW
n <- length(t)
log_ratio <- log(pdf_MKNFW(par = LMKNFW$par, y) / pdf_KWW(par = LKWW$par, y))

RL1 <- (1 / sqrt(n)) * sum(log_ratio)
W1  <- (1 / n) * sum(log_ratio^2)
S1  <- ((1 / n) * sum(log_ratio))^2
RLNN1 <- RL1 / (W1 - S1)

RLNN1

#MKNFW vs WW
log_ratio1 <- log(pdf_MKNFW(par = LMKNFW$par, y) / pdf_WW(par = LWW$par, y))

RL2 <- (1 / sqrt(n)) * sum(log_ratio1)
W2  <- (1 / n) * sum(log_ratio1^2)
S2  <- ((1 / n) * sum(log_ratio1))^2
RLNN2 <- RL2 / (W2 - S2)

RLNN2

#MKNFW vs BW
log_ratio2 <- log(pdf_MKNFW(par = LMKNFW$par, y) / pdf_BW(par = LBW$par, y))

RL3 <- (1 / sqrt(n)) * sum(log_ratio2)
W3  <- (1 / n) * sum(log_ratio2^2)
S3  <- ((1 / n) * sum(log_ratio2))^2
RLNN3 <- RL3 / (W3 - S3)

RLNN3

