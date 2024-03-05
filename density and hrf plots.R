################################################################################
#######################density plot of the MKNFW################################
t <- seq(0.001, 10, 0.01)
MKNFW <- function(lambda, beta, alpha, phi) {
  g <- (alpha/phi) * (t/phi)^(alpha-1) * exp(-(t/phi)^alpha)
  G <- 1 - exp(-(t/phi)^alpha)
  h = g*(1-G)^G*(G/(1-G) - log(1-G))
  H = 1 - (1 - G)^G
  f = lambda*beta*h*(H)^(beta-1)*(1-H)^(-(beta+1))*exp(-lambda*(H/(1-H))^beta)
  return(f)
}

p <- plot(c(0,4.5), c(0, 0.55), type = "n", xlab = "x", ylab = "pdf")
lambda <- c(0.2, 0.5, 0.2, 0.5, 1.1)
beta   <- c(0.2, 0.9, 0.1, 0.3, 0.2)
alpha  <- c(1.9, 0.9, 1.8, 2.1, 4.5)
phi    <- c(0.9, 0.9, 0.8, 1.2, 3.0)
cor1 <- c('black', 'darkred', 'darkblue', 'darkgreen', "darkorange")
linhas <- c(1:5)
for (i in 1:length(lambda)) {
  p <- lines(t, MKNFW(lambda[i], beta[i], alpha[i], phi[i]), lwd = c(2, 2, 2, 2, 2), lty = linhas[i], col = cor1[i])
}
legend(2.5,0.57,
       expression(
         paste(lambda, "=0.2; ", beta, "=0.2; ", alpha, "=1.9; ", tau, "=0.9"),
         paste(lambda, "=0.5; ", beta, "=0.9; ", alpha, "=0.9; ", tau, "=0.9"),
         paste(lambda, "=0.2; ", beta, "=0.1; ", alpha, "=1.8; ", tau, "=0.8"),
         paste(lambda, "=0.5; ", beta, "=0.3; ", alpha, "=2.1; ", tau, "=1.2"),
         paste(lambda, "=1.1; ", beta, "=0.2; ", alpha, "=4.5; ", tau, "=3.0")),
       lty = c(1:5), lwd = c(2, 2, 2, 2, 2), col = cor1, bty = "n")

################################################################################
############################hrf plot of the MKNFW###############################
t <- seq(0.001, 10, 0.01)
hMKNFW <- function(lambda, beta, alpha, phi) {
  g <- (alpha/phi) * (t/phi)^(alpha-1) * exp(-(t/phi)^alpha)
  G <- 1 - exp(-(t/phi)^alpha)
  f <- lambda * beta * g * (1-G)^(-beta*G) * (1 - (1 - G)^G)^(beta-1) * (G/(1-G) - log(1-G)) 
  return(f)
}

p <- plot(c(0,3.5), c(0, 0.55), type = "n", xlab = "x", ylab = "hrf")
lambda <- c(0.2, 0.5, 0.2, 0.1, 0.8)
beta   <- c(0.2, 0.1, 0.1, 0.1, 0.5)
alpha  <- c(1.6, 0.9, 1.8, 1.6, 6.0)
phi    <- c(0.9, 0.9, 0.8, 0.9, 4.0)
cor1 <- c('black', 'darkred', 'darkblue', 'darkgreen', "darkorange")
linhas <- c(1:5)
for (i in 1:length(lambda)) {
  p <- lines(t, hMKNFW(lambda[i], beta[i], alpha[i], phi[i]), lwd = c(2, 2, 2, 2, 2), lty = linhas[i], col = cor1[i])
}
legend(0.7,0.55,
       expression(
         paste(lambda, "=0.2; ", beta, "=0.2; ", alpha, "=1.6; ", tau, "=0.9"),
         paste(lambda, "=0.5; ", beta, "=0.1; ", alpha, "=0.9; ", tau, "=0.9"),
         paste(lambda, "=0.2; ", beta, "=0.1; ", alpha, "=1.8; ", tau, "=0.8"),
         paste(lambda, "=0.1; ", beta, "=0.1; ", alpha, "=1.6; ", tau, "=0.9"),
         paste(lambda, "=0.8; ", beta, "=0.5; ", alpha, "=6.0; ", tau, "=4.0")),
       lty = c(1:5), lwd = c(2, 2, 2, 2, 2), col = cor1, bty = "n")

################################################################################
#######################density plot of the MKNFK################################
################################################################################
rm(list = ls())
t<-seq(0.001,1,.01)
MKNFK = function(lambda,beta,b,c){
  g = b*c*t^(b-1)*(1 - t^b)^(c-1)
  G = 1 - (1 - t^b)^c
  f <- lambda * beta * g * (1-G)^(-beta*G) * (1 - (1 - G)^G)^(beta-1) * (G/(1-G) - log(1-G)) *
    exp(-lambda * ((1-G)^(-G)-1)^beta)
  return(f)
}

p=plot(c(0,1), c(0,3), type = "n", xlab="x", ylab="pdf")
lambda = c(0.5,0.6,0.9,1.3,0.9)
beta   = c(0.2,1.1,0.8,0.2,0.25)
b      = c(1.2,0.9,0.9,0.5,2.6)
c      = c(3.5,0.9,3.5,0.5,8.5)
cor1=c('black','darkred','darkblue','darkgreen',"darkorange")
linhas=c(1:5)
for( i in 1:length(lambda)){
  p[i]=lines(t,MKNFK(lambda[i],beta[i],b[i],c[i]),lwd=c(2,2,2,2,2),lty=linhas[i],col=cor1[i])}
legend(0.3,3, expression(
  paste(lambda,"=0.5; ",beta,"=0.2; ", a,"=1.2; ",b,"=3.5"),    
  paste(lambda,"=0.6; ",beta,"=1.1; ", a,"=0.9; ",b,"=0.9"),
  paste(lambda,"=0.9; ",beta,"=0.8; ", a,"=0.9; ",b,"=3.5"),
  paste(lambda,"=1.3; ",beta,"=0.2; ", a,"=0.5; ",b,"=0.5"),
  paste(lambda,"=0.6; ",beta,"=0.25; ", a,"=2.6; ",b,"=9.5")),
  lty=c(1:5), lwd=c(2,2,2,2,2),col=cor1,
  bty="n")
################################################################################
###########################hrf plot of the MKNFK################################
################################################################################
rm(list = ls())
t<-seq(0.001,1,.01)
hMKNFK = function(lambda,beta,b,c){
  g = b*c*t^(b-1)*(1 - t^b)^(c-1)
  G = 1 - (1 - t^b)^c
  f <- lambda * beta * g * (1-G)^(-beta*G) * (1 - (1 - G)^G)^(beta-1) * (G/(1-G) - log(1-G))
  return(f)
}

p=plot(c(0,1), c(0,1.0), type = "n", xlab="x", ylab="hrf")
lambda = c(0.1,0.2,0.2,1.2,0.2)
beta   = c(0.1,1.2,1.9,0.1,0.2)
b      = c(2.5,0.9,0.9,0.5,1.6)
c      = c(3.5,0.9,0.5,0.5,5.5)
cor1=c('black','darkred','darkblue','darkgreen',"darkorange")
linhas=c(1:5)
for( i in 1:length(lambda)){
  p[i]=lines(t,hMKNFK(lambda[i],beta[i],b[i],c[i]),lwd=c(2,2,2,2,2),lty=linhas[i],col=cor1[i])}
legend(0.13,1.0, expression(
  paste(lambda,"=0.1; ",beta,"=0.1; ", a,"=2.5; ",b,"=3.5"),    
  paste(lambda,"=0.2; ",beta,"=1.2; ", a,"=0.9; ",b,"=0.9"),
  paste(lambda,"=0.2; ",beta,"=1.9; ", a,"=0.9; ",b,"=0.5"),
  paste(lambda,"=1.2; ",beta,"=0.1; ", a,"=0.5; ",b,"=0.5"),
  paste(lambda,"=0.2; ",beta,"=0.2; ", a,"=1.6; ",b,"=5.5")),
  lty=c(1:5), lwd=c(2,2,2,2,2),col=cor1,
  bty="n")
################################################################################
#######################density plot of the MKNFN################################
################################################################################
rm(list = ls())
t<-seq(-10,10,.01)
MKNFN = function(lambda,beta,mu,sigma){
  g = dnorm(t,mu, sigma)
  G = pnorm(t, mu, sigma)
  f <- lambda * beta * g * (1-G)^(-beta*G) * (1 - (1 - G)^G)^(beta-1) * (G/(1-G) - log(1-G)) *
    exp(-lambda * ((1-G)^(-G) - 1)^beta)
  return(f)
}

p=plot(c(-1.7,4.5), c(0,0.40), type = "n", xlab="x", ylab="pdf")
lambda = c(1.1,1.6,1.1,0.9,1.2)
beta   = c(0.1,0.1,0.5,0.1,0.6)
mu     = c(0.6,1.8,2.8,1.4,0.1)
sigma  = c(0.5,0.7,1.0,0.7,1.5)
cor1=c('black','darkred','darkblue','darkgreen',"darkorange")
linhas=c(1:5)
for( i in 1:length(lambda)){
  p[i]=lines(t,MKNFN(lambda[i],beta[i],mu[i],sigma[i]),lwd=c(2,2,2,2,2),lty=linhas[i],col=cor1[i])}
legend("topleft", expression(
  paste(lambda,"=1.1; ",beta,"=0.1; ", mu,"=0.6; ",sigma,"=0.5"),    
  paste(lambda,"=1.6; ",beta,"=0.1; ", mu,"=0.8; ",sigma,"=0.7"),
  paste(lambda,"=1.1; ",beta,"=0.5; ", mu,"=2.8; ",sigma,"=1.0"),
  paste(lambda,"=0.9; ",beta,"=0.1; ", mu,"=0.4; ",sigma,"=0.7"),
  paste(lambda,"=1.2; ",beta,"=0.6; ", mu,"=0.1; ",sigma,"=1.5")),
  lty=c(1:5), lwd=c(2,2,2,2,2),col=cor1,
  bty="n")

################################################################################
############################hrf plot of the MKNFN###############################
################################################################################
rm(list = ls())
t<-seq(-10,10,.01)
hMKNFN = function(lambda,beta,mu,sigma){
  g = dnorm(t,mu, sigma)
  G = pnorm(t, mu, sigma)
  f <- lambda * beta * g * (1-G)^(-beta*G) * (1 - (1 - G)^G)^(beta-1) * (G/(1-G) - log(1-G)) 
  return(f)
}

p=plot(c(-2.7,2.7), c(0,0.8), type = "n", xlab="x", ylab="hrf")
lambda = c(1.9,1.6,2.0,0.9,1.2)
beta   = c(0.1,0.1,0.7,0.1,1.6)
mu     = c(0.6,0.8,2.8,0.4,0.1)
sigma  = c(0.7,0.7,1.0,0.7,1.5)
cor1=c('black','darkred','darkblue','darkgreen',"darkorange")
linhas=c(1:5)
for( i in 1:length(lambda)){
  p[i]=lines(t,hMKNFN(lambda[i],beta[i],mu[i],sigma[i]),lwd=c(2,2,2,2,2),lty=linhas[i],col=cor1[i])}
legend("topleft", expression(
  paste(lambda,"=1.9; ",beta,"=0.1; ", mu,"=0.6; ",sigma,"=0.7"),    
  paste(lambda,"=1.6; ",beta,"=0.1; ", mu,"=0.8; ",sigma,"=0.7"),
  paste(lambda,"=2.0; ",beta,"=0.7; ", mu,"=2.8; ",sigma,"=1.0"),
  paste(lambda,"=0.9; ",beta,"=0.1; ", mu,"=0.4; ",sigma,"=0.7"),
  paste(lambda,"=1.2; ",beta,"=1.6; ", mu,"=0.1; ",sigma,"=1.5")),
  lty=c(1:5), lwd=c(2,2,2,2,2),col=cor1,
  bty="n")