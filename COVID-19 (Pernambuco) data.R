
setwd("D:\\Sync\\Sync\\Doutorado UFPE\\Tese\\tese_escrita_inglês\\tese_MO-NFGF\\Aplicação_MONFBXII\\dados")

require("csv")
require("xlsx")
library(utils)
require(survival)
require(AdequacyModel)
require(zipfR)

dados = read.csv("CVPE3.csv", 1)
t = dados$tempo1
hist(t)
ekm=survfit(Surv(t)~1)
descriptive(t)
########################################################################################
########################################MKNFW###########################################
########################################################################################
pdf_MKNFW = function(par,x){
  lambda = par[1]
  beta   = par[2]
  alpha  = par[3]
  tau    = par[4]
  g = alpha*tau^(-alpha)*x^(alpha-1)*exp(-(x/tau)^alpha)
  G = 1 - exp(-(x/tau)^alpha)
  f = lambda*beta*g*(1-G)^(-beta*G)*(1 - (1 - G)^G)^(beta-1)*(G/(1-G) - log(1-G))*
    exp(-lambda*((1 - G)^(-G) - 1)^beta)
}

# EPI acumulada

cdf_MKNFW = function(par,x){
  lambda = par[1]
  beta   = par[2]
  alpha  = par[3]
  tau    = par[4]
  G = 1 - exp(-(x/tau)^alpha)
  1 - exp(-lambda*((1 - G)^(-G) - 1)^beta)
}

# EPI sobrevivencia

sdf_MKNFW = function(par,x){
  lambda = par[1]
  beta   = par[2]
  alpha  = par[3]
  tau    = par[4]
  G = 1 - exp(-(x/tau)^alpha)
  exp(-lambda*((1 - G)^(-G) - 1)^beta)
}
results_MKNFW = goodness.fit(pdf = pdf_MKNFW, cdf = cdf_MKNFW, 
                             starts = c(1.553910, 0.164072,  3.834597, 46.604641), data = t, 
                             method = "B", domain = c(0, 100000),mle = NULL, lim_inf = c(0,0,0,0),
                             lim_sup = c(20,20,20,20), S = 250, prop=0.1, N=500)

results_MKNFW

########################################################################################
###############################KW_W#####################################################
########################################################################################
pdf_KWW <- function(par,x){
  a     = par[1]
  b     = par[2]
  alpha = par[3]
  phi   = par[4]
  g = alpha*phi^(-alpha)*x^(alpha - 1)*exp(-(x/phi)^alpha)
  G = 1 - exp(-(x/phi)^alpha)
  a*b*g*G^(a-1)*(1 - G^a)^(b-1)
  
}

# KWW acumulada

cdf_KWW <- function(par,x){
  a     = par[1]
  b     = par[2]
  alpha = par[3]
  phi   = par[4]
  G = 1 - exp(-(x/phi)^alpha)
  1 - (1 - G^a)^b
}

# KWW sobrevivencia

sdf_KWW <- function(par,x){
  a     = par[1]
  b     = par[2]
  alpha = par[3]
  phi   = par[4]
  G     = 1 - exp(-(x/phi)^alpha)
  (1 - G^a)^b
}
results_KWW <- goodness.fit(pdf = pdf_KWW, cdf = cdf_KWW, 
                            starts = c(0.1486102,  1.9614828,  7.2077348, 81.4755704), data = t, 
                            method = "B", mle = NULL, domain = c(0, Inf),lim_inf = c(0,0,0,0,0),
                            lim_sup = c(Inf,Inf,Inf,Inf), S = 250, prop=0.1, N=50)

results_KWW
########################################################################################
#################################BW#####################################################
########################################################################################
pdf_BW <- function(par,x){
  a     = par[1]
  b     = par[2]
  alpha = par[3]
  phi   = par[4]
  g = alpha*phi^(-alpha)*x^(alpha - 1)*exp(-(x/phi)^alpha)
  G = 1 - exp(-(x/phi)^alpha)
  BWpdf = (1/(beta(a,b)))*g*G^(a-1)*(1 - G)^(b-1)
  return(BWpdf)
  
}

# KWW acumulada

cdf_BW <- function(par,x){
  a     = par[1]
  b     = par[2]
  alpha = par[3]
  phi   = par[4]
  G = 1 - exp(-(x/phi)^alpha)
  MCWcdf = pbeta(G,a,b)
  return(MCWcdf)
}

# KWW sobrevivencia

sdf_BW <- function(par,x){
  a     = par[1]
  b     = par[2]
  alpha = par[3]
  phi   = par[4]
  G = 1 - exp(-(x/phi)^alpha)
  MCWcdf = pbeta(G,a,b)
  return(1 - MCWcdf)
}
results_BW <- goodness.fit(pdf = pdf_BW, cdf = cdf_BW, 
                           starts = c(0.3845483,  0.7184505,  2.6329727, 43.1315192), data = t, 
                           method = "B", mle = NULL, domain = c(0, Inf),lim_inf = c(0,0,0,0,0),
                           lim_sup = c(2,2,2,2,2), S = 250, prop=0.1, N=50)

results_BW
########################################################################################
########################################WW##############################################
########################################################################################
pdf_WW <- function(par,x){
  a     = 1
  lambda= par[1]
  beta  = par[2]
  alpha = par[3]
  phi   = par[4]
  g = alpha*phi^(-alpha)*x^(alpha - 1)*exp(-(x/phi)^alpha)
  G = 1 - exp(-(x/phi)^alpha)
  f = a*lambda*beta*g*G^(a*beta-1)*(1-G)^(-(a*beta+1))*exp(-lambda*(G/(1-G))^(a*beta))
  return(f)
  
}

# bmw acumulada

cdf_WW <- function(par,x){
  a     = 1
  lambda= par[1]
  beta  = par[2]
  alpha = par[3]
  phi   = par[4]
  G = 1 - exp(-(x/phi)^alpha)
  cdf = 1 - exp(-lambda*(G/(1-G))^(a*beta))
  return(cdf)
}

# bmw sobrevivencia

sdf_WW <- function(par,x){
  a     = 1
  lambda= par[1]
  beta  = par[2]
  alpha = par[3]
  phi   = par[4]
  G = 1 - exp(-(x/phi)^alpha)
  sf =  exp(-lambda*(G/(1-G))^(a*beta))
  return(sf)
}
results_WW <- goodness.fit(pdf = pdf_WW, cdf = cdf_WW, 
                           starts = c(1,.01,.1,1), data = t, 
                           method = "B", mle = NULL, domain = c(0, Inf),lim_inf = c(0,0,0,0,0),
                           lim_sup = c(100,100,100,100,100), S = 250, prop=0.1, N=50)

results_WW
#####################################################################################################
#################################=====G_We=====######################################################
#####################################################################################################
pdf_GWe = function(par,x){
  a     = par[1]
  alpha = par[2]
  beta  = par[3]
  g = dweibull(x, shape = alpha, scale = beta)
  G = pweibull(x, shape = alpha, scale = beta)
  (g/gamma(a))*(-log(1 - G))^(a-1)
}

cdf_GWe = function(par, x){
  a     = par[1]
  alpha = par[2]
  beta  = par[3]
  G = pweibull(x, shape = alpha, scale = beta)
  (Igamma(a,-log(1-G)))/gamma(a)
}

sdf_GWe = function(par, x){
  a     = par[1]
  alpha = par[2]
  beta  = par[3]
  G = pweibull(x, shape = alpha, scale = beta)
  1 - (Igamma(a,-log(1-G)))/gamma(a)
}

results_GWe <- goodness.fit(pdf = pdf_GWe, cdf = cdf_GWe, 
                            starts = c(1,.1,1), data = t, 
                            method = "B", domain = c(0, Inf), lim_inf = c(0,0,0),
                            lim_sup = c(2,2,2), mle = NULL, S = 250, prop=0.1,N=50)

results_GWe
################################################################################################
###############################distribuição weibull#############################################
################################################################################################

pdf_weibull <- function(par,x){
  alpha = par[1]
  lambda = par[2]
  lambda*(alpha^lambda)*(x^(lambda-1))*exp(-(alpha*x)^lambda)
}
cdf_weibull <- function(par,x){
  alpha = par[1]
  lambda = par[2]
  1 - exp(-(alpha*x)^lambda)
}

sdf_weibull = function(par,x){
  alpha = par[1]
  lambda = par[2]
  G = 1 - exp(-(alpha*x)^lambda)
  1 - (G)
}
results_WE = goodness.fit(pdf = pdf_weibull, cdf = cdf_weibull,
                          starts = c(0.03162729, 1.38327486), data = t,
                          method = "BFGS", domain = c(0,Inf), mle = NULL, lim_inf = c(0,0), lim_sup = c(10,10),
                          N = 100, S = 250)

results_WE
########################################################################################
########################################LW##############################################
########################################################################################
pdf_LW <- function(par,x){
  lambda= par[1]
  beta  = par[2]
  alpha = par[3]
  phi   = par[4]
  g = alpha*phi^(-alpha)*x^(alpha - 1)*exp(-(x/phi)^alpha)
  G = 1 - exp(-(x/phi)^alpha)
  f = ((lambda*beta^lambda*g)/((1 - G)^2))*(beta + (G)/(1-G))^(-(lambda+1)) 
  return(f)
  
}

# bmw acumulada

cdf_LW <- function(par,x){
  lambda= par[1]
  beta  = par[2]
  alpha = par[3]
  phi   = par[4]
  G     = 1 - exp(-(x/phi)^alpha)
  cdf   = 1 - beta^lambda*(beta + G/(1 - G))^(-lambda)
  return(cdf)
}

# bmw sobrevivencia

sdf_LW <- function(par,x){
  a     = 1
  lambda= par[1]
  beta  = par[2]
  alpha = par[3]
  phi   = par[4]
  G     = 1 - exp(-(x/phi)^alpha)
  sf    =  beta^lambda*(beta + G/(1 - G))^(-lambda)
  return(sf)
}
results_LW <- goodness.fit(pdf = pdf_LW, cdf = cdf_LW, 
                           starts = c(1,.01,.1,1), data = t, 
                           method = "B", mle = NULL, domain = c(0, Inf),lim_inf = c(0,0,0,0),
                           lim_sup = c(100,100,100,100), S = 250, prop=0.1, N=50)

results_LW
################################################################################
################################################################################
################################################################################
##grafico do histograma da distribuição

x = seq(0.7,80, by = 0.01)
hist(t, probability = TRUE, xlab = "x", ylab = "pdf", main = "",border = "black",col = "white")
lines(x, pdf_MKNFW(par = results_MKNFW$mle, x), col = 'darkred',lwd = 2)
lines(x, pdf_KWW(par = results_KWW$mle, x), col = 'darkblue',lwd = 2, lty = 2)
#lines(x, pdf_BW(par = results_BW$mle, x), col = "darkgreen", lwd = 2, lty = 4)
#lines(x, pdf_WW(par = results_WW$mle, x),col = "darkgreen", lwd = 2, lty = 1)
lines(x, pdf_GWe(par = results_GWe$mle, x),col = "darkgreen", lwd = 2, lty = 4)
#lines(x, pdf_LW(par = results_LW$mle, x),col = "blue", lwd = 2, lty = 1)
#lines(x, pdf_weibull(par = results_WE$mle, x),col = "green", lwd = 2, lty = 1)
#lines(x, pdf_KwWeibull(par = results_KW_W$mle, x), col = "green")
#lines(x, pdf_BOLLW(par = results_BOLLW$mle, x), col = 'purple', lwd = 1)
legend("topright",lwd=c(2),lty=c(1,2,4),
       c("MKFW","KwW","GW"),
       col=c("darkred","darkblue","darkgreen"),
       bty="n")

##gr?fico de fun??o de sobreviv?ncia

#x = seq(0,80, by = 0.01)
#plot(ekm,conf.int=F,xlab = "x" , ylab = "Survival function",main = "")
#lines(x, sdf_MKNFW(par = results_MKNFW$mle, x), col = 'darkred',lwd = 2,lty = 1)
#lines(x, sdf_WW(par = results_WW$mle, x), col = 'darkblue',lwd = 2, lty = 2)
#lines(x, sdf_BMW(par = results_BMW$mle, x),col = "darkgreen", lwd = 2, lty = 4)
#lines(x, sdf_MCW(par = results_MCW$mle, x),col = "green", lwd = 2, lty = 1)
#lines(x, sdf_GWe(par = results_GWe$mle, x),col = "darkgreen", lwd = 2, lty = 4)
#lines(x, sdf_EWLL(par = results_EWLL$mle, x), col = 1,lwd = 1, lty = 4)
#lines(x, sdf_MW(par = results_MW$mle, x), col = "purple")
#lines(x, sdf_BOLLW(par = results_BOLLW$mle, x), col = 'purple',lwd = 2,lty = 1)
#lines(x, sdf_KWW(par = results_KWW$mle, x), col = "darkblue",lwd=2,lty=2)
#legend("topright",lwd=c(2),lty=c(1,2,4),
#       c("MKNFW","KW","GW"),
#       col=c("darkred","darkblue","darkgreen"),
#       bty="n")


##grafico da função empírica
plot(ecdf(t), lty=1, lwd=1, do.points=FALSE, verticals=T, 
     ylim=c(0.0,1.0), xlim = c(0.3,76), ylab="cdf", main="", col.01line="white")
x = seq(0,90, by = 0.01)
#plot(ekm$time, 1-ekm$surv, xlab = "x", ylab = "cdf",type = "l", lwd=1, col=1, lty=1)
lines(x, cdf_MKNFW(par = results_MKNFW$mle, x), col = 'darkred',lwd = 2)
lines(x, cdf_KWW(par = results_KWW$mle, x), col = 'darkblue',lwd = 2, lty = 2)
lines(x, cdf_GWe(par = results_GWe$mle, x),col = "darkgreen", lwd = 2, lty = 4)
#lines(x, cdf_MOFW(par = results_MOFW$mle, x),col = "green", lwd = 1, lty = 1)
#lines(x, cdf_FW(par = results_FW$mle, x),col = "blue", lwd = 1, lty = 1)
#lines(x, cdf_KwBXII(par = results_KwBXII$mle, x),col = "green", lwd = 1, lty = 1)
#lines(x, cdf_BFW(par = results_BFW$mle, x), col = 1,lwd = 1, lty = 4)
#lines(x, cdf_MW(par = results_MW$mle, x), col = "purple")
#lines(x, cdf_BOLLW(par = results_BOLLW$mle, x), col = 'purple', lwd = 2,lty = 1)
#lines(x, cdf_KwWeibull(par = results_KW_W$mle, x), col = "green")
#legend("bottomright",lwd=c(2),lty=c(1,2,4),
#       c("MKFW","KwW","GW"),
#       col=c("darkred","darkblue","darkgreen"),
#       bty="n")


## Comparando
MKNFW_estimativas = c(results_MKNFW$mle)
MKNFW_erro = c(results_MKNFW$Erro)
KW_estimativa = results_KWW$mle
KW_erro = results_KWW$Erro
BW_estimativa = results_BW$mle
BW_erro = results_BW$Erro
WW_estimativas = c(results_WW$mle)
WW_erro = c(results_WW$Erro)
GW_estimativas = c(results_GWe$mle)
GW_erro = c(results_GWe$Erro)
LW_estimativas = c(results_LW$mle)
LW_erro = c(results_LW$Erro)
WE_estimativas = c(results_WE$mle)
WE_erro = c(results_WE$Erro)

tabela1 =  rbind(MKNFW_estimativas,MKNFW_erro,KW_estimativa, KW_erro ,BW_estimativa, BW_erro, WW_estimativas, WW_erro,
                 LW_estimativas,LW_erro,GW_estimativas,GW_erro,WE_estimativas,WE_erro)                  
colnames(tabela1) = c("a","b","alpha","beta","gama")
stargazer::stargazer(round(tabela1,6))
tabela1

#####################################################################################################################################
#####################################################################################################################################
A.est <- c(results_MKNFW$A,results_KWW$A, results_BW$A, results_WW$A, results_LW$A, results_GWe$A,results_WE$A)
W.est <- c(results_MKNFW$W,results_KWW$W, results_BW$W, results_WW$W, results_LW$W, results_GWe$W,results_WE$W)
CAIC.est <- c(results_MKNFW$`CAIC `,results_KWW$`CAIC `, results_BW$`CAIC `, results_WW$`CAIC `, results_LW$`CAIC `,
              results_GWe$`CAIC `,results_WE$`CAIC `) 
AIC.est <- c(results_MKNFW$AIC,results_KWW$AIC, results_BW$AIC, results_WW$AIC, results_LW$AIC,
             results_GWe$AIC,results_WE$AIC) 
BIC.est <- c(results_MKNFW$BIC,results_KWW$BIC, results_BW$BIC, results_WW$BIC, results_LW$BIC,
             results_GWe$BIC,results_WE$BIC)
HQIC.est <- c(results_MKNFW$HQIC,results_KWW$HQIC, results_BW$HQIC, results_WW$HQIC, results_LW$HQIC,
              results_GWe$HQIC,results_WE$HQIC)
#KS.est <- c(results_MKOLLW$KS,results_KWOLLW$KS, results_MCW$KS, results_BMW$KS, results_KMB$KS,
# results_KWW$KS,  results_BW$KS)
tabela <- cbind(W.est, A.est, AIC.est, CAIC.est, BIC.est, HQIC.est)
rownames(tabela) <- c("MKNFW","KW", "BW", "WW","LW","GW","WE")
tabela
stargazer::stargazer(round(tabela,6))


results_WE$KS
##############################################################################################################
#teste da razao de verossimilhan?as generalizado

#MONFBXII vs NFBXII
RL1   = (1/sqrt(length(t)))*sum((log(pdf_MKNFW(par = results_MKNFW$mle, t)/pdf_KWW(par = results_KWW$mle, t))))
W1    = (1/length(t))*sum((log(pdf_MKNFW(par = results_MKNFW$mle, t)/pdf_KWW(par = results_KWW$mle, t)))^2)
S1    = ((1/length(t))*sum((log(pdf_MKNFW(par = results_MKNFW$mle, t)/pdf_KWW(par = results_KWW$mle, t)))))^2
RLNN1 = (RL1)/(W1-S1)
RLNN1

#MKOLLW vs MCW
RL2   = (1/sqrt(length(t)))*sum((log(pdf_MKNFW(par = results_MKNFW$mle, t)/pdf_BW(par = results_BW$mle, t))))
W2    = (1/length(t))*sum((log(pdf_MKNFW(par = results_MKNFW$mle, t)/pdf_BW(par = results_BW$mle, t)))^2)
S2    = ((1/length(t))*sum((log(pdf_MKNFW(par = results_MKNFW$mle, t)/pdf_BW(par = results_BW$mle, t)))))^2
RLNN2 = (RL2)/(W2-S2)
RLNN2

#MONFBXII vs BMW
RL3   = (1/sqrt(length(t)))*sum((log(pdf_MKNFW(par = results_MKNFW$mle, t)/pdf_WW(par = results_WW$mle, t))))
W3    = (1/length(t))*sum((log(pdf_MKNFW(par = results_MKNFW$mle, t)/pdf_WW(par = results_WW$mle, t)))^2)
S3    = ((1/length(t))*sum((log(pdf_MKNFW(par = results_MKNFW$mle, t)/pdf_WW(par = results_WW$mle, t)))))^2
RLNN3 = (RL3)/(W3-S3)
RLNN3

#MKOLLW vs KMW
RL4   = (1/sqrt(length(t)))*sum((log(pdf_MKNFW(par = results_MKNFW$mle, t)/pdf_LW(par = results_LW$mle, t))))
W4    = (1/length(t))*sum((log(pdf_MKNFW(par = results_MKNFW$mle, t)/pdf_LW(par = results_LW$mle, t)))^2)
S4    = ((1/length(t))*sum((log(pdf_MKNFW(par = results_MKNFW$mle, t)/pdf_LW(par = results_LW$mle, t)))))^2
RLNN4 = (RL4)/(W4-S4)
RLNN4

#MKOLLW vs KMW
RL5   = (1/sqrt(length(t)))*sum((log(pdf_MKNFW(par = results_MKNFW$mle, t)/pdf_GWe(par = results_GWe$mle, t))))
W5    = (1/length(t))*sum((log(pdf_MKNFW(par = results_MKNFW$mle, t)/pdf_GWe(par = results_GWe$mle, t)))^2)
S5    = ((1/length(t))*sum((log(pdf_MKNFW(par = results_MKNFW$mle, t)/pdf_GWe(par = results_GWe$mle, t)))))^2
RLNN5 = (RL5)/(W5-S5)
RLNN5

#MKOLLW vs KMW
RL6   = (1/sqrt(length(t)))*sum((log(pdf_MKNFW(par = results_MKNFW$mle, t)/pdf_weibull(par = results_WE$mle, t))))
W6    = (1/length(t))*sum((log(pdf_MKNFW(par = results_MKNFW$mle, t)/pdf_weibull(par = results_WE$mle, t)))^2)
S6    = ((1/length(t))*sum((log(pdf_MKNFW(par = results_MKNFW$mle, t)/pdf_weibull(par = results_WE$mle, t)))))^2
RLNN6 = (RL6)/(W6-S6)
RLNN6

