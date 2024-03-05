
require("csv")
require("xlsx")
library(utils)
require(survival)
require(AdequacyModel)
require(zipfR)

#book weibull models page 231 data set 12.2 \ better
t = c(0.13, 0.62, 0.75, 0.87, 1.56, 2.28, 3.15, 3.25, 3.55, 4.49,
      4.50, 4.61, 4.79, 7.17, 7.31, 7.43, 7.84, 8.49, 8.94, 9.40,
      9.61, 9.84, 10.58, 11.18, 11.84, 13.28, 14.47, 14.79, 15.54, 16.90,
      17.25, 17.37, 18.69, 18.78, 19.88, 20.06, 20.10, 20.95, 21.72, 23.87)

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
                             starts = c(0.9080219,  0.1420226,  3.6867290, 12.1536463), data = t, 
                             method = "B", domain = c(0, Inf),mle = NULL, lim_inf = c(0,0,0,0),
                             lim_sup = c(20,20,20,20), S = 250, prop=0.1, N=500)

results_MKNFW

########################################################################################
################################KFW#####################################################
########################################################################################
pdf_KFW <- function(par,x){
  a     = par[1]
  b     = par[2]
  alpha = par[3]
  beta  = par[4]
  g = (alpha + beta/x^2)*(exp(alpha*x - beta/x))*(exp(-exp(alpha*x - beta/x)))
  G = 1 - exp(-exp(alpha*x - beta/x))
  a*b*g*G^(a-1)*(1 - G^a)^(b-1)
  
}

# KWW acumulada

cdf_KFW <- function(par,x){
  a     = par[1]
  b     = par[2]
  alpha = par[3]
  beta  = par[4]
  G = 1 - exp(-exp(alpha*x - beta/x))
  1 - (1 - G^a)^b
}

# KWW sobrevivencia

sdf_KFW <- function(par,x){
  a     = par[1]
  b     = par[2]
  alpha = par[3]
  beta  = par[4]
  G     = 1 - exp(-exp(alpha*x - beta/x))
  (1 - G^a)^b
}
results_KFW <- goodness.fit(pdf = pdf_KFW, cdf = cdf_KFW, 
                            starts = c(0.2393011, 0.1197107, 0.1507195, 1.8208803), data = t, 
                            method = "B", mle = NULL, domain = c(0, Inf),lim_inf = c(0,0,0,0),
                            lim_sup = c(2,2,2,2), S = 250, prop=0.1, N=50)

results_KFW

#####################################################################################################
#########################################GOLLFW######################################################
#####################################################################################################
pdf_GOLLW = function(par,x){
  a     = par[1]
  b     = par[2]
  alpha = par[3]
  tau   = par[4]
  g = alpha*tau^(-alpha)*x^(alpha-1)*exp(-(x/tau)^alpha)
  G = 1 - exp(-(x/tau)^alpha)
  (a*b*g*G^(a*b-1)*(1 - G^b)^(a-1))/(G^(a*b) + (1 - G^b)^a)^2
}

cdf_GOLLW = function(par, x){
  a     = par[1]
  b     = par[2]
  alpha = par[3]
  tau   = par[4]
  G = 1 - exp(-(x/tau)^alpha)
  G^(a*b)/(G^(a*b) + (1 - G^b)^a)
}

sdf_GOLLW = function(par, x){
  a     = par[1]
  b     = par[2]
  alpha = par[3]
  tau   = par[4]
  G = 1 - exp(-(x/tau)^alpha)
  1 - G^(a*b)/(G^(a*b) + (1 - G^b)^a)
}

results_GOLLW <- goodness.fit(pdf = pdf_GOLLW, cdf = cdf_GOLLW, 
                              starts = c(1,1,.01,.1), data = t, 
                              method = "B", domain = c(0, Inf), lim_inf = c(0,0,0,0),
                              lim_sup = c(2,2,2,2), mle = NULL, S = 250, prop=0.1,N=50)

results_GOLLW


########################################################################################
########################################WMOW############################################
########################################################################################
pdf_WMOW = function(par,x){
  lambda = par[1]
  beta   = par[2]
  alpha  = par[3]
  tau    = par[4]
  g = alpha*tau^(-alpha)*x^(alpha-1)*exp(-(x/tau)^alpha)
  G = 1 - exp(-(x/tau)^alpha)
  H = (lambda*(1 - G))
  H1 = 1 - (1 - lambda)*(1 -G)
  f = (beta*g)/((1-G)*H1)*(-log(H/H1))^(beta-1)*exp(-(-log(H/H1))^beta)
}

# EPI acumulada

cdf_WMOW = function(par,x){
  lambda = par[1]
  beta   = par[2]
  alpha  = par[3]
  tau    = par[4]
  G = 1 - exp(-(x/tau)^alpha)
  H = (lambda*(1 - G))
  H1 = 1 - (1 - lambda)*(1 -G)
  1 - exp(-(-log(H/H1))^beta)
}

# EPI sobrevivencia

sdf_WMOW = function(par,x){
  lambda = par[1]
  beta   = par[2]
  alpha  = par[3]
  tau    = par[4]
  G = 1 - exp(-(x/tau)^alpha)
  H = (lambda*(1 - G))
  H1 = 1 - (1 - lambda)*(1 -G)
  exp(-(-log(H/H1))^beta)
}
results_WMOW = goodness.fit(pdf = pdf_WMOW, cdf = cdf_WMOW, 
                            starts = c(1,.1,1,1), data = t, 
                            method = "B", domain = c(0, Inf),mle = NULL, lim_inf = c(0,0,0,0),
                            lim_sup = c(20,20,20,20), S = 250, prop=0.1, N=500)

results_WMOW

########################################################################################
######################Extend WLL########################################################
########################################################################################

pdf_EWLL <- function(par,x){
  alpha = par[1]
  theta = par[2]
  beta  = par[3]
  g = beta*theta*alpha*x^(alpha - 1)*(1 + x^alpha)^(theta - 1)*
    ((1 + x^alpha)^theta - 1)^(beta - 1)*exp(-((1 + x^alpha)^theta - 1)^beta) 
}

# bw acumulada

cdf_EWLL <- function(par,x){
  alpha = par[1]
  theta = par[2]
  beta  = par[3]
  G = 1 - exp(-((1 + x^alpha)^theta - 1)^beta)
}

# bw sobrevivencia

sdf_EWLL <- function(par,x){
  alpha = par[1]
  theta = par[2]
  beta  = par[3]
  G = 1 - exp(-((1 + x^alpha)^theta - 1)^beta)
  1-G
}
results_EWLL <- goodness.fit(pdf = pdf_EWLL, cdf = cdf_EWLL, 
                             starts = c(1,1,1), data = t, 
                             method = "S", mle = NULL, domain = c(0, Inf),lim_inf = c(0,0,0,0),
                             lim_sup = c(2,2,2,2), S = 250, prop=0.1, N=50)

results_EWLL
################################################################################
################################################################################
################################################################################
##grafico do histograma da distribuição

x = seq(0,25, by = 0.01)
hist(t, probability = TRUE, xlab = "x", ylab = "pdf", main = "",border = "black",col = "white",ylim= c(0,0.1))
lines(x, pdf_MKNFW(par = results_MKNFW$mle, x), col = 'darkred',lwd = 2)
#lines(x, pdf_KWW(par = results_KWW$mle, x), col = 'darkblue',lwd = 2, lty = 2)
#lines(x, pdf_WMOW(par = results_WMOW$mle, x), col = "darkgreen", lwd = 2, lty = 4)
#lines(x, pdf_WW(par = results_WW$mle, x),col = "darkgreen", lwd = 2, lty = 1)
lines(x, pdf_KFW(par = results_KFW$mle, x),col = "darkgreen", lwd = 2, lty = 4)
lines(x, pdf_GOLLW(par = results_GOLLW$mle, x),col = "darkblue", lwd = 2, lty = 2)
#lines(x, pdf_weibull(par = results_WE$mle, x),col = "green", lwd = 2, lty = 1)
#lines(x, pdf_KwWeibull(par = results_KW_W$mle, x), col = "green")
#lines(x, pdf_BOLLW(par = results_BOLLW$mle, x), col = 'purple', lwd = 1)
legend("topright",lwd=c(2),lty=c(1,2,4),
       c("MKFW","GOLLW","KFW"),
       col=c("darkred","darkblue","darkgreen"),
       bty="n")

##gr?fico de fun??o de sobreviv?ncia

#x = seq(0,25, by = 0.01)
#plot(ekm,conf.int=F,xlab = "x" , ylab = "Survival function",main = "")
#lines(x, sdf_MKNFW(par = results_MKNFW$mle, x), col = 'darkred',lwd = 2,lty = 1)
#lines(x, sdf_KWW(par = results_KWW$mle, x), col = 'darkblue',lwd = 2, lty = 2)
#lines(x, sdf_BMW(par = results_BMW$mle, x),col = "darkgreen", lwd = 2, lty = 4)
#lines(x, sdf_WW(par = results_WW$mle, x),col = "green", lwd = 2, lty = 1)
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
     ylim=c(0.0,1.0), xlim = c(0.5,23.5), ylab="cdf", main="", col.01line="white")
x = seq(0,25, by = 0.01)
#plot(ekm$time, 1-ekm$surv, xlab = "x", ylab = "cdf",type = "l", lwd=1, col=1, lty=1)
lines(x, cdf_MKNFW(par = results_MKNFW$mle, x), col = 'darkred',lwd = 2,lty=1)
lines(x, cdf_GOLLW(par = results_GOLLW$mle, x),col = "darkblue", lwd = 2, lty = 2)
#lines(x, cdf_MOFW(par = results_MOFW$mle, x),col = "green", lwd = 1, lty = 1)
lines(x, cdf_KFW(par = results_KFW$mle, x),col = "darkgreen", lwd = 2, lty = 4)
#lines(x, cdf_KwBXII(par = results_KwBXII$mle, x),col = "green", lwd = 1, lty = 1)
#lines(x, cdf_BFW(par = results_BFW$mle, x), col = 1,lwd = 1, lty = 4)
#lines(x, cdf_MW(par = results_MW$mle, x), col = "purple")
#lines(x, cdf_BOLLW(par = results_BOLLW$mle, x), col = 'purple', lwd = 2,lty = 1)
#lines(x, cdf_KwWeibull(par = results_KW_W$mle, x), col = "green")
legend("bottomright",lwd=c(2),lty=c(1,2,4),
       c("MKFW","GOLLW","KFW"),
       col=c("darkred","darkblue","darkgreen"),
       bty="n")


## Comparando
MKNFW_estimativas = c(results_MKNFW$mle)
MKNFW_erro = c(results_MKNFW$Erro)
KFW_estimativa = results_KFW$mle
KFW_erro = results_KFW$Erro
GOLLW_estimativa = results_GOLLW$mle
GOLLW_erro = results_GOLLW$Erro
WMOW_estimativas = c(results_WMOW$mle)
WMOW_erro = c(results_WMOW$Erro)
EWLL_estimativas = c(results_EWLL$mle)
EWLL_erro = c(results_EWLL$Erro)


tabela1 =  rbind(MKNFW_estimativas,MKNFW_erro,KFW_estimativa, KFW_erro ,GOLLW_estimativa, GOLLW_erro, WMOW_estimativas, WMOW_erro,
                 EWLL_estimativas,EWLL_erro)                  
colnames(tabela1) = c("a","b","alpha","beta","gama")
stargazer::stargazer(round(tabela1,6))
tabela1

#####################################################################################################################################
#####################################################################################################################################
A.est <- c(results_MKNFW$A,results_KFW$A, results_GOLLW$A, results_WMOW$A, results_EWLL$A)
W.est <- c(results_MKNFW$W,results_KFW$W, results_GOLLW$W, results_WMOW$W, results_EWLL$W)
CAIC.est <- c(results_MKNFW$`CAIC `,results_KFW$`CAIC `, results_GOLLW$`CAIC `, results_WMOW$`CAIC `, results_EWLL$`CAIC `) 
AIC.est <- c(results_MKNFW$AIC,results_KFW$AIC, results_GOLLW$AIC, results_WMOW$AIC, results_EWLL$AIC) 
BIC.est <- c(results_MKNFW$BIC,results_KFW$BIC, results_GOLLW$BIC, results_WMOW$BIC, results_EWLL$BIC)
HQIC.est <- c(results_MKNFW$HQIC,results_KFW$HQIC, results_GOLLW$HQIC, results_WMOW$HQIC, results_EWLL$HQIC)
#KS.est <- c(results_MKOLLW$KS,results_KWOLLW$KS, results_MCW$KS, results_BMW$KS, results_KMB$KS,
# results_KWW$KS,  results_BW$KS)
tabela <- cbind(W.est, A.est, AIC.est, CAIC.est, BIC.est, HQIC.est)
rownames(tabela) <- c("MKNFW","KFW", "GOLLW", "WMOW","EWLL")
tabela
stargazer::stargazer(round(tabela,6))


results_WE$KS
##############################################################################################################
#teste da razao de verossimilhan?as generalizado

#MONFBXII vs NFBXII
RL1   = (1/sqrt(length(t)))*sum((log(pdf_MKNFW(par = results_MKNFW$mle, t)/pdf_KFW(par = results_KFW$mle, t))))
W1    = (1/length(t))*sum((log(pdf_MKNFW(par = results_MKNFW$mle, t)/pdf_KFW(par = results_KFW$mle, t)))^2)
S1    = ((1/length(t))*sum((log(pdf_MKNFW(par = results_MKNFW$mle, t)/pdf_KFW(par = results_KFW$mle, t)))))^2
RLNN1 = (RL1)/(W1-S1)
RLNN1

#MKOLLW vs MCW
RL2   = (1/sqrt(length(t)))*sum((log(pdf_MKNFW(par = results_MKNFW$mle, t)/pdf_GOLLW(par = results_GOLLW$mle, t))))
W2    = (1/length(t))*sum((log(pdf_MKNFW(par = results_MKNFW$mle, t)/pdf_GOLLW(par = results_GOLLW$mle, t)))^2)
S2    = ((1/length(t))*sum((log(pdf_MKNFW(par = results_MKNFW$mle, t)/pdf_GOLLW(par = results_GOLLW$mle, t)))))^2
RLNN2 = (RL2)/(W2-S2)
RLNN2

#MONFBXII vs BMW
RL3   = (1/sqrt(length(t)))*sum((log(pdf_MKNFW(par = results_MKNFW$mle, t)/pdf_WMOW(par = results_WMOW$mle, t))))
W3    = (1/length(t))*sum((log(pdf_MKNFW(par = results_MKNFW$mle, t)/pdf_WMOW(par = results_WMOW$mle, t)))^2)
S3    = ((1/length(t))*sum((log(pdf_MKNFW(par = results_MKNFW$mle, t)/pdf_WMOW(par = results_WMOW$mle, t)))))^2
RLNN3 = (RL3)/(W3-S3)
RLNN3

#MKOLLW vs KMW
RL4   = (1/sqrt(length(t)))*sum((log(pdf_MKNFW(par = results_MKNFW$mle, t)/pdf_EWLL(par = results_EWLL$mle, t))))
W4    = (1/length(t))*sum((log(pdf_MKNFW(par = results_MKNFW$mle, t)/pdf_EWLL(par = results_EWLL$mle, t)))^2)
S4    = ((1/length(t))*sum((log(pdf_MKNFW(par = results_MKNFW$mle, t)/pdf_EWLL(par = results_EWLL$mle, t)))))^2
RLNN4 = (RL4)/(W4-S4)
RLNN4

