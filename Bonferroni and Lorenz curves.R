#########curvas de Lorenz##########

rm(list = ls())

require(GoFKernel)

rMKNFW = function (n = 1, f, lower = -Inf, upper = Inf, kind = "density") 
{
  if (!is.finite(f(lower))) 
    lower <- lower + .Machine$double.eps
  if (!is.finite(f(upper))) 
    upper <- upper - .Machine$double.eps
  if (!is.numeric(lower) || !is.numeric(upper) || lower >= 
      upper) 
    stop("lower < upper is not fulfilled")
  cumulative <- function(f, lower) {
    function(z) integrate(f, lower, upper = z, rel.tol = 1e-10, 
                          subdivisions = 1000)$value
  }
  if (kind == "density") {
    f.distr <- cumulative(f = f, lower = lower)
    limits <- support.facto(f, lower = lower, upper = upper)
    lower <- max(c(lower, limits[1]))
    upper <- min(c(upper, limits[2]))
  }
  else {
    f.distr <- f
  }
  cum.inverse <- inverse(f.distr, lower = lower, upper = upper)
  sample.unif <- n
  sample.function <- sapply(sample.unif, cum.inverse)
  return(sample.function)
}

f1 = function(x){
  G <- 1 - exp(-(x/phi)^alpha)
  H = 1 - (1 - G)^G
  cdf = 1 - exp(-lambda*(H/(1-H))^beta)
  return(cdf)
}

f <- function(x){
  g <- (alpha/phi) * (x/phi)^(alpha-1) * exp(-(x/phi)^alpha)
  G <- 1 - exp(-(x/phi)^alpha)
  h = g*(1-G)^G*(G/(1-G) - log(1-G))
  H = 1 - (1 - G)^G
  pdf = lambda*beta*h*(H)^(beta-1)*(1-H)^(-(beta+1))*exp(-lambda*(H/(1-H))^beta)
  return(x*pdf)
}

lambda = 0.5; beta = 9.0; phi = 0.2; alpha = 0.3

u=seq(0.00000000001,0.999999,0.001)
q=rMKNFW(u, f1, lower = 0, upper = 1, "cumulative")
n=length(q)
l=rep(0,n)
for(i in 1:n){
  b1=integrate(f,0,Inf)$value
  b2=integrate(f,0,q[i])$value
  l[i]=b2/b1
}

f2 = function(x){
  G <- 1 - exp(-(x/phi1)^alpha1)
  H = 1 - (1 - G)^G
  cdf = 1 - exp(-lambda1*(H/(1-H))^beta1)
  return(cdf)
}

lambda1 = 1.0; beta1 = 9.2; phi1 = 0.2; alpha1 = 0.3
q1=rMKNFW(u, f2, lower = 0, upper = 1, "cumulative")

n=length(q1)
l1=rep(0,n)
for(i in 1:n){
  b11=integrate(f,0,Inf)$value
  b21=integrate(f,0,q1[i])$value
  l1[i]=b21/b11
}

f3 = function(x){
  G <- 1 - exp(-(x/phi2)^alpha2)
  H = 1 - (1 - G)^G
  cdf = 1 - exp(-lambda2*(H/(1-H))^beta2)
  return(cdf)
}

lambda2 = 1.5; beta2 = 9.4; phi2 = 0.2; alpha2 = 0.3
q2=rMKNFW(u, f3, lower = 0, upper = 1, "cumulative")
n=length(q2)
l2=rep(0,n)
for(i in 1:n){
  b12=integrate(f,0,Inf)$value
  b22=integrate(f,0,q2[i])$value
  l2[i]=b22/b12
}

f4 = function(x){
  G <- 1 - exp(-(x/phi3)^alpha3)
  H = 1 - (1 - G)^G
  cdf = 1 - exp(-lambda3*(H/(1-H))^beta3)
  return(cdf)
}

lambda3 = 2.0; beta3 = 9.6; phi3 = 0.2; alpha3 = 0.3
q3=rMKNFW(u, f4, lower = 0, upper = 1, "cumulative")

n=length(q3)
l3=rep(0,n)
for(i in 1:n){
  b13=integrate(f,0,Inf)$value
  b23=integrate(f,0,q3[i])$value
  l3[i]=b23/b13
}

plot(c(0,1), c(0,1), type="n", xlab=expression(paste(nu)), ylab=expression(paste(L(nu))))

lines(sort(u, decreasing = F),l,col='black', lty=1, lwd=2)
lines(sort(u, decreasing = F),l1,col='darkred', lty=2, lwd=2)
lines(sort(u, decreasing = F),l2,col='darkblue', lty=3, lwd=2)
lines(sort(u, decreasing = F),l3,col='darkgreen', lty=4, lwd=2)


legend("topleft", expression(paste(lambda,"=0.5; " , beta,"=9.0"), 
                             paste(lambda,"=1.0; " , beta,"=9.2"),
                             paste(lambda,"=1.5; " , beta,"=9.4"), 
                             paste(lambda,"=2.0; " , beta,"=9.6")),
       lty=c(1,2,3,4), lwd=c(2,2,2,2),col=c("black",'darkred','darkblue','darkgreen'), 
       bty="n", cex=1)


######################curvas de Bonferroni########################
rm(list = ls())

rMKNFW = function (n = 1, f, lower = -Inf, upper = Inf, kind = "density") 
{
  if (!is.finite(f(lower))) 
    lower <- lower + .Machine$double.eps
  if (!is.finite(f(upper))) 
    upper <- upper - .Machine$double.eps
  if (!is.numeric(lower) || !is.numeric(upper) || lower >= 
      upper) 
    stop("lower < upper is not fulfilled")
  cumulative <- function(f, lower) {
    function(z) integrate(f, lower, upper = z, rel.tol = 1e-10, 
                          subdivisions = 1000)$value
  }
  if (kind == "density") {
    f.distr <- cumulative(f = f, lower = lower)
    limits <- support.facto(f, lower = lower, upper = upper)
    lower <- max(c(lower, limits[1]))
    upper <- min(c(upper, limits[2]))
  }
  else {
    f.distr <- f
  }
  cum.inverse <- inverse(f.distr, lower = lower, upper = upper)
  sample.unif <- n
  sample.function <- sapply(sample.unif, cum.inverse)
  return(sample.function)
}

f1 = function(x){
  G <- 1 - exp(-(x/phi)^alpha)
  H = 1 - (1 - G)^G
  cdf = 1 - exp(-lambda*(H/(1-H))^beta)
  return(cdf)
}

f <- function(x){
  g <- (alpha/phi) * (x/phi)^(alpha-1) * exp(-(x/phi)^alpha)
  G <- 1 - exp(-(x/phi)^alpha)
  h = g*(1-G)^G*(G/(1-G) - log(1-G))
  H = 1 - (1 - G)^G
  pdf = lambda*beta*h*(H)^(beta-1)*(1-H)^(-(beta+1))*exp(-lambda*(H/(1-H))^beta)
  return(x*pdf)
}

lambda = 0.5; beta = 9.0; phi = 0.2; alpha = 0.3

u=seq(0.00000000001,0.999999,0.001)
q=rMKNFW(u, f1, lower = 0, upper = 1, "cumulative")
n=length(q)
l=rep(0,n)
for(i in 1:n){
  b1=integrate(f,0,Inf)$value
  b2=integrate(f,0,q[i])$value
  l[i]=b2/(u[i]*b1)
}

f2 = function(x){
  G <- 1 - exp(-(x/phi1)^alpha1)
  H = 1 - (1 - G)^G
  cdf = 1 - exp(-lambda1*(H/(1-H))^beta1)
  return(cdf)
}

lambda1 = 1.0; beta1 = 9.2; phi1 = 0.2; alpha1 = 0.3
q1=rMKNFW(u, f2, lower = 0, upper = 1, "cumulative")

n=length(q1)
l1=rep(0,n)
for(i in 1:n){
  b11=integrate(f,0,Inf)$value
  b21=integrate(f,0,q1[i])$value
  l1[i]=b21/(u[i]*b11)
}

f3 = function(x){
  G <- 1 - exp(-(x/phi2)^alpha2)
  H = 1 - (1 - G)^G
  cdf = 1 - exp(-lambda2*(H/(1-H))^beta2)
  return(cdf)
}

lambda2 = 1.5; beta2 = 9.4; phi2 = 0.2; alpha2 = 0.3
q2=rMKNFW(u, f3, lower = 0, upper = 1, "cumulative")

n=length(q2)
l2=rep(0,n)
for(i in 1:n){
  b12=integrate(f,0,Inf)$value
  b22=integrate(f,0,q2[i])$value
  l2[i]=b22/(u[i]*b12)
}

f4 = function(x){
  G <- 1 - exp(-(x/phi3)^alpha3)
  H = 1 - (1 - G)^G
  cdf = 1 - exp(-lambda3*(H/(1-H))^beta3)
  return(cdf)
}

lambda3 = 2.0; beta3 = 9.6; phi3 = 0.2; alpha3 = 0.3
q3=rMKNFW(u, f4, lower = 0, upper = 1, "cumulative")


n=length(q3)
l3=rep(0,n)
for(i in 1:n){
  b13=integrate(f,0,Inf)$value
  b23=integrate(f,0,q3[i])$value
  l3[i]=b23/(u[i]*b13)
}

plot(c(0.0,1), c(0,1), type="n", xlab=expression(paste(nu)), ylab=expression(paste(B(nu))))

lines(sort(u, decreasing = F),l,col='black', lty=1, lwd=2)
lines(sort(u, decreasing = F),l1,col='darkred', lty=2, lwd=2)
lines(sort(u, decreasing = F),l2,col='darkblue', lty=3, lwd=2)
lines(sort(u, decreasing = F),l3,col='darkgreen', lty=4, lwd=2)


legend("topleft", expression(paste(lambda,"=0.5; " , beta,"=9.0"), 
                             paste(lambda,"=1.0; " , beta,"=9.2"),
                             paste(lambda,"=1.5; " , beta,"=9.4"), 
                             paste(lambda,"=2.0; " , beta,"=9.6")),
       lty=c(1,2,3,4), lwd=c(2,2,2,2),col=c('black','darkred','darkblue','darkgreen'), 
       bty="n", cex=1)


################################################################################
