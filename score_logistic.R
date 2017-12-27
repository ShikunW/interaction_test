# logistic ---------
# 1. generate the data
logit.inv = function(x) 1/(1+exp(-x))
n=1000
x1 = rnorm(n,1,1)
x2 = rnorm(n,2,3)
y = apply(matrix(1 + 2*x1 +3*x2,ncol=1),1,function(x) rbinom(1,1,logit.inv(x)))

mod = glm(y~x1 + x2 + I(x1*x2),family = binomial)
mod.null = glm(y~x1 + x2,family = binomial)

# wald -----

b_fit = matrix(coef(mod),ncol=1)
mu = logit.inv(x%*%b_fit)
res = y - fitted(mod)
x = matrix(cbind(rep(1,n),x1,x2,x1*x2),ncol=4,byrow=F)

O = matrix(0,nrow=4,ncol=4)
B_mb = matrix(0,nrow=4,ncol=4)
B_san = matrix(0,nrow=4,ncol=4)
for (i in 1:n){
  xi = matrix(x[i,],nrow=1)
  O = O + mu[i] * (1-mu[i]) * t(xi) %*% xi
  B_mb = B_mb + mu[i] * (1-mu[i]) * t(xi) %*% xi
  B_san = B_san + t(xi)%*%xi * res[i]^2
}
Oinv = solve(O)
V_mb = Oinv %*% B_mb %*% Oinv
V_san = Oinv %*% B_san %*% Oinv
Wald_mb = b_fit[4]^2/V_mb[4,4]
Wald_san = b_fit[4]^2/V_san[4,4]
1-pchisq(Wald_mb,df = 1)
coef(summary(mod))[4,4]
1-pchisq(Wald_san,df = 1)
coeftest(mod,df=n-4,vcov. = sandwich)['I(x1 * x2)',4]

# score -------
x0 = x[,1:3]
x3 = x[,4]
b0_fit = matrix(coef(mod.null),ncol=1)
mu0_hat = logit.inv(x0%*%b0_fit)
res = residuals(mod.null) # y - mu_hat
U3 = sum(x[,4]*res)
D_san = matrix(0,ncol=4,nrow=4)
A = matrix(0,ncol=4,nrow=4)
t1 = matrix(0,nrow=1,ncol=3);t2 = matrix(0,nrow=3,ncol=3)
D_mb = matrix(0,ncol=4,nrow=4);D_san = matrix(0,ncol=4,nrow=4)
for (i in 1:n) {
  x0i = matrix(x0[i,],nrow=1)
  t1 = t1 + x3[i] * mu0_hat[i] * (1-mu0_hat[i]) * x0i
  t2 = t2 + mu0_hat[i] * (1-mu0_hat[i]) * t(x0i) %*% x0i
  D_mb = D_mb + x[i,] %*% t(x[i,]) * mu0_hat[i]*(1-mu0_hat[i]) ####### this is incorrect??
  D_san = D_san + x[i,] %*% t(x[i,]) * res[i]^2
}

A = matrix(c(-t1 %*% solve(t2),1),nrow=1)
vs_mb = A %*% D_mb %*% t(A)
vs_san = A %*% D_san %*% t(A)

Score_mb = U3^2/vs_mb
1 - pchisq(Score_mb,df = 1) # model-based score test
anova(mod,mod.null,test='Rao')[2,'Pr(>Chi)'] # model-based score test by package

vs_san = A %*% D_san %*% t(A)
Score_san = U3^2/vs_san
1 - pchisq(Score_san,df = 1) #sandwich score test
