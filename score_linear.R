# ------------------------------------------------------
# 1. generate the data
n=1000
x1 = rnorm(n,1,1)
x2 = rnorm(n,2,3)
y = (1 + 2*x1 +3*x2) + rnorm(n)

mod = glm(y~x1 + x2 + I(x1*x2))
mod.null = glm(y~x1 + x2)

# wald -----

B_san = matrix(0,nrow=4,ncol=4)
b_fit = matrix(coef(mod),ncol=1)
res = residuals(mod) # res = y - matrix(x%*%b_fit,ncol=1)
Vi = sum(res^2)/(n-4)
x = matrix(cbind(rep(1,n),x1,x2,x1*x2),ncol=4,byrow=F)
Oinv = solve(t(x) %*% x / Vi)
V_mb = Oinv

for (i in 1:n){
  xi = matrix(x[i,],nrow=1)
  B_san = B_san + t(xi)%*%xi * res[i]^2 / Vi^2
}

V_san = Oinv %*% B_san %*% Oinv
Wald_mb = b_fit[4]^2/V_mb[4,4]
Wald_san = b_fit[4]^2/V_san[4,4]
1-pchisq(Wald_mb,df = 1)
coef(summary(mod))[4,4]
1-pchisq(Wald_san,df = 1)
coeftest(mod,df=n-4,vcov. = sandwich)['I(x1 * x2)',4]

# score -------
x0 = matrix(cbind(rep(1,n),x1,x2),ncol=3)
x3 = matrix(x1*x2,ncol=1)
b0_fit = matrix(coef(mod.null),ncol=1)
res0 = y - matrix(x0%*%b0_fit,ncol=1)
Vi = sum(res0^2)/(n-4)
U3 = sum(x3*res0) / Vi

t1 = matrix(0,nrow=1,ncol=3);t2 = matrix(0,nrow=3,ncol=3)
D_mb = matrix(0,ncol=4,nrow=4);D_san = matrix(0,ncol=4,nrow=4)
for(i in 1:n){
  x0i = matrix(x0[i,],nrow=1)
  x3i = x3[i]
  xi = matrix(c(x0i,x3i),nrow=1)
  t1 =  t1 + x3i * x0i / Vi
  t2 = t2 + t(x0i) %*% x0i / Vi
  D_mb = D_mb + t(xi) %*% xi / Vi
  D_san = D_san + t(xi) %*% xi * res0[i]^2 / Vi^2
}
D_mb = D_mb
A = matrix(c(-t1 %*% solve(t2),1),nrow=1)

vs_mb = A %*% D_mb %*% t(A)
Score_mb = U3^2/vs_mb
1 - pchisq(Score_mb,df = 1) # model-based score test
anova(mod,mod.null,test='Rao')[2,'Pr(>Chi)'] # model-based score test by package

vs_san = A %*% D_san %*% t(A)
Score_san = U3^2/vs_san
1 - pchisq(Score_san,df = 1) # sandwich score test