library(pracma)
library(Matrix)
library(tidyverse)

n = 10000
d = 50
m = 500
lambda = 0.5

gamma = d/n
xi = m/n
invmom_1 = (1-gamma)/(xi-gamma)
invmom_2 = (1-gamma)*(gamma^2+xi-2*gamma*xi)/(xi-gamma)^3

S = t(gramSchmidt(matrix(rnorm(n*m),n,m))$Q)
S2 = t(gramSchmidt(matrix(rnorm(n*m),n,m))$Q)

A = matrix(rnorm(n*d),n,d)
A.svd = svd(A)
SU = S %*% A.svd$u
#invUSSU = solve(t(SU)%*%SU)
#mean(diag(invUSSU))
#mean(diag(t(invUSSU)%*%invUSSU))

sv = A.svd$d
inv_rr = solve( t(SU)%*%SU + diag(lambda/sv^2) )
#mean(diag( inv_rr ))
mean(diag( inv_rr %*% diag(1+lambda/sv^2) ))
#mean(diag( t(inv_rr)%*%inv_rr ))
mean(diag( t(inv_rr)%*%inv_rr %*% diag((1+lambda/sv^2)^2) ))

S2U = S2 %*% A.svd$u
inv_rr2 = solve( t(S2U)%*%S2U + diag(lambda/sv^2) )

Q1 = diag(d) - inv_rr %*% diag((1+lambda/sv^2)^2)
Q2 = diag(d) - inv_rr2 %*% diag((1+lambda/sv^2)^2)

mean(diag( t(Q1) %*% t(Q2) %*% Q2 %*% Q1 ))
mean(diag( t(Q2) %*% Q2 )) * mean(diag( t(Q1) %*% Q1 ))




gamma = 0.05
xi = 0.1
lambda = 0.5

Ns = c(1,5,10,20,50)
QQ = matrix(0, 10, 5)
Q_Q = matrix(0, 10, 5)
for (j in 1:length(Ns))
{
  N = Ns[j]
  n = 100*N
  d = gamma*n
  m = xi*n
  
  A = matrix(rnorm(n*d),n,d)
  A.svd = svd(A)
  sv = A.svd$d
  
  for (i in 1:10)
  {
    S = t(gramSchmidt(matrix(rnorm(n*m),n,m))$Q)
    S2 = t(gramSchmidt(matrix(rnorm(n*m),n,m))$Q)
    
    SU = S %*% A.svd$u
    S2U = S2 %*% A.svd$u
    
    Q1 = diag(d) - solve(t(SU)%*%SU + diag(lambda/sv^2)) %*% diag((1+lambda/sv^2)^2)
    Q2 = diag(d) - solve(t(S2U)%*%S2U + diag(lambda/sv^2)) %*% diag((1+lambda/sv^2)^2)
    
    QQ[i,j] = mean(diag( t(Q1) %*% t(Q2) %*% Q2 %*% Q1 ))
    Q_Q[i,j] = mean(diag( t(Q2) %*% Q2 )) * mean(diag( t(Q1) %*% Q1 ))
  }
}

rmse = sqrt(colMeans((QQ-Q_Q)^2))
se = apply(abs(QQ-Q_Q),2,sd)
err_df = data.frame(n=Ns*100,
                    rmse=rmse)
get_rmse = function(i){rmse[i]}
se_ub = function(i){rmse[i]+se[i]}
se_lb = function(i){pmax(rmse[i]-se[i], 0)}
err_df %>%
  ggplot(aes(x=n, y=rmse)) +
  geom_line(size=1, color="blue") +
  geom_point(size=3, color="blue") +
  stat_summary(geom="ribbon", aes(y=1:5),
               fun=get_rmse, fun.min=se_lb, fun.max=se_ub, alpha=0.3, fill="blue") +
  labs(y="RMSE") +
  theme_bw()
