library(pracma)
library(Matrix)

n = 1000
d = 500
m = 600
lambda = 50

gamma = d/n
xi = m/n
invmom_1 = (1-gamma)/(xi-gamma)
invmom_2 = (1-gamma)*(gamma^2+xi-2*gamma*xi)/(xi-gamma)^3

S = t(gramSchmidt(matrix(rnorm(n*m),n,m))$Q)

A = matrix(rnorm(n*d),n,d)
A.svd = svd(A)
SU = S %*% A.svd$u
invUSSU = solve(t(SU)%*%SU)
mean(diag(invUSSU))
mean(diag(t(invUSSU)%*%invUSSU))

sv = A.svd$d
inv_rr = solve(t(SU)%*%SU+lambda*diag(sv^(-2)))
mean(diag(inv_rr))
mean(diag(t(inv_rr)%*%inv_rr))

1+lambda/sv^2
