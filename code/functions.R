library(pracma)

form_hadamard = function(p)
{
  n = 2^p
  H = matrix(0, n, n)
  H[n,n] = 1
  if (p > 0)
  {
    Hi = H[n,n]
    i = 1
    for (j in 1:p)
    {
      i = i*2
      H[(n-i+1):(n-i/2),(n-i+1):(n-i/2)] = Hi
      H[(n-i+1):(n-i/2),(n-i/2+1):n] = Hi
      H[(n-i/2+1):n,(n-i+1):(n-i/2)] = Hi
      H[(n-i/2+1):n,(n-i/2+1):n] = -Hi
      Hi = H[(n-i+1):n,(n-i+1):n]
    }
  }
  return(H)
}

SRHT = function(X, m, p)
{
  n = dim(X)[1]
  H = form_hadamard(p)
  B = c()
  while (length(B) == 0)
  {
    B = (1:n)[rbernoulli(n=n, p=m/n)]
  }
  return((H %*% (X[sample(1:n),]*sample(c(1,-1),replace=T,size=n)))[B,]/sqrt(n))
}

sketch_SU = function(U, sketch, m, p)
{
  n = dim(U)[1]
  if (sketch=="haar")
  {
    S = t( gramSchmidt(matrix(rnorm(n*m),n,m))$Q )
    return(S %*% U)
  } else if (sketch=="gaussian") {
    S = matrix(rnorm(n*m,sd=sqrt(1/m)), m, n)
    return(S %*% U)
  } else {
    return(SRHT(U, m, p))
  }
}