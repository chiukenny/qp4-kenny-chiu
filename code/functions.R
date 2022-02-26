# This file contains the user-defined functions used in simulation.R.
# This file should not be modified if only reproducing the work.
# ---------------------------------------------------------


library(pracma)


# Form the n x n Walsh-Hadamard matrix
#   p: the exponent in n=2^p
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

# Generate a SRHT sketch
#   X: the matrix to sketch
#   m: the sketch size
#   p: the exponent in n=2^p
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

# Generate a sketch
#   U: the matrix to sketch
#   sketch: one of {"haar", "gaussian", "srht"}
#   m: the sketch size
#   p: the exponent in n=2^p
sketch_SU = function(U, sketch, m, p)
{
  n = dim(U)[1]
  if (sketch=="haar")
  {
    # Apply Gram-Schmidt to an i.i.d. standard Gaussian matrix
    S = t( gramSchmidt(matrix(rnorm(n*m),n,m))$Q )
    return(S %*% U)
  } else if (sketch=="gaussian") {
    S = matrix(rnorm(n*m,sd=sqrt(1/m)), m, n)
    return(S %*% U)
  } else {
    return(SRHT(U, m, p))
  }
}