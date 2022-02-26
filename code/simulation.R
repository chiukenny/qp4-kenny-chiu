# This file contains the R code used to generate the plots shown
# in the simulation section of the report.
# Only the parameters in the following section should be modified
# when running this file.
# ---------------------------------------------------------


# Simulation and configuration parameters
# ---------------------------------------------------------
# Simulation parameters
ps = 8:11      # Vector of integers p that determines sizes n=2^p to simulate
sims = 50      # Number of simulations to compute RMSE/SE over at each point
save_plot = F  # Save output plots (T/F)

# Configuration parameters
gamma = 0.05  # Ratio d/n
xi = 0.1      # Ratio m/n
lambda = 0.1  # Regularization penalty 

sketch = "haar"  # Sketch type: one of {"haar","srht","gaussian"}
num_Q = 3        # Number of matrices to decouple (number of iterations in algorithm)


# Initialization
# ---------------------------------------------------------
library(tidyverse)
source("functions.R")

set.seed(1)

# Set colours for plotting
S_COLS = list("srht"="blue", "haar"="darkgreen", "gaussian"="red")

# Calculate IHS-optimal step size for convenience
inv_mom_1 = (1-gamma)/(xi-gamma)
inv_mom_2 = (1-gamma)*(gamma^2+xi-2*gamma*xi)/(xi-gamma)^3
a = inv_mom_1 / inv_mom_2


# Simulation
# ---------------------------------------------------------
# Initialize matrices to hold computed traces
trace_QQ = matrix(0, sims, length(ps))
trace_Q_Q = matrix(1, sims, length(ps))

# Run simulation
for (j in 1:length(ps))
{
  p = ps[j]
  
  n = 2^p
  d = ceil(gamma*n)
  m = ceil(xi*n)
  
  A = matrix(rnorm(n*d),n,d)
  A.svd = svd(A)
  A.sv = A.svd$d
  
  for (i in 1:sims)
  {
    ridge = diag(lambda/A.sv^2)
    ridge_fac = diag((1+lambda/A.sv^2)^2)
    
    Qs = diag(d)
    for (k in 1:num_Q)
    {
      SU = sketch_SU(A.svd$u, sketch, m, p)
      Q = diag(d) - a*solve(t(SU)%*%SU+ridge, ridge_fac)
      trace_Q_Q[i,j] = trace_Q_Q[i,j] * sum(Q^2)/d
      Qs = Q %*% Qs
    }
    trace_QQ[i,j] = sum(Qs^2) / d
  }
}

trace_err = trace_QQ - trace_Q_Q
rmse = sqrt(colMeans(trace_err^2))
se = apply(abs(trace_err),2,sd)
err_df = data.frame(n=2^ps, rmse=rmse)

get_rmse = function(i){rmse[i]}
se_ub = function(i){rmse[i]+se[i]}
se_lb = function(i){pmax(rmse[i]-se[i], 0)}

err_df %>%
  ggplot(aes(x=n, y=rmse)) +
  geom_line(size=1, color=S_COLS[[sketch]]) +
  geom_point(size=2, color=S_COLS[[sketch]]) +
  stat_summary(geom="ribbon", aes(y=1:length(ps)),
               fun=get_rmse, fun.min=se_lb, fun.max=se_ub,
               alpha=0.3, fill=S_COLS[[sketch]]) +
  labs(y="RMSE") +
  theme_bw()
if (save_plot)
{
  ggsave(paste("plots/rmse_sketch",sketch,"_xi",xi,".png",sep=""), width=6, height=4, units="cm")
}
