# 01_fct_simulations/mct_helpers.R
# ================================
# Author: Joe Wan
# Helper functions for calculating and visualizing MCT metrics.

# Calculate alphas from given MCT metrics
get.alphas <- function(rho, fr, Y1=1, Y2=1, focal.sp=1) {
  # Ensure that focal.sp is 1 or 2; if it is 2, recalculate fr with 1 as focal species
  if (focal.sp==2) {
    fr <- 1/fr
  } else if (focal.sp!=1) {
    stop('Focal species must be 1 or 2')
  }
  
  # Yield is the inverse of intraspecific competition, so:
  alpha11 <- 1/Y1; alpha22 <- 1/Y2
  # We use the following identity to calculate interspecific competition: 
  #   rho * f1/f2 = a21 / a11
  # The ratio determines invasion success and corresponds to Carroll (2011)'s S_2
  alpha21 <- alpha11 * rho * fr
  alpha12 <- alpha22 * rho / fr
  
  return(list(alpha11=alpha11, alpha12=alpha12,
              alpha21=alpha21, alpha22=alpha22))
}

# Given alphas, calculate equilibrium
get.equilibrium <- function(alpha11, alpha12, alpha21, alpha22, censor=T) {
  if (typeof(alpha11) == 'list') return(with(alpha11, get.equilibrium(alpha11, alpha12, alpha21, alpha22, censor)))
  
  # Calculate equilibrium. One way to do this is Cramer's rule,
  # using matrix determinants: given the matrix A with A_i,j = alpha_i,j
  # we calculate the equilibrium N_i = det(A_i)/det(A) where A_i is obtained
  # by replacing column i of A with ones.
  
  # Writing out the determinants gives:
  dt <- alpha11*alpha22 - alpha12*alpha21
  N1 <- (alpha22-alpha12)/dt
  N2 <- (alpha11-alpha21)/dt
  
  # If the "censor" option is false, just return the calculation
  if (!censor) return(list(N1=N1, N2=N2))
  
  # Otherwise check if equilibrium is stable/feasible, then return the actual outcome
  if (dt<=0) return(list(N1=as.double(NA), N2=as.double(NA))) # Unstable (ASS)
  if (N1>0 & N2<0) return(list(N1=1/alpha11, N2=0)) # N1 wins
  if (N1<0 & N2>0) return(list(N1=0, N2=1/alpha22)) # N2 wins
  return(list(N1=N1, N2=N2)) # Coexistence
}