#!/usr/bin/env Rscript
# 04_multispecies/multispecies_params.R
# ======================
# Author: Joe Wan
# Defines global parameters and metaparameters for the multispecies model 
# simulations.

# General parameters
n.sp <- 20      # Number of species to simulate
R0 <- 100       # Initial resource concentration

# Metaparameters for random generation of traits
base.cv <- 0.1      # Coefficient of variation for all traits
mu.mean <- 0.5      # Mean mortality rate
up.mean <- 2        # Mean resource uptake ability
eff.mean <- 0.05    # Mean resource use efficiency
B.inter.mean <- 1   # Mean interspecific interference
B.intra.mean <- 5   # Mean intraspecific interference

# Combine into a single list for convenience
metaparams <- list(
  base.cv = base.cv,
  mu.mean = mu.mean,
  up.mean = up.mean,
  eff.mean = eff.mean,
  B.inter.mean = B.inter.mean,
  B.intra.mean = B.intra.mean
)