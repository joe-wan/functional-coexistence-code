#!/usr/bin/env Rscript
# 01_fct_simulations/process_params.R
# ===================================
# Author: Joe Wan
# Shared parameters for analyses of the effect of changing MCT/FCT components.

# Colors and ranges
stab.colors <- blues[4:1]
stab.nd <- seq(0.025, 0.5, by=resolution)
stab.fr <- seq(0.93, 1.23, len=4)

eq.colors <- oranges[5:2]
eq.nd <- seq(0.15, 0.36, len=4)
eq.fr <- seq(0.6, 1.7, by=resolution)

fu.colors <- greens[2:4]
fu.nd <- eq.nd[3] # Was 0.3, but it's so close to one of the above values
fu.fr <- stab.fr[3] # Was 1.1, but it's so close to one of the above values
fu.yr <- seq(-0.3, +0.3, by=resolution) + Y1/Y2

# Yield ratio choices
delta.yr <- 0.25
yr.examples <- c(-delta.yr, 0, delta.yr) + Y1/Y2