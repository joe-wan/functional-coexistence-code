# 03_cedar_creek/analysis_params.R
# ================================
# Author: Joe Wan
# Shared options for the Cedar Creek analysis.

library(stringr)

outpath <- 'outputs/'
# Ensure the path ends with a slash
if (!str_ends(outpath, '/')) outpath <- paste0(outpath, '/')
# Ensure output directory exists
dir.create(outpath, showWarnings = FALSE)

# Only show plots if we're in an interactive session
show.plots <- interactive()

nitrogen.values <- seq(0,1500,by=0.1)