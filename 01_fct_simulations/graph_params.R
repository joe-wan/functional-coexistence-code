# 01_fct_simulations/graph_params.R
# =================================
# Author: Joe Wan
# Shared parameters for graphing in the FCT and MCT simulations.

library(stringr)

# Output options
# --------------
show.plots <- interactive() # Will suppress plot output if not REPL
outpath <- 'outputs/'       # Directory to save plots
out.width <- 4              # Width of main plots
out.height <- 3             # Height of main plots
# Font, etc. will be added by this theme:
th <- theme(text=element_text(family="Helvetica"))

# Ensure the path ends with a slash
if (!str_ends(outpath, '/')) outpath <- paste0(outpath, '/')
# Ensure output directory exists
dir.create(outpath, showWarnings = FALSE)

#  Text labels
# ------------
biomass.label <- expression(paste('total biomass, ', Sigma * N))
# total biomass, ΣN
fr.label <- expression(paste('fitness ratio, ', f[1]/f[2]))
# 'fitness ratio, f₁/f₂'
nd.label <- expression(paste('niche difference, ', - log * ' ' * rho))
# 'niche difference, −log ρ'
yield.label <- expression(paste('intrinsic yield, ', K[2]))
fr.prefix <- 'f[1]/f[2] ==' # 'f₁/f₂ ='
nd.prefix <- '-log * rho ==' # '−log(ρ) ='
yr.prefix <- 'K[1]/K[2] ==' # 'K₁/K₂ ='

# Colors
# ------
oranges <- c('#ffc27d', '#ffa154', '#f37329', '#cc3b02', '#a62100')
greens <- c('#d1ff82', '#9bdb4d', '#68b723', '#3a9104', '#206b00')
blues <- c('#8cd5ff', '#64baff', '#3689e6', '#0d52bf', '#002e99')
reds <- c('#ff8c82', '#ed5353', '#c6262e', '#a10705', '#7a0000')
coex.col <- 'gray65'
ass.col <- 'gray90'

# Line thicknesses, etc.
# ----------------------
axes.lt <- 9
highlight.lt <- 8
line.size <- 0.75

arrow.size <- 1
custom.arrow <- arrow(length = unit(0.02, "npc"))
arrow.avoid <- 0.025