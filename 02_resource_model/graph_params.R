# 02_resource_model/graph_params.R
# ================================
# Author: Joe Wan
# Options for graphing.

library(stringr)

outpath <- 'outputs/'
# Ensure the path ends with a slash
if (!str_ends(outpath, '/')) outpath <- paste0(outpath, '/')
# Ensure output directory exists
dir.create(outpath, showWarnings = FALSE)

oranges <- c('#ffc27d', '#ffa154', '#f37329', '#cc3b02', '#a62100')
greens <- c('#d1ff82', '#9bdb4d', '#68b723', '#3a9104', '#206b00')
blues <- c('#8cd5ff', '#64baff', '#3689e6', '#0d52bf', '#002e99')

coex.col <- 'gray65'
ass.col <- 'gray90'

# Only show plots if we're in an interactive session
show.plots <- interactive()

# Output appearance
out.width <- 4
out.height <- 3
th <- theme(text=element_text(family="Helvetica"))

# Label text
fr.label <- expression(paste('fitness ratio, ', f[1]/f[2]))
# 'fitness ratio, f₁/f₂'
nd.label <- expression(paste('niche difference, ', 1 - rho))
resource.label <- 'resource level, R'
