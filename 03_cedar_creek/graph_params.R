# 03_cedar_creek/graph_params.R
# =============================
# Author: Joe Wan
# Shared parameters for the Cedar Creek data analysis.

# Color for Poa and Agropyron
oranges <- c('#ffc27d', '#ffa154', '#f37329', '#cc3b02', '#a62100')
greens <- c('#d1ff82', '#9bdb4d', '#68b723', '#3a9104', '#206b00')
blues <- c('#8cd5ff', '#64baff', '#3689e6', '#0d52bf', '#002e99')
color.agr <- oranges[3]
color.poa <- blues[2]

coex.col <- 'gray65'
ass.col <- 'gray90'
toy.col <- greens[2]
toy.alpha <- 0.1

# Thicken lines a bit
line.size <- 0.75
# These x-limits contain all data
data.xlim <- c(60,1400)
# General theme
my.theme <- theme_cowplot() +
    theme(text=element_text(family="Helvetica"))
# General theme
th <- theme(text=element_text(family="Helvetica"))

out.width <- 4*4/3
out.height <- 3
composite.width <- 4
composite.base.height <- 2.5

# Labels
biomass.label <- expression(paste("biomass (", g/m^2, ")"))
nitrogen.label <- expression(paste("soil nitrogen (", "mg"/"kg", ")"))
nd.label <- expression(paste('niche difference, ', 1 - rho))
fr.label <- expression(paste('fitness ratio, ', 
    f[italic(paste("Poa"))]/f[italic(paste("Agr."))]))
# fr.label <- expression(atop('fitness ratio, ', 
#     f[italic(Poa)]/f[italic(Agropyron)]))
delta.label <- expression(paste("biomass diff. (", g/m^2, ")"))
# delta.label <- expression(paste("difference in biomass (", g/m^2, ")"))

