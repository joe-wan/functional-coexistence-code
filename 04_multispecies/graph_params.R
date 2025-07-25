#!/usr/bin/env Rscript
# 04_multispecies/graph_params.R
# ===============================
# Author: Joe Wan
# Shared plotting parameters for multispecies analysis.

# Color palettes
oranges <- c('#ffc27d', '#ffa154', '#f37329', '#cc3b02', '#a62100')
greens  <- c('#d1ff82', '#9bdb4d', '#68b723', '#3a9104', '#206b00')
blues   <- c('#8cd5ff', '#64baff', '#3689e6', '#0d52bf', '#002e99')
coex.col <- 'gray65'

# Default line width for colored lines
colored.line.width <- 0.75

# Axis labels and math expressions
nd.label <- expression(paste('niche difference, ', - log * ' ' * rho))
median.nd.label <- expression(paste('median niche difference, ', - log * ' ' * rho))

fd.label <- expression(paste('fitness difference, ', log(f[1] / f[2])))
fr.label <- expression(paste('fitness ratio, ', f[1] / f[2]))
median.fr.label <- expression(paste('median fitness ratio, ', f[1] / f[2]))

K.label <- expression(paste('minimum intrinsic yield, ', K[min]))
K.diff.label <- expression(paste('yield range, ', K[max] - K[min]))

resource.label <- "resource lvl., R"
abundance.label <- expression(paste("abund., ", N[i]))

# Panel axis titles
biomass.label <- expression(paste("total biomass, ", Sigma * N))
biomass.label.short <- expression(paste("tot. biomass, ", Sigma * N))

stabilization.label <- "Stabilization"
equalization.label <- "Equalization (fitness–function)"
prodeq.label <- "Functional equalization"

ref.symbol <- "★"   # Symbol denoting the reference parameters
normal.symbol <- "" # Symbol denoting modified, back-calculated parameters

# Output settings
show.plots <- interactive() # If running in a REPL, show plots
outpath <- 'outputs/'       # Directory to save outputs
if (!stringr::str_ends(outpath, '/')) outpath <- paste0(outpath, '/')
dir.create(outpath, showWarnings = FALSE)

out.width <- 4   # Width of main plots (inches)
out.height <- 3  # Height of main plots (inches)
main.width <- out.width
main.height <- out.height
sads.width <- out.width
sads.height <- out.height * 2 / 3

# Font and theme
library(ggplot2)
library(cowplot)

# Base theme
th <- theme(text = element_text(family = "Helvetica LT Std"))

# Miscellaneous
nudge.factor <- 0.15 # Used for positioning labels in plots
sad.extra <- scale_x_continuous(breaks = seq(5, 15, by = 5))
