#!/usr/bin/env Rscript
# 01_fct_simulations/heatmap_params.R
# ===================================
# Author: Joe Wan
# Parameters for heatmap visualization: color palettes, axis settings, and resolution.

# Color palette functions (interpolated in Lab space)
green.fun  <- get.palette.fun(c('white', greens[2], greens[4]), entire.L = TRUE)  # Gradient for biomass (green shades)
orange.fun <- get.palette.fun(c('white', oranges[2], oranges[4]), entire.L = TRUE) # Gradient for negative SE/CE
blue.fun   <- get.palette.fun(c('white', blues[2], blues[4]), entire.L = TRUE)    # Gradient for positive SE/CE
gray.fun   <- get.palette.fun(c('white', 'gray90'), entire.L = TRUE)             # Neutral gray gradient

# Heatmap grid resolution
heatmap.res <- 100  # Number of divisions along FD axis

# x-axis label for ND
nd.label <- expression(paste('niche difference, ', -log * ' ' * rho))

# Axis transformations
xtran <- I   # Identity transform for ND (no scaling applied)
ytran <- exp # Exponential transform for FD (convert log scale back to original values)

# Axis tick and label settings
break.inc    <- 0.1     # Spacing for axis tick marks
label.inc    <- 0.5     # Spacing for axis labels
label.format <- "%1.1f" # Format for numeric labels (1 decimal)

# Line styles for annotations
axis.lt <- "dotted" # Baseline (FD = 0)
roy.lt  <- "solid"  # ROY boundary line
toy.lt  <- "dashed" # TOY boundary line
coex.lt <- "dashed" # Coexistence boundaries

# Contour and boundary settings
contour.res <- 4 # Number of steps between Y1 and Y2 for biomass contour lines
fudge <- 1e-4    # Small offset to avoid overlapping contour calculations

# Partition component color scaling
q <- 7 # Number of steps in diverging color scale for SE and CE
k <- 1 - col.to.lab(oranges[2])[1] / 100  # Lightness scaling factor for color interpolation