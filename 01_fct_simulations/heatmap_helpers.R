#!/usr/bin/env Rscript
# 01_fct_simulations/heatmap_helpers.R
# ====================================
# Author: Joe Wan
# Helper functions for color interpolation, transformations, and finalizing heatmap aesthetics.
# These utilities support generating smooth color palettes in perceptual color space (Lab),
# transforming log scales, and applying consistent plot styling for heatmaps.

# Convert a color from sRGB (hex or name) to CIE Lab color space
col.to.lab <- function(col) {
  # Convert hex or named color to RGB (0–1 scale)
  rgb <- as.vector(col2rgb(col)) / 255
  # Convert RGB to Lab color space for perceptual uniformity
  return(convertColor(rgb, from = 'sRGB', to = 'Lab'))
}

# Convert a color from CIE Lab back to sRGB and return hex representation
lab.to.col <- function(lab) {
  # Convert Lab to RGB
  rgb <- convertColor(lab, from = 'Lab', to = 'sRGB')
  # Convert RGB vector to hex color code
  return(rgb(rgb[1], rgb[2], rgb[3], maxColorValue = 1))
}

# Create a palette interpolation function in Lab space
get.palette.fun <- function(cols, entire.L = FALSE) {
  # Convert input colors to Lab color space
  cols.lab <- lapply(cols, col.to.lab)
  
  # Extract individual L, a, and b components
  Ls <- sapply(cols.lab, function(x) x[1])
  as <- sapply(cols.lab, function(x) x[2])
  bs <- sapply(cols.lab, function(x) x[3])

  # Interpolation functions for a and b channels (L used as x)
  a.fun <- splinefun(Ls, as, method = 'natural')
  b.fun <- splinefun(Ls, bs, method = 'natural')

  # Return a function to compute interpolated colors for x ∈ [0, 1]
  get.col <- function(x) {
    # Handle vector input by recursion
    if (is.vector(x) && length(x) > 1) 
      return(sapply(x, get.col))

    # Determine L value based on position and interpolation mode
    if (entire.L) {
      L <- 100 * (1 - x)   # Full L range [0,100]
    } else {
      L <- max(Ls) - x * diff(range(Ls))  # Subrange of L from palette
    }

    # Interpolate a and b for given L
    a <- a.fun(L)
    b <- b.fun(L)

    # Convert back to RGB and return hex
    rgb <- convertColor(c(L, a, b), from = 'Lab', to = 'sRGB')
    return(rgb(rgb[1], rgb[2], rgb[3], maxColorValue = 1))
  }

  return(get.col)
}

# Transformation: Convert log base e (natural log) to log base 10
# Useful for ggplot2 axis transformations
log_to_log10_trans <- function() {
  trans_new("log_to_log10", 
    transform = function(x) x * log(10), 
    inverse   = function(x) x / log(10))
}

# Finalize heatmap styling: add reference lines, axis formatting, and theme adjustments
finish.heatmap <- function(x, lt = coex.lt) {
  x +
    # Add reference lines: horizontal (FD = 0) and coexistence boundaries
    geom_line(aes(x = xtran(ND), y = ytran(0)), linetype = axis.lt) +
    geom_line(aes(x = xtran(ND), y = ytran(ND)), linetype = lt) +
    geom_line(aes(x = xtran(ND), y = ytran(-ND)), linetype = lt) +

    # Use log scaling on y-axis (functional ratios)
    coord_trans(
      x = "identity", y = "log",
      xlim = c(0, xtran(max.nd.heatmap)),
      ylim = c(ytran(log(min.fr)), ytran(log(max.fr))),
      expand = FALSE
    ) +

    # Apply theme and axis settings
    theme_cowplot() +
    scale_x_continuous(
      breaks = breaks_width(break.inc), expand = c(0, 0),
      labels = get.label(c(0, xtran(max.nd.heatmap)), label.inc, label.format)
    ) +
    scale_y_continuous(
      breaks = breaks_width(break.inc), expand = c(0, 0),
      labels = get.label(c(ytran(log(min.fr)), ytran(log(max.fr))), label.inc, label.format)
    ) +

    # Axis labels
    xlab(nd.label) + ylab(fr.label) +

    # Adjust legend ticks and plot aspect ratio
    theme(
      legend.ticks.length = unit(c(-3, 0), 'pt'),
      legend.ticks = element_line(color = 'black'),
      aspect.ratio = heatmap.aspect
    )
}