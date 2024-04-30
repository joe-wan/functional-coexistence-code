#!/usr/bin/env Rscript
# 01_fct_simulations/03_extended_fct_space.R
# ==========================================
# Author: Joe Wan
# Generates FCT space figures with additional boundaries (Figures 3d, S2.1).

# Data manipulation
library(dplyr)
library(tidyr)

# Plotting
library(ggplot2)
library(cowplot)
library(scales)


##### 0. Load helpers and params #####
# Helper functions for MCT calculations
source("mct_helpers.R")

# Define parameters of system
source("mct_params.R")

# Define visual language and graph options
source("graph_params.R")

# Helper functions for plotting
source("graph_helpers.R")


##### 1. Implement the math for the new form of overyielding #####
relative.oy.boundaries <- function(yr, rho) {
  B <- (1+yr)*(rho+1/rho)/4
  H <- sqrt(B^2-yr)
  return(B + c(-H, H))
}

# Uses the general "get.fct.df" function to calculate the boundaries, then adds
# the new boundaries
get.roy.df <- function(nds, Y1=NULL, Y2=NULL, log.version=F) {
  result <- get.fct.df(nds, Y1=Y1, Y2=Y2, log.version=T) %>%
    rowwise %>% mutate(
      roy.min=log(relative.oy.boundaries(Y1/Y2, exp(-ND))[1]),
      roy.max=log(relative.oy.boundaries(Y1/Y2, exp(-ND))[2])
    ) %>% ungroup
  # If needed, undo the log transformation
  if (!log.version) result <- result %>% 
    mutate(across(-ND, ~ exp(.x))) %>% 
    mutate(ND=1-exp(-ND))
  return(result)
}

plot.roy.df <- function(df, ...) {
  result <- plot.fct.df(df, ...) +
    geom_ribbon(aes(x=ND, ymin=roy.min, 
      ymax=if_else(!is.na(toy.min), toy.min, coex.max)), fill=greens[3]) +
    geom_line(aes(x=ND, y=if_else(ND>0, roy.min, NA)), linetype='dashed')
  # Lines must appear on top of ribbons
  result <- bring.geoms.to.top(result, c('GeomLine', 'GeomVline', 'GeomHline'))
  return(result)
}


##### 2a. Make FCT conceptual figure #####
# Get data frame with boundaries
fct.data <- get.roy.df(
  seq(min.nd, max.nd, by=resolution), Y1=Y1, Y2=Y2)

# Generate plot
fct.plot <- plot.roy.df(fct.data, 
    xlab=nd.label, ylab=fr.label,
    ylim=c(min.fr, max.fr), 
    xlim=1-exp(-c(min.nd, max.nd))) +
  theme_cowplot() + 
  theme(aspect.ratio=ratio)
if (show.plots) print(fct.plot)


##### 2b. Make wider FCT conceptual figure #####
# Get data frame with boundaries
fct.data.wide <-  get.roy.df(
  seq(min.nd, max.nd.wide, by=resolution), Y1=Y1, Y2=Y2)

# Generate plot
fct.plot.wide <- plot.roy.df(fct.data.wide,  
    xlab=nd.label, ylab=fr.label,
    ylim=c(min.fr, max.fr), 
    xlim=1-exp(-c(min.nd, max.nd.wide))) +
  cowplot::theme_cowplot() + 
  theme(aspect.ratio=ratio.wide)
if (show.plots) print(fct.plot.wide)


##### 3. Extra figures for appendix #####
fct.approx <- fct.data.wide %>% get.log.version() %>%
  mutate(rho=exp(-ND)) %>%
  mutate(roy.approx=log(Y1*2/(Y1+Y2)*rho)) %>%
  select(-rho) %>% undo.log.version()

fct.approx.plot <- fct.plot.wide +
  geom_line(aes(x=ND, y=if_else(ND<0, roy.min, NA)),
    data=fct.approx, 
    linetype='dashed') +
  geom_line(aes(x=ND, y=roy.max),
    data=fct.approx, 
    linetype='dashed') +
  geom_line(aes(x=ND, y=roy.approx),
    data=fct.approx, 
    linetype='dotted') +
  geom_hline(yintercept=1, linetype='dotted')
if (show.plots) print(fct.approx.plot)


##### 4. Export plots #####
# All the panels have the same height, so we can keep them in proportion
# by exporting with the same dimensions
ggsave(paste0(outpath, '03_d_fct_extended_space.pdf'), fct.plot + th, 
       device=cairo_pdf, width=out.width*1.5, height=out.height*1.5)
ggsave(paste0(outpath, '03_d_fct_extended_space_wide.pdf'), fct.plot.wide + th, 
       device=cairo_pdf, width=out.width*1.5, height=out.height*1.5)
ggsave(paste0(outpath, 'S2_01_fct_extended_space.pdf'), fct.approx.plot + th, 
       device=cairo_pdf, width=out.width*1.5, height=out.height*1.5)