#!/usr/bin/env Rscript
# 01_fct_simulations/03_extended_fct_space.R
# ==========================================
# Author: Joe Wan
# Generates FCT space figures with additional boundaries (Figure S1).

# Data manipulation
library(dplyr)
library(ggplot2)
library(tidyr)
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
# the new relative overyielding boundaries
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

# Define custom color gradient for relative overyielding regions
more.greens <- c("#FFFFFF", "#ECF9D7", "#D8F2B1", "#C4EB8C", "#AFE36B", "#9BDB4D", "#89D233", "#78C81A", "#69BE00")

plot.roy.df <- function(df, ...) {
  result <- plot.fct.df(df, ...) +
    geom_ribbon(aes(x=ND, ymin=roy.min, 
      ymax=if_else(!is.na(toy.min), toy.min, coex.max)), fill=more.greens[4]) +
    geom_line(aes(x=ND, y=if_else(ND>0, roy.min, NA)), linetype='dashed')
  # Lines must appear on top of ribbons
  result <- bring.geoms.to.top(result, c('GeomLine', 'GeomVline', 'GeomHline'))
  return(result)
}

##### 2a. Make FCT conceptual figure #####
# Get data frame with new boundaries
fct.data <- get.roy.df(seq(min.nd, max.nd, by=resolution), Y1=Y1, Y2=Y2) %>%
  mutate(ND=trad.to.log(ND))

# Generate plot
fct.plot <- plot.roy.df(fct.data,
    log.version=T,
    xlab=nd.label, ylab=fr.label,
    ylim=c(min.fr, max.fr),
    xlim=trad.to.log(1-exp(-c(min.nd, max.nd)))) +
  theme_cowplot() +
  theme(aspect.ratio=ratio)
if (show.plots) print(fct.plot)

##### 2b. Make wider FCT conceptual figure #####
# Get data frame with extended ND range
fct.data.wide <- get.roy.df(seq(min.nd, max.nd.wide, by=resolution), Y1=Y1, Y2=Y2) %>%
  mutate(ND=trad.to.log(ND))

# Generate wider plot
fct.plot.wide <- plot.roy.df(fct.data.wide,
    log.version=T,
    xlab=nd.label, ylab=fr.label,
    ylim=c(min.fr, max.fr),
    xlim=trad.to.log(1-exp(-c(min.nd, max.nd.wide)))) +
  theme_cowplot() +
  theme(aspect.ratio=ratio.wide)
if (show.plots) print(fct.plot.wide)

##### 3. Extra figures for appendix #####
# Compute approximate relative overyielding boundaries and overlay on plot
fct.approx <- fct.data.wide %>%
  mutate(rho=exp(-ND)) %>%
  mutate(roy.approx=(Y1*2/(Y1+Y2)*rho)) %>%
  select(-rho)

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
ggsave(paste0(outpath, 'S01_fct_extended_space.pdf'), fct.approx.plot + th, 
       device=cairo_pdf, width=out.width*1.5, height=out.height*1.5)