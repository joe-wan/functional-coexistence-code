#!/usr/bin/env Rscript
# 01_fct_simulations/01_fct_space.R
# =================================
# Author: Joe Wan
# Generates conceptual figures for the MCT and FCT spaces (Figure 1a-b).


# Data manipulation
library(dplyr)
library(tidyr)

# Plotting
library(ggplot2)
library(cowplot)
library(scales)


##### 0a. MCT helper functions #####
source("mct_helpers.R")


##### 0b. Define parameters of system and ranges for visualizations #####
source("mct_params.R")


##### 0c. Define visual language #####
source("graph_params.R")
source("graph_helpers.R")


##### 1. Make MCT conceptual figure (Figure 1a) #####
mct.data <- get.mct.df(seq(-max.nd.lo, max.nd.lo, by=resolution))
mct.plot <- plot.mct.df(mct.data, 
  ylim=c(log(min.fr), log(max.fr)),
  xlab='niche difference', ylab='fitness difference')
if (show.plots) print(mct.plot)

##### 2. Make FCT conceptual figure (Figure 1b) #####
fct.data <- get.fct.df(seq(min.nd, max.nd, by=resolution),
  Y1=Y1, Y2=Y2, log.version=F)
fct.plot <- plot.fct.df(fct.data, ylim=fct.y, xlim=fct.x, 
      zero.linetype=axes.lt, optimal.linetype=highlight.lt,
      xlab=nd.label, ylab=fr.label) +
    theme_cowplot() + theme(aspect.ratio=ratio)
if (show.plots) print(fct.plot)


##### 2b. Make wider FCT conceptual figure #####
fct.data.wide <- get.fct.df(seq(min.nd, max.nd.wide, by=resolution),
  Y1=Y1, Y2=Y2, log.version=F)
fct.plot.wide <- plot.fct.df(fct.data.wide, ylim=fct.y, xlim=fct.x.wide, 
      zero.linetype=axes.lt, optimal.linetype=highlight.lt,
      xlab=nd.label, ylab=fr.label) +
    theme_cowplot() + theme(aspect.ratio=ratio.wide)
if (show.plots) print(fct.plot.wide)


##### 3. Export plots #####
# All the panels have the same height, so we can keep them in proportion
# by exporting with the same dimensions
ggsave(paste0(outpath, '01_a_mct_space.pdf'), mct.plot + th, 
       device=cairo_pdf, width=out.width*1.5, height=out.height*1.5)
ggsave(paste0(outpath, '01_b_fct_space.pdf'), fct.plot + th, 
       device=cairo_pdf, width=out.width*1.5, height=out.height*1.5)
ggsave(paste0(outpath, '01_b_fct_space_wide.pdf'), fct.plot.wide + th, 
       device=cairo_pdf, width=out.width*1.5, height=out.height*1.5)