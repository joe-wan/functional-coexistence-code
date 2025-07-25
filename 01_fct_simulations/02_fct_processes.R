#!/usr/bin/env Rscript
# 01_fct_simulations/02_fct_processes.R
# =====================================
# Author: Joe Wan
# Generates figures illustrating stabilization, equalization, and functional
# equalization (Figure 2a-c).

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

##### 0d. Specific parameters and helpers for varying MCT processes #####
source("process_params.R")
source("process_helpers.R")

##### 1. Make FCT conceptual figure (just for the inset later) #####
# Generate base FCT plot and faded overlay for inset visuals
fct.data <- get.fct.df(seq(min.nd, max.nd, by=resolution),
  Y1=Y1, Y2=Y2, log.version=F) %>%
  mutate(ND=trad.to.log(ND))

fct.plot <- plot.fct.df(fct.data, ylim=fct.y, xlim=trad.to.log(fct.x),
      log.version=T,
      xlab=nd.label, ylab=fr.label) +
    theme_cowplot() + theme(aspect.ratio=ratio)

# Add a fade overlay to highlight inset regions
fct.plot.fade <- fct.plot + fade.box(xlim=trad.to.log(fct.x), ylim=fct.y)

##### 2. Illustrate the three new mechanisms #####

##### 2a. Stabilization (increasing ND) #####
# From process_params.R:
# stab.colors, stab.nd, stab.fr

# Generate simulated data for stabilization mechanism
stab.df <- get.process.df(stab.nd, stab.fr, Y1, Y2) %>%
  mutate(group=fr, group.name=paste(fr.prefix, fr))

# Plot biomass as we vary niche difference
stab.biomass <- plot.process.df(stab.df, expr(trad.to.log(nd)),
    show.Y1=T, show.Y2=T,
    legend.name='fitness ratio', colors=stab.colors,
    xlab=nd.label, ylab=biomass.label) +
  scale_x_continuous(breaks=breaks_width(0.2), expand=c(0,0),
    limits=c(0.025, 0.6+1e-4)) +
  theme_cowplot() + theme(aspect.ratio=ratio)
if (show.plots) print(stab.biomass)

# Plot inset showing where we are in FCT space
stab.inset <- plot.inset(stab.df, fct.plot.fade, expr(nd), colors=stab.colors)
if (show.plots) print(stab.inset)

##### 2b. Equalization (increasing FR) #####
# From process_params.R:
# eq.colors, eq.nd, eq.fr

# Generate simulated data for equalization mechanism
eq.df <- get.process.df(eq.nd, eq.fr, Y1, Y2) %>%
  mutate(group=nd, group.name=paste(nd.prefix, round(trad.to.log(nd),2)))

# Plot biomass as we vary fitness ratio
eq.biomass <- plot.process.df(eq.df, expr(fr),
    show.Y1=T, show.Y2=T,
    legend.name='niche difference', colors=eq.colors,
    xlab=fr.label, ylab=biomass.label) +
  # Add vertical line at optimal FR
  geom_vline(xintercept=sqrt(Y1/Y2), linetype=highlight.lt) +
  scale_x_continuous(breaks=breaks_width(0.1), expand=c(0,0),
                     labels=custom_label('%1.1f', c(0.5, 1, 1.5, 2))) +
  theme_cowplot() + theme(aspect.ratio=ratio)
if (show.plots) print(eq.biomass)

# Plot inset showing where we are in FCT space
eq.inset <- plot.inset(eq.df, fct.plot.fade, expr(fr), colors=eq.colors)
if (show.plots) print(eq.inset)

##### 2c. Functional equalization #####
# From process_params.R:
# fu.colors, fu.nd, fu.fr, fu.yr
# delta.yr, yr.examples

# Generate simulated data for functional equalization
fu.df <- get.process.df(fu.nd, fu.fr, Y1, Y1/fu.yr) %>%
  mutate(yr=Y1/Y2, group=Y2, group.name=paste(yr.prefix, yr))

# Now specify points to highlight (special for this plot)
delta.yr <- 0.25
yr.examples <- c(-delta.yr, 0, delta.yr) + Y1/Y2

# Plot biomass as we vary yield ratio
fu.biomass <- plot.process.df(fu.df, expr(Y2),
    show.Y1=T, show.Y2=F, show.oy=T,
    legend.name='yield ratio', colors=fu.colors,
    highlight.condition=expr(yr %in% yr.examples),
    xlab=yield.label, ylab=biomass.label) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(breaks=breaks_width(0.1)) +
  theme_cowplot() + theme(aspect.ratio=ratio)
if (show.plots) print(fu.biomass)

# Inset: recalculate overyield boundaries at different yield ratios
new.toy <- expand_grid(nd=seq(0,1-exp(-max.nd), by=resolution), yr=yr.examples) %>%
  mutate(rho=1-nd, Y1=Y1, Y2=Y1/yr,
         toy.min=case_when(rho<=1/sqrt(Y1/Y2)~Y1/Y2*rho),
         toy.max=if_else(yr==max(yr.examples), 1/rho, pmin(1/rho, (yr+delta.yr)*rho)))

# Plot inset showing where we are in FCT space
fu.inset <- fct.plot.fade +
  geom_ribbon(aes(x=trad.to.log(nd), ymin=toy.min, ymax=toy.max, fill=as.factor(yr)), data=new.toy) +
  geom_text(aes(x=x, y=y), label='★', color='white', size=4, data=data.frame(x=fu.nd, y=fu.fr)) +
  scale_fill_manual(values=fu.colors, labels=parse_format()) +
  theme_void() + theme(legend.position="none", aspect.ratio=ratio)
if (show.plots) print(fu.inset)

##### 3. Export plots #####
# All panels have consistent height; export with same dimensions for proportionality
ggsave(paste0(outpath, '02_a_stabilization_biomass.pdf'), stab.biomass + th,
       device=cairo_pdf, width=out.width, height=out.height)
ggsave(paste0(outpath, '02_a_stabilization_inset.pdf'), stab.inset + th,
       device=cairo_pdf, width=out.width/2, height=out.height/2)
ggsave(paste0(outpath, '02_b_equalization_biomass.pdf'), eq.biomass + th,
       device=cairo_pdf, width=out.width, height=out.height)
ggsave(paste0(outpath, '02_b_equalization_inset.pdf'), eq.inset + th,
       device=cairo_pdf, width=out.width/2, height=out.height/2)
ggsave(paste0(outpath, '02_c_functional_biomass.pdf'), fu.biomass + th,
       device=cairo_pdf, width=out.width, height=out.height)
ggsave(paste0(outpath, '02_c_functional_inset.pdf'), fu.inset + th,
       device=cairo_pdf, width=out.width/2, height=out.height/2)
