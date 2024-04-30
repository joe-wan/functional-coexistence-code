#!/usr/bin/env Rscript
# 03_cedar_creek/03_appendix_plots.R
# ==================================
# Author: Joe Wan
# Generates extra appendix figures for the Cedar Creek data (Figures S7.2-3).

# Data manipulation libraries
library(dplyr)
library(tidyr)
library(stringr)

# Plotting libraries
library(ggplot2)
library(cowplot)


##### 0a. Load params and helpers #####
source('analysis_params.R')
source('graph_params.R')


##### 0b. Load data #####
param.path <- paste0(outpath, '05_param_df.rds')
raw.param.path <- paste0(outpath, '05_raw_params.rds')
monoculture.path <- paste0(outpath, '05_monoculture_data.rds')
competition.path <- paste0(outpath, '05_competition_data.rds')

# Check if data exists and regenerate if not
paths <- c(param.path, raw.param.path, monoculture.path, competition.path)
exists <- file.exists(paths)
if (any(!exists)) {
    cat("File(s) not found: ")
    cat(paste(paths[!exists], collapse=", "))
    cat("\nRegenerating...\n")
    with(list(), source("01_fit_cedar_creek.R"))
}
param.fits <- readRDS(param.path)
raw.params <- readRDS(raw.param.path)
monoculture <- readRDS(monoculture.path)
competition <- readRDS(competition.path)

competition.fits.wide.path <- paste0(outpath, '05_competition_fits_wide.rds')
if (!file.exists(competition.fits.wide.path)) {
  cat("File not found: ", competition.fits.wide.path, "\nRegenerating...\n")
  with(list(), source("02_visualize_fits.R"))
}
competition.fits.wide <- readRDS(competition.fits.wide.path)

# Set up the order of species for plotting
biomass.order <- c('Poa', 'Agropyron')
to.biomass.order <- function(x) factor(x, levels=biomass.order)


##### 1. Plot alpha values (Figure S7.2)  #####
line.size <- 0.75
species <- function(x) if_else(x=='1', 'Poa', 'Agropyron')
alpha.plot <- param.fits %>%
  select(soil.nitrogen, Rstar1, Rstar2, starts_with('a'), starts_with('alpha')) %>%
  gather('key', 'value', starts_with('a'), starts_with('alpha')) %>%
  mutate(
    param=str_sub(key, 1, -3),
    donor=species(str_sub(key, -1, -1)),
    recipient=species(str_sub(key, -2, -2)),
    type=if_else(donor==recipient, 'conspecific', 'heterospecific')) %>%
  arrange(param) %>% # "a" sorts first, so will be the left panel
  mutate(
    value=if_else(param=='alpha', value*sqrt((soil.nitrogen-Rstar1)*(soil.nitrogen-Rstar2)), value),
    param=if_else(param=='alpha', 
      "paste('realized (scaled ', paste(alpha[ij], ')'))", 
      "paste('resource-independent (', b[ij], ')')")) %>% 
  mutate(param=factor(param, levels=unique(param))) %>%
  ggplot(aes(x=soil.nitrogen, y=value)) +
    geom_line(aes(group=key, color=recipient, linetype=type),
      linewidth=line.size) +
    scale_color_manual(name='recipient species, i', 
      values=c(Poa=color.poa, Agropyron=color.agr),
      labels=c(Poa=expr(italic('Poa')), Agropyron=expr(italic('Agropyron')))) +
    scale_linetype_manual(name='interaction type',
      values=c(conspecific='solid', heterospecific='dotted'),
      labels=c(
        expr(paste('conspecific, ', i == j)), 
        expr(paste('heterospecific, ', i != j)))) +
    facet_wrap(~param, labeller=label_parsed, scales='free_y') +
    xlab(nitrogen.label) + ylab('competition strength') +
    coord_cartesian(xlim=data.xlim, 
      ylim=c(3.75, 8.25), expand=F) +
    theme_cowplot() + theme(
      # legend.position='bottom', legend.justification='center',
      strip.background=element_blank(), 
      strip.clip='off',
      strip.text=element_text(size=14),
      legend.title=element_text(size=14),
      legend.text=element_text(size=14))
if (show.plots) print(alpha.plot)

ggsave(paste0(outpath, 'S7_02_alphas.pdf'), 
  alpha.plot + th, 
  device=cairo_pdf, width=out.width*3/2, height=out.height)


##### 2. Plot selection and complementarity (S7.3)  #####
# Helpers for calculating the AP
pop.sd <- function(x) sqrt(mean((x - mean(x, na.rm=T))^2, na.rm=T))
pop.cov <- function(x, y) mean((x - mean(x, na.rm=T))*(y - mean(y, na.rm=T)), na.rm=T)
pmean <- function(x, y) (x + y)/2

# Calculate the AP
ap.df <- competition.fits.wide %>% 
  rowwise %>% mutate(
    K1=1/alpha11, K2=1/alpha22, Kbar=(K1+K2)/2,
    RY1=N1/K1, RY2=N2/K2,
    CE=(RY1+RY2-1)*Kbar,
    SE=2*pop.cov(c(RY1, RY2), c(K1, K2)),
    ROY=N1+N2-Kbar, TOY=N1+N2-pmax(K1, K2),
    SE.toy=2*pop.cov(c(RY1-1, RY2), c(K1, K2)))

# Plot the AP
ap.plot <- ap.df %>% expand_grid(relative=c(T, F)) %>%
  transmute(soil.nitrogen, relative,
    CE=CE, 
    SE=if_else(relative, SE, SE.toy),
    OY=if_else(relative, ROY, TOY),
  ) %>% gather('component', 'value', -c(soil.nitrogen, relative)) %>%
  mutate(type=if_else(relative, "relative", "transgressive"),
    facet=if_else(relative, "paste('overyielding relative to ', bar(K))", 
    "paste('transgressive overyielding')")) %>%
  mutate(component=factor(component, levels=c('CE', 'SE', 'OY'))) %>%
ggplot(aes(x=soil.nitrogen)) +
  geom_hline(yintercept=0) +
  geom_line(aes(y=value, color=type, linetype=component), 
    linewidth=line.size) +
  facet_wrap(~facet, scales='free', ncol=1, labeller=label_parsed) +
  scale_linetype_manual(NULL,
    values=c(CE='dotted', SE='dashed', OY='solid'),
    labels=c(CE='complementarity', SE='selection', OY='total')) +
  scale_color_manual(values=c(relative=greens[4], transgressive=greens[2])) +
  coord_cartesian(xlim=data.xlim) +
  scale_x_continuous(expand=c(0,0)) +
  xlab(nitrogen.label) + ylab(expression(paste("biomass effect (", g/m^2, ")"))) +
  # Show the linetype but not the color legend
  guides(color=guide_none()) + 
  theme_cowplot() + theme(
    legend.position='bottom', legend.justification='center',
    legend.title=element_text(size=14),
    legend.text=element_text(size=14),
    legend.key.width=unit(28, 'pt'),
    strip.background=element_blank(),
    strip.text=element_text(size=14))
if(show.plots) print(ap.plot)

ggsave(paste0(outpath, 'S7_03_additive_partition.pdf'), 
  ap.plot + th, 
  device=cairo_pdf, width=out.width, height=out.height*2)