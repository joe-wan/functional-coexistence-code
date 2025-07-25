#!/usr/bin/env Rscript
# 03_cedar_creek/02_visualize_fits.R
# ==================================
# Author: Joe Wan
# Visualizes the fits of the Cedar Creek data (Figure 8a-c, S5).

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
source('../01_fct_simulations/mct_helpers.R') # Only for get.equilibrium()


##### 0b. Load data #####
param.path <- paste0(outpath, '08_param_df.rds')
raw.param.path <- paste0(outpath, '08_raw_params.rds')
monoculture.path <- paste0(outpath, '08_monoculture_data.rds')
competition.path <- paste0(outpath, '08_competition_data.rds')

# Check if data exists and regenerate if not
paths <- c(param.path, raw.param.path, monoculture.path, competition.path)
exists <- file.exists(paths)
if (any(!exists)) {
    cat("File(s) not found: ")
    cat(paste(paths[!exists], collapse=", "))
    cat("\nRegenerating...\n")
    with(list(), source("01_fit_cedar_creek.R"))
}

# Load the data
param.fits <- readRDS(param.path)
raw.params <- readRDS(raw.param.path)
monoculture <- readRDS(monoculture.path)
competition <- readRDS(competition.path)

# Set up the order of species for plotting
biomass.order <- c('Poa', 'Agropyron')
to.biomass.order <- function(x) factor(x, levels=biomass.order)


##### 1. Calculate outcomes  #####
# Calculate monoculture fits
monoculture.fits <- param.fits %>% 
  crossing(data.frame(species.id=c(1, 2))) %>%
  mutate(species=case_when(
    species.id==1 ~ 'Poa',
    species.id==2 ~ 'Agropyron')) %>%
  gather(key, alphaii, alpha11:alpha22) %>%
  filter(key == paste0('alpha', species.id, species.id)) %>%
  transmute(soil.nitrogen, species, species.id,
    biomass=1/alphaii, K=biomass) %>%
  mutate(species=to.biomass.order(species))

# Wide version of the monoculture fits
monoculture.fits.wide <- monoculture.fits %>%
  transmute(soil.nitrogen, key=paste0("K", species.id), biomass) %>%
  spread(key, biomass)

# Wide version of the empirical data
competition.wide <- competition %>%
  transmute(soil.nitrogen, plot, 
    key=paste0("N", species.id), biomass) %>%
  spread(key, biomass)
# Add fitted monoculture predictions, if needed
competition.calcs.wide <- competition.wide %>%
  mutate(params=as.data.frame(raw.params)) %>%
  unnest(params) %>% 
  mutate(
    K1=(soil.nitrogen-Rstar1)/(f11+g11*soil.nitrogen), 
    K2=(soil.nitrogen-Rstar2)/(f22+g22*soil.nitrogen)) %>%
  mutate(delta=N1+N2-pmax(K1, K2))


##### 2. Visualize #####
# Set up scales for species and outcomes
species.color.scale <- scale_color_manual('',
  values=c(Agropyron=color.agr,  Poa=color.poa), 
  breaks=biomass.order,
  na.value=NA)
outcome.fill.scale <- scale_fill_manual('',
  values=c(coexistence=coex.col, overyielding=toy.col, exclusion='white'),
  guide=guide_legend(override.aes=list(color='black', linetype='dashed')))
linetype.scale <- scale_linetype_manual("",
  values=c(competition='solid',  monoculture='dashed'),
  # Force to show all breaks
  na.value=NA)

# Theme for the plots
bottom.legend.theme <- theme(legend.position='bottom', 
  legend.justification="center",
  legend.text=element_text(size=12, margin=margin(l=-2)), 
  legend.spacing.x=unit(8, "pt"))


##### 1a. Visualize the fits #####
# Visualize the competition fits
# ------------------------------
# Calculate the equilibrium biomass under competition
competition.fits.wide <- param.fits %>% rowwise %>%
  mutate(equilibrium=as.data.frame(
    get.equilibrium(alpha11, alpha12, alpha21, alpha22)
  )) %>% ungroup %>%
  unnest(equilibrium)
competition.fits <- competition.fits.wide %>%
  select(soil.nitrogen, N1, N2) %>%
  gather(species, biomass, N1:N2) %>%
  mutate(species=case_when(
    species=='N1' ~ 'Poa', 
    species=='N2' ~ 'Agropyron'))

# Save the competitive fits
for (suff in c('', '_wide')) {
  fits <- if (suff == '_wide') competition.fits.wide else competition.fits
  saveRDS(fits, paste0(outpath, '08_competition_fits', suff, '.rds'))
  write.csv(fits, file=paste0(outpath, '08_competition_fits', suff, '.csv'), row.names=F)
}

# Plot the competition fits
competition.plot <- ggplot(competition.fits, 
  aes(x=soil.nitrogen, y=biomass, color=species)) +
  geom_line(size=line.size, linetype="dashed",
    data=monoculture.fits) +
  geom_line(size=line.size) + 
  geom_point(data=competition, size=1) +
  species.color.scale +
  xlim(data.xlim) + # Remove out of x range data
  coord_cartesian(xlim=data.xlim, ylim=c(0,225), expand=F, clip="off") +
  xlab(nitrogen.label) + ylab(biomass.label) +
  my.theme
if (show.plots) print(competition.plot)

# Visualize the monoculture fits (for appendix)
# ---------------------------------------------
monoculture.plot <- ggplot(monoculture.fits, aes(x=soil.nitrogen, color=species)) +
  geom_line(aes(y=biomass), size=line.size) +
  geom_point(aes(y=biomass), data=monoculture) +
  species.color.scale +
  xlim(data.xlim) + # Remove out of x range data
  coord_cartesian(xlim=data.xlim, ylim=c(0,225), expand=F, clip="off") +
  xlab(nitrogen.label) + ylab(biomass.label) +
  my.theme
if (show.plots) print(monoculture.plot)


##### 1b. Visualize the total fits #####
# Make a box around the coexistence region
delta.range <- c(-15,5)
coex.range <- competition.fits.wide %>% filter(N1>0 & N2>0) %>%
  .$soil.nitrogen %>% range

# Calculate the difference between competition and best monoculture
competition.delta.plot <- competition.fits.wide %>%
  mutate(soil.nitrogen, biomass=N1+N2) %>%
  left_join(monoculture.fits.wide, 'soil.nitrogen') %>%
  mutate(best.mono=pmax(K1, K2), delta=biomass-best.mono) %>%
  ggplot(aes(x=soil.nitrogen, y=delta)) +
  geom_ribbon(aes(
    ymin=delta, ymax=if_else(N1>0&N2>0, 0, NA)), 
    fill=coex.col, alpha=toy.alpha) +
  geom_line(aes(y=K1-best.mono), size=line.size, linetype="dashed", color=color.poa) +
  geom_line(aes(y=K2-best.mono), size=line.size, linetype="dashed", color=color.agr) +
  geom_line() +
  # Dummy geom to force legend
  geom_line(aes(y=0, linetype=lt), 
    data=data.frame(soil.nitrogen=0:1, lt=c("competition", "monoculture"))) +
  linetype.scale +
  xlim(data.xlim) + # Remove out of x range data
  coord_cartesian(xlim=data.xlim, ylim=delta.range, expand=F) +
  xlab(nitrogen.label) + ylab(delta.label) +
  my.theme
if (show.plots) print(competition.delta.plot)


##### 1c. Visualize fitness ranges for different outcomes #####
fct.fits <- param.fits %>%
  mutate(
    rho=sqrt(alpha12*alpha21/alpha11/alpha22),
    fr=sqrt(alpha21*alpha22/alpha12/alpha11),
    K1=1/alpha11, K2=1/alpha22, yr=K1/K2) %>%
  mutate(
    coex.min=case_when(rho<1/rho ~ rho), 
    coex.max=case_when(rho<1/rho ~ 1/rho),
    toy.min=if_else(yr>1, yr*coex.min, coex.min),
    toy.max=if_else(yr>1, coex.max, yr*coex.max)) %>%
  mutate(
    toy.possible=toy.min<=coex.max & toy.max>=coex.min,
    toy.min=if_else(toy.possible, toy.min, NA),
    toy.max=if_else(toy.possible, toy.max, NA))

fr.range <- c(0.7, 1.3)
fct.plot <- ggplot(fct.fits, aes(x=soil.nitrogen)) +
  geom_ribbon(aes(ymin=fr.range[1], ymax=coex.min, fill="exclusion"), alpha=0) +
  geom_ribbon(aes(ymin=coex.max, ymax=fr.range[2], fill="exclusion"), alpha=0) +
  geom_ribbon(aes(ymin=coex.max, ymax=coex.max, fill="coexistence")) +
  geom_ribbon(aes(ymin=coex.min, ymax=coex.max, fill="coexistence")) +
  geom_ribbon(aes(ymin=toy.min, ymax=toy.max, fill="overyielding")) +
  geom_line(aes(y=coex.min), linetype="dashed") +
  geom_line(aes(y=coex.max), linetype="dashed") +
  geom_line(aes(y=yr*coex.min), linetype="dashed") +
  geom_line(aes(y=fr)) +
  outcome.fill.scale +
  coord_cartesian(xlim=data.xlim, ylim=fr.range, expand=F) +
  geom_hline(yintercept=1, linetype='dotted') +
  xlab(nitrogen.label) + ylab(fr.label) +
  my.theme
if (show.plots) print(fct.plot)


##### Part 4: Combining and saving plots #####
# Ensure output directory exists
dir.create('outputs/')

# Save individual panels
th <- my.theme + bottom.legend.theme
ggsave(paste0(outpath, 'S05_monoculture_fits.pdf'), 
  monoculture.plot + th, 
  device=cairo_pdf, width=out.width, height=out.height)
ggsave(paste0(outpath, '08_a_competition_fits.pdf'),
  competition.plot + th, 
  device=cairo_pdf, width=out.width, height=out.height)
ggsave(paste0(outpath, '08_b_competition_delta_fits.pdf'),
  competition.delta.plot + th, 
  device=cairo_pdf, width=out.width, height=out.height)
ggsave(paste0(outpath, '08_c_fct_outcomes.pdf'),
  fct.plot + th, 
  device=cairo_pdf, width=out.width, height=out.height)