#!/usr/bin/env Rscript
# 04_multispecies/03_misc_analyses.R
# ===================================
# Author: Joe Wan
# Additional multispecies analyses for supplement: time series for different 
# species pools and invasion growth rate analysis (Figures 5b, S04, and S12).

##### 0. Helper functions, etc. #####
# Load required libraries and sources
library(deSolve)
library(tidyr)
library(dplyr)
library(stringr)
library(tibble)
library(broom)

library(ggplot2)
library(cowplot)
library(ggpubr)

source("species_pool_utils.R")
source("multispecies_params.R")
source("multispecies_helpers.R")
source("graph_params.R")
source("graph_helpers.R")

##### 1a. Define general parameters #####
# General parameters are in multispecies_params.R

##### 1b. Generate system using defined random seed #####
# Set the seed
set.seed(20240204)
traits <- get.species.traits(metaparams, n.sp)
traits$R0 <- R0
params <- get.model.params(traits, R0)

# Output path
outpath <- "outputs/"
if (!str_ends(outpath, '/')) outpath <- paste0(outpath, '/')
dir.create(outpath, showWarnings = FALSE)

##### 2a. Example time series for reference system #####
# Numbers of species to simulate
example.ns <- c(1,2,20)
# Time series length
example.timeseries.length <- 100

# Simulate the time series
example.ts <- tibble(dummy = TRUE) %>%
  # Add traits
  mutate(traits = list(traits)) %>%
  # Join with all numbers of species to simulate
  right_join(data.frame(n.sp = example.ns, dummy = TRUE), by = "dummy") %>%
  select(-dummy) %>%
  # Subset traits for the current number of species
  rowwise() %>%
  mutate(traits = list(subset.system(traits, 1:n.sp))) %>%
  # Simulate time series for each subset
  add.simulated.timeseries(seq(0, example.timeseries.length, by = 0.1)) %>%
  # Unnest the nested list column "timeseries" into rows
  unnest(timeseries) %>%
  # Add the resource use efficiency for plotting
  left_join(
    data.frame(species = 1:n.sp, eff = traits$effs),
    by = "species"
  )

# Calculate resource level, etc.
example.ts.R <- example.ts %>% 
    group_by(time, n.sp) %>%
    summarize(consumed=sum(N / eff), Ntot=sum(N)) %>%
    mutate(R=traits$R0-consumed) %>%
    ungroup

# Grapical options
ts.line.width <- 0.5
 # Theme for all panels
tts <- th +
    # Increase margin on the right
    theme(plot.margin=margin(t=7, b=7, l=7, r=14))

# Plot the time series by looping through number of species
example.ts.plots <- list()
example.ts.R.plots <- list()
example.ts.composites <- list()
for (n in unique(example.ts$n.sp)) {
    # Subset needed data
    cur.ts <- filter(example.ts, n.sp==n)
    cur.R <- filter(example.ts.R, n.sp==n)
    max.n <- max(cur.ts$N)

    # Plot the time series 
    example.ts.plots[[n]] <- cur.ts %>%
        ggplot(aes(x=time)) +
        geom_line(aes(y=N, color=as.factor(species)), linewidth=ts.line.width) +
        expand_limits(y=c(0, 1, max.n + 0.01)) +
        scale_y_continuous(breaks=seq(0,10,by=0.5)) +
        coord_cartesian(expand=F) +
        ylab(abundance.label) +
        theme_cowplot() + tts +
        theme(legend.position="none")

    # Plot resource level
    example.ts.R.plots[[n]] <- cur.R %>%
        ggplot(aes(x=time)) +
        geom_line(aes(y=R), linewidth=ts.line.width) +
        expand_limits(y=c(0, traits$R0)) +
        coord_cartesian(expand=F) +
        ylab(resource.label) +
        theme_cowplot() + tts +
        theme(legend.position="none")
    
    # Add more lines to plots depending on number of species
    if (n > 5) example.ts.R.plots[[n]] <- example.ts.R.plots[[n]] +
        geom_line(aes(y=consumed), 
        linewidth=ts.line.width, linetype="dashed")
    # For fewer, we can show each species' consumption
    else example.ts.R.plots[[n]] <- example.ts.R.plots[[n]] +
        geom_line(aes(y=N/eff, color=as.factor(species)), 
        linewidth=ts.line.width, linetype="dashed",
            data=cur.ts)
    
    # Save
    ggsave(paste0(outpath, 'S04_example_timeseries_n=', n, '.pdf'), 
        example.ts.plots[[n]],
        device=cairo_pdf, width=main.width, height=main.height*0.75)
    ggsave(paste0(outpath, 'S04_example_timeseries_n=', n, '_R.pdf'), 
        example.ts.R.plots[[n]],
        device=cairo_pdf, width=main.width, height=main.height*0.75)
}

# S04_example_timeseries_n=20.pdf and S04_example_timeseries_n=20_R are also Figure 5 in the main text
# Copy them:
file.copy(
    paste0(outpath, 'S04_example_timeseries_n=20.pdf'), 
    paste0(outpath, '05_b_example_timeseries.pdf'))
file.copy(
    paste0(outpath, 'S04_example_timeseries_n=20_R.pdf'), 
    paste0(outpath, '05_b_example_timeseries_R.pdf'))

##### 2b. Invasion growth rate analysis #####
get.n0s <- function(n.sp, i, n0=1e-5, focal=0) {
    # Get the initial abundances for the invasion
    n0s <- rep(n0, n.sp)
    n0s[i] <- 0
    return(n0s)
}

del.sp <- function(equilibrium, i) {
    # Get the equilibrium abundance of species i
    result <- equilibrium
    # Remove the species
    result[i] <- 0
    return(result)
}

ref.df <- data.frame(x=NA) %>%
    rowwise %>%
    mutate(traits=list(traits)) %>%
    mutate(equilibrium=list(
        get.equilibrium(resource.model, get.model.params(traits), resource.growth)
    )) %>%
    select(-x)
ref.eq <- ref.df$equilibrium[[1]]

get.subset <- function(n.sp, i) {
    result <- 1:n.sp
    result <- result[-i]
    return(result)
}

expand.equilibrium <- function(equilibrium, i) {
    # Expand the equilibrium to include the invader
    if (i == 1) head <- c()
    else head <- equilibrium[1:(i-1)]
    if (i == length(equilibrium)+1) tail <- c()
    else tail <- equilibrium[i:length(equilibrium)]
    return(c(head, 0, tail))
}

invasion.df <- data.frame(i=1:n.sp) %>%
    rowwise %>%
    mutate(sub.traits=list(
        subset.system(traits, get.subset(n.sp, i))
    )) %>% mutate(sub.equilibrium=list(
        get.equilibrium(resource.model, get.model.params(sub.traits), resource.growth)
    )) %>% mutate(equilibrium=list(
        expand.equilibrium(sub.equilibrium, i)
    )) %>% mutate(
        igr=resource.growth(equilibrium, get.model.params(traits))[i],
        r=resource.growth(rep(0, n.sp), get.model.params(traits))[i],
        K=get.Ks(traits)[i]
    ) %>% ungroup

full.df <- ref.df %>%
    crossing(i=1:n.sp) %>%
    rowwise %>%
    mutate(N=equilibrium[i]) %>% 
    ungroup %>%
    left_join(invasion.df, by="i")

# Estimate the slope of the N ~ FK relationship
FK.slope <- full.df %>% filter(igr>0) %>%
    transmute(FK=K*igr/r, N) %>%
    lm(N ~ 0 + FK, data=.) %>%
    tidy() %>%
    filter(term=='FK') %>% pull(estimate)
cat("Slope of N ~ FK relationship: ", FK.slope, "\n")

FK.plot <- ggplot(full.df) +
    geom_point(aes(x=igr/r * K, y=N)) +
    geom_line(aes(x=x, y=pmax(0, FK.slope*x)), 
        data=data.frame(x=seq(-0.05,0.1, len=1000)),
        linetype='dashed') +
    # geom_hline(yintercept=0) + 
    # geom_vline(xintercept=0) +
    xlab(expression(paste('fitness' %*% 'intrinsic yield, ', F[i]~K[i]))) + 
    ylab(expression(paste('equilbrium abundance, ', N[i]))) +
    scale_x_continuous(expand=c(0,0)) +
    cowplot::theme_cowplot() +
    theme(legend.position="none")

# Save the plot
ggsave(paste0(outpath, 'S12_invasion_growth.pdf'), 
    FK.plot + th, 
    device=cairo_pdf, width=main.width, height=main.height)