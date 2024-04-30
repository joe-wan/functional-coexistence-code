#!/usr/bin/env Rscript
# 01_fct_simulations/04_additive_partition.R
# ==========================================
# Author: Joe Wan
# Generates figures relating FCT to the additive partition (Figures 3a-c).

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
complementarity.label <- "complementarity"
selection.label <- "selection"


##### 0d. Specific parameters and helpers for varying MCT processes #####
source("process_params.R")
source("process_helpers.R")


##### 1. Define helpers to implement additive partition math #####
# The additive partition is:
# ∆Y = n mean(∆RY) mean(K) + n cov(∆RY, K)
#    = CE                  + SE
# where ∆RY is the difference between observed relative yield (biomass divided by intrinsic yield)
# and expected (simply 1/n in this case).
pop.cov <- function(x, y) sum((x-mean(x)) * (y-mean(y))) / length(x)
calculate.additive.partition <- function(abundance, intrinsic.yield, expected.ry=NULL) {
    if (length(abundance) != length(intrinsic.yield)) stop("Lengths of abundance and intrinsic yield must match")
    if (is.null(expected.ry)) expected.ry <- rep(1/length(abundance), length(abundance))
    # Calculate the additive partition
    delta.ry <- abundance/intrinsic.yield - 1/2
    complementarity <- 2 * mean(delta.ry) * mean(intrinsic.yield)
    selection <- 2 * pop.cov(delta.ry, intrinsic.yield)
    total <- sum(abundance) - sum(expected.ry * intrinsic.yield)
    return(list(complementarity=complementarity, selection=selection, total=total))
}

add.additive.partition <- function(df) {
    result <- df %>%
        as_tibble() %>% rowwise() %>% mutate(
            result=as.data.frame(calculate.additive.partition(c(N1, N2), c(Y1, Y2)))
        ) %>% ungroup() %>% 
        unnest(result) %>%
        select(-biomass) %>%
        as.data.frame() %>% gather(component, biomass, 
            c(complementarity, selection, total))
    return(result)
}


##### 2a. Stabilization (increasing ND) #####
# From process_params.R:
# stab.colors, stab.nd, stab.fr

# Generate simulated data
stab.df <- get.process.df(stab.nd, stab.fr, Y1, Y2) %>%
  mutate(group=fr, 
    # group.name=paste(fr.prefix, fr)
    group.name=paste(group),
    ) %>%
  add.additive.partition()

# Define a helper to modify labels
modify.labels <- function(df, 
        ce.label=complementarity.label,
        se.label=selection.label) {
    df <- df %>% mutate(component=case_when(
        component == 'complementarity' ~ ce.label,
        component == 'selection' ~ se.label,
        TRUE ~ component))
    return(df)
}

# Define a helper to plot the additive partition
plot.additive.partition <- function(df, x.expr, show.total=F, 
        xlab=NULL, show.baseline=F,
        ce.label=complementarity.label,
        se.label=selection.label, 
        parse.facet.labels=F, ...) {
    if (!show.total) df <- df %>% filter(component != 'total')
    df <- df %>% modify.labels(ce.label, se.label)
    result <- df %>% plot.process.df(x.expr, facet=expr(component),
            show.Y1=F, show.oy=F, ylab="effect on biomass",
            ...)
    if (parse.facet.labels) result <- result +
        facet_grid(~component, labeller=label_parsed)
    else result <- result +
        facet_wrap(~component)
    result <- result +
        xlab(xlab) +
        theme_cowplot() + theme(
            # aspect.ratio=ratio,
            strip.background = element_blank(),
            strip.clip="off",
            strip.text = element_text(size=14, margin=margin(b=7)))
    if (show.baseline) result <- result + geom_hline(yintercept=0, linetype=axes.lt)
    return(result)
}

# Plot biomass as we vary niche difference
stab.biomass <- plot.additive.partition(stab.df, 
            expr(nd), legend.name='fitness ratio', 
            xlab=nd.label, add.scale=T,
            colors=stab.colors, show.baseline=T,
            show.total=F) +
    scale_x_continuous(breaks=breaks_width(0.1), expand=c(0,0),
                labels=custom_label('%1.1f', seq(0,1,by=0.2)))
    # scale_x_continuous(breaks=breaks_width(0.1), expand=c(0,0)) # +
    # theme(aspect.ratio=ratio)
if (show.plots) print(stab.biomass)


##### 2b. Equalization (increasing FR) #####
# From process_params.R:
# eq.colors, eq.nd, eq.fr

# Generate simulated data
eq.df <- get.process.df(eq.nd, eq.fr, Y1, Y2) %>%
  mutate(group=nd, 
    # group.name=paste(nd.prefix, nd)
    group.name=paste(group),
    ) %>%
  add.additive.partition()

# Plot biomass as we vary FR
eq.biomass <- plot.additive.partition(eq.df, 
            expr(fr), legend.name='niche difference', 
            xlab=fr.label,
            colors=eq.colors, show.baseline=T,
            show.total=F) +
    # scale_x_continuous(breaks=breaks_width(0.1), expand=c(0,0)) +
      # Add vertical line at optimal FR
    geom_vline(xintercept=sqrt(Y1/Y2), linetype=highlight.lt) +
    scale_x_continuous(breaks=breaks_width(0.1), expand=c(0,0),
                     labels=custom_label('%1.1f', c(0.5, 1, 1.5, 2))) # + # Same as FCT plot's y axis
    # theme(aspect.ratio=ratio)
if (show.plots) print(eq.biomass)


##### 2c. Functional equalization #####
# From process_params.R:
# fu.colors, fu.nd, fu.fr, fu.yr
# delta.yr, yr.examples

# We need to change the simulation strategy
fu.yr.extra <- seq(1, max(fu.yr), by=resolution)
# Ensure that the examples are included in the simulation
fu.yr.extra <- unique(c(fu.yr.extra, yr.examples))

# Generate simulated data
fu.df <- get.process.df(fu.nd, stab.fr, Y1, Y1/fu.yr.extra) %>%
  mutate(yr=Y1/Y2, group=fr, 
    # group.name=paste(fr.prefix, fr),
    group.name=paste(group),
    yr.name=paste(yr.prefix, yr)) %>%
    add.additive.partition()

# Now find points to highlight (special for this plot)
examples.df <- filter(fu.df, 
    yr %in% yr.examples, fr==fu.fr,
    component!='total') %>% modify.labels()

# Plot biomass as we vary fitness
fu.biomass <- plot.additive.partition(fu.df, 
            expr(Y2), legend.name="fitness ratio", 
            colors=stab.colors, show.baseline=T,
            xlab=yield.label,
            show.total=F, 
            # highlight.condition=expr(yr %in% yr.examples)
            ) +
    #  Add green points indicating original highlighted points
    geom_point(data=examples.df, 
        aes(x=Y2, y=biomass, fill=as.factor(yr.name)), 
        shape=21, color="white", size=3) +
    get.scale(colors=fu.colors, legend.name='yield ratio', parse.legend.labels=T, is.fill=T) +
    # scale_x_continuous(expand=c(0,0)) # +
    scale_x_continuous(breaks=breaks_width(0.1), expand=c(0,0),
                labels=custom_label('%1.1f', seq(0.7,1,by=0.2)))
    # theme(aspect.ratio=ratio)
if (show.plots) print(fu.biomass)


##### 3. Export individual plots #####
ti <- theme(legend.position="none")
ggsave(paste0(outpath, '03_a_additive_partition_stabilization.pdf'), 
    stab.biomass + ti,
    device=cairo_pdf, width=out.width, height=out.height)
ggsave(paste0(outpath, '03_b_additive_partition_equalization.pdf'), 
    eq.biomass + ti,
    device=cairo_pdf, width=out.width, height=out.height)
ggsave(paste0(outpath, '03_c_additive_partition_functional_equalization.pdf'), 
    fu.biomass + ti,
    device=cairo_pdf, width=out.width, height=out.height)