#!/usr/bin/env Rscript
# 02_resource_model/01_resource_model.R
# =======================================
# Author: Joe Wan
# Generates the theoretical result figure for resource model (Figure 4b).

# Data manipulation libraries
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)

# Plotting libraries
library(ggplot2)
library(cowplot)
library(scales)


##### 0a. Define graph options #####
source('graph_params.R')


##### 0b. User-defined parameters #####
# Range of values of R
R.i <- 200
R.f <- 1000

# Actual value of niche overlap
rho <- 0.8 # Coexistence boundary: 0.8 < f1/f2 < 1.25
# Asymptotic value of fitness ratio f1/f2 when R is large
FR.asy <- 1.5
# Ratio of threshold fitness for OY to actual fitness
thresh.ratio <- 1/0.9 # = a12/a11
# Fitness ratio at minimum R
FR.i <- 1/2 # Should be < rho
# Yield of species 2 at minimum R
K2.i <- 22
# Rstar of species 2
Rstar2 <- 2 # Less than R.i


##### 1. Calculate parameter values #####
# Solve for Rstars and a22
q.i <- FR.i/FR.asy
Rstar1 <- (1-q.i)*R.i + q.i*Rstar2
a22 <- (R.i-Rstar2)/K2.i
# Solve for remaining alphas
S1.asy <- rho/FR.asy # = a12/a22
S2.asy <- rho*FR.asy # = a21/a11
a12 <- S1.asy*a22
a11 <- a12/thresh.ratio
a21 <- S2.asy*a11


##### 1a. Export to CSV #####
# Parameter values
params <- c(a11=a11, a12=a12, a21=a21, a22=a22, Rstar1=Rstar1, Rstar2=Rstar2)

# Biological parameters corresponding to the chosen parameter values
params.bio <- c(params, eff1=0.9, eff2=0.5, mu1=5, mu2=1)
params.bio <- c(params.bio, with(as.list(params.bio), c(
    v1=mu1/(eff1*Rstar1), v2=mu2/(eff2*Rstar2),
    beta11=(a11-eff1^-1)/Rstar1, beta22=(a22-eff2^-1)/Rstar2,
    beta12=(a12-eff2^-1)/Rstar1, beta21=(a21-eff1^-1)/Rstar2
)))


##### 2. Plot coexistence and overyielding boundaries #####
# Values of resource on x axis
Rs <- seq(R.i, R.f, by=0.01)

# Calculate the fitness ratio and overyielding conditions
gradient.data <- data.frame(R=Rs) %>%
    mutate(rho=sqrt(a12*a21/(a11*a22)),
        q=(R-Rstar1)/(R-Rstar2),
        FR=q*sqrt(a21*a22/(a12*a11)),
        YR=q*a22/a11)

# Pick out some examples
gradient.examples <- data.frame(
    # R=c(270, 360, 470, 630, 850)) %>%
    # R=c(270, 470, 850)) %>%
    R=c(270, 380, 580, 850)) %>%
    mutate(id=row_number()) %>%
    left_join(gradient.data, by='R')

# Plot the effect of changing resource
resource.effect <- ggplot(gradient.data, aes(x=R)) +
    geom_ribbon(aes(
            ymin=rho, 
            ymax=1/rho), 
        fill=coex.col) +
    # There is a minor edge case when YR<1 (i.e. the other species is actually superior)
    geom_ribbon(aes(
            ymin=pmax(pmin(YR*rho,1/rho), rho), 
            ymax=if_else(YR>1, 1/rho, NA)),
        fill=greens[2]) +
    # Fill in a tiny spurious gap in the green area
    geom_ribbon(aes(
            ymin=if_else(abs(YR-1)<0.005, rho, NA),
            ymax=1/rho),
        fill=greens[2]) +
    geom_ribbon(aes(
            ymin=if_else(YR<=1, rho, NA),
            ymax=pmin(pmax(YR*1/rho, rho), 1/rho)),
        fill=greens[2], color=greens[2]) +
    geom_hline(yintercept=1, linetype='dotted') +
    geom_line(aes(y=FR)) +
    geom_line(aes(y=rho*YR), linetype='dashed') +
    geom_line(aes(y=1/rho), linetype='dashed') +
    geom_line(aes(y=rho), linetype='dashed') +
    geom_point(aes(y=FR), data=gradient.examples) +
    # Add text labels below and horizontally centered
    geom_text(aes(y=FR, label=id), nudge_y=-0.05, 
        hjust=0.5, data=gradient.examples) +
    # geom_line(aes(y=sqrt(YR)), linetype='dashed') +
    # scale_y_log10() +
    # Add untransformed x-axis
    scale_x_continuous(limits=c(225, 900), expand=c(0,0)) +
    scale_y_continuous(limits=c(0.6, 1.45), expand=c(0,0), breaks=seq(0,2,by=0.2)) +
    xlab(resource.label) + ylab(fr.label) +
    theme_cowplot()
if (show.plots) print(resource.effect)


##### 3. Plot a series of FCT space diagrams #####
# Define the range of values
min.FD <- -0.7; max.FD <- 0.7
min.ND <- -0.2; max.ND <- 1.2
ratio <- 1
# Helper to transform MCT's "niche difference" to the "linear" ND (= -log(rho))
to_nd_trans <- function() trans_new("to_nd", function(x) -log(1-x), function(x) 1-exp(-x))

# Define the niche/fitness space
mct.data <- data.frame(ND=seq(min.ND, max.ND, by=0.001)) %>%
    mutate(rho=exp(-ND), ND=1-exp(-ND))

# Plot the FCT space for different examples
fct.strip <- gradient.examples %>%
    mutate(dummy=T) %>%
    left_join(mutate(mct.data, dummy=T), by='dummy',
        relationship='many-to-many',
        suffix=c('.point', '')) %>%
ggplot(aes(x=ND)) +
    geom_ribbon(aes(ymin=if_else(ND<0, 1/rho, NA), 
        ymax=rho), fill=ass.col) +
    geom_ribbon(aes(ymin=if_else(ND>0, rho, NA), 
        ymax=1/rho), fill=coex.col) +
    geom_ribbon(
        aes(ymin=if_else((rho*YR<1/rho) & (rho*YR>rho), rho*YR, NA), 
        ymax=1/rho), fill=greens[2]) +
    geom_vline(aes(xintercept=1-rho), linetype='dashed',
        data=gradient.examples) +
    geom_point(aes(x=1-rho, y=FR), data=gradient.examples) +
    geom_text(aes(x=1-rho, y=FR, label=id), 
        nudge_x=-0.075, vjust=0.5,
        data=gradient.examples) +
    coord_trans(x="to_nd", y="log",
                ylim=c(0.7,1.4), 
                xlim=1-exp(-log(1.4/0.7)/ratio*c(-1/4,3/4)),
                expand=F) +
    scale_x_continuous(breaks=c(-0.3,0,0.3)) +
    facet_wrap(~id, nrow=1, scales='free') +
    xlab(nd.label) +
    theme_cowplot() + 
    theme(aspect.ratio=ratio, 
        strip.background=element_blank(),
        strip.text=element_blank(),
        axis.title.y = element_blank())
if (show.plots) print(fct.strip)

# Save individual plots
# ---------------------
ggsave(paste0(outpath, '04_b_resource_effect.pdf'), resource.effect + th, 
       device=cairo_pdf, width=out.width*4/3, height=out.height*1)
ggsave(paste0(outpath, '04_b_fct_strip.pdf'), fct.strip + th, 
       device=cairo_pdf, width=out.width*4/3, height=out.height*2/3)