#!/usr/bin/env Rscript
# 04_multispecies/01_multispecies.R
# ==================================
# Author: Joe Wan
# Purpose: Generate primary multispecies model figures (Figures 06a-c and related outputs).

library(deSolve)

library(tidyr)
library(dplyr)
library(stringr)
library(tibble)
library(broom)

library(ggplot2)
library(cowplot)
library(cmocean)

##### 0. Helper functions, etc. #####
# Load functions from utility files
source("species_pool_utils.R")
source("multispecies_params.R")
source("multispecies_helpers.R")
source("graph_helpers.R")
source("graph_params.R")

##### 1. Generate the system using a defined random seed #####
set.seed(20240204)
# General parameters and traits
traits <- get.species.traits(metaparams, n.sp)
traits$R0 <- R0
params <- get.model.params(traits, R0)
# Use the original NFDs and Ks from the model defined above
nfds <- get.nfds(traits)
Ks <- get.Ks(traits)

##### 2a. Set up back-calculation visualization #####
# We'll vary NFD by changing the pairwise metrics, then back-calculate model parameters.

##### 2b. Simulate scenarios #####
# Panel (a): Stabilizing 
# ----------------------
# Stabilizing i.e. decreasing niche overlap corresponds to increasing rho. 
# To do so, we will add some factor d to log(rho_ij), equivalent to multiplying by exp(d).

# Define transformation functions
transform.stabilizing <- function(rhos, frs, Ks, d) {
    new.rhos <- exp(log(rhos) + d)
    diag(new.rhos) <- 1
    return(list(rhos = new.rhos, frs = frs, Ks = Ks))
}

transform.equalizing <- function(rhos, frs, Ks, d) {
    n.sp <- length(Ks)
    signs <- sign(diag(Ks) %*% matrix(1, n.sp, n.sp) - matrix(1, n.sp, n.sp) %*% diag(Ks))
    new.frs <- exp(log(frs) + signs * d)
    return(list(rhos = rhos, frs = new.frs, Ks = Ks))
}

transform.prodeq <- function(rhos, frs, Ks, d) {
    ref <- max(Ks)
    new.Ks <- exp(log(ref) + (log(Ks) - log(ref)) * exp(d))
    return(list(rhos = rhos, frs = frs, Ks = new.Ks))
}

rho.ds <- seq(-0.15, 0.6, by = 0.01)
cat("Simulating stabilizing scenarios...\n")
stabilizing.df <- get.process.df(rho.ds, nfds, Ks, transform.stabilizing) %>%
    # Simulate to equilibrium and add summary stats
    add.simulated.equilibria(time = 100, verbose = FALSE, refine = TRUE) %>% 
    add.equilibrium.summary()

stabilizing.examples <- stabilizing.df %>%
    get.examples.df(rev = TRUE, num.sads = 3, pad = 1, pad.orig = 0) %>%
    mutate(facet.title.short = paste0(if_else(is.ref, ref.symbol, normal.symbol), " ", id),
           facet.title = paste0(facet.title.short, "\n", "-log ρ = ", sprintf("%.2f", -log(median.rho))))

stabilizing.examples.sp <- stabilizing.examples %>%
    get.example.sp.df()

# Panel (b): (De)equalization 
# ---------------------------
# We want to increase fitness in favor of the more productive species.
# To do so, we will add some factor sgn(K_i-K_j)*d to log(f_i/f_j), equivalent to multiplying/dividing by exp(d).

fr.ds <- seq(-0.9, 0.6, by = 0.02)
cat("Simulating equalizing scenarios...\n")
equalizing.df <- get.process.df(fr.ds, nfds, Ks, transform.equalizing) %>%
    add.simulated.equilibria(time = 100, verbose = FALSE, refine = TRUE) %>%
    add.equilibrium.summary()

equalizing.examples <- equalizing.df %>%
    get.examples.df(rev = FALSE, num.sads = 3, pad = 1, pad.orig = 1) %>%
    mutate(facet.title.short = paste0(if_else(is.ref, ref.symbol, normal.symbol), " ", id),
           facet.title = paste0(facet.title.short, "\n", "f₁/f₂ = ", sprintf("%.2f", median.ordered.fr)))

equalizing.examples.sp <- equalizing.examples %>%
    get.example.sp.df()

# Panel (c): Functional (productivity) equalization
# -------------------------------------------------
# Increase or decrease variation in productivity by scaling Ks.
# We scale the distribution of Ks around the most productive species.

K.ds <- seq(-1.5, 1.5, length = 76)
cat("Simulating productivity equalization scenarios...\n")
prodeq.df <- get.process.df(K.ds, nfds, Ks, transform.prodeq) %>%
    add.simulated.equilibria(time = 100, verbose = FALSE, refine = TRUE) %>%
    add.equilibrium.summary()

prodeq.examples <- prodeq.df %>%
    get.examples.df(rev = TRUE, num.sads = 3, pad = 5, pad.orig = 1) %>%
    mutate(facet.title.short = paste0(if_else(is.ref, ref.symbol, normal.symbol), " ", id),
           facet.title = paste0(facet.title.short, "\n", "Kₘᵢₙ = ", sprintf("%.2f", min.K)))

prodeq.examples.sp <- prodeq.examples %>%
    get.example.sp.df()

##### 3a. Set up plotting #####
outpath <- "outputs/"
if (!str_ends(outpath, '/')) outpath <- paste0(outpath, '/')
dir.create(outpath, showWarnings = FALSE)

##### 3b. Plot the scenarios #####
# Panel (a): Stabilizing
# ----------------------
stabilizing.y <- c(2.85, 4.15)
stabilizing.plot <- stabilizing.df %>%
    filter(median.rho > 0.52) %>%
    plot.process(stabilizing.examples, aes(x = -log(median.rho), y = total.N),
                 color = blues[2], ref.color = NULL,
                 show.min = FALSE, expand.y = stabilizing.y, nudge = nudge.factor * 1) +
    xlab(median.nd.label) + ylab(biomass.label)
stabilizing.sads <- plot.sads(stabilizing.examples.sp, color = blues[2], max.x = 16)

# Panel (b): Equalizing
# ---------------------
optimal.fitness <- equalizing.df %>% arrange(desc(total.N)) %>% head(1) %>% pull(median.ordered.fr)
equalizing.plot <- plot.process(equalizing.df, equalizing.examples,
                                aes(x = median.ordered.fr, y = total.N),
                                color = oranges[3], ref.color = NULL,
                                show.min = TRUE, ylim = c(1.25, 4.25), nudge = nudge.factor * -2.25) +
    geom_vline(xintercept = optimal.fitness, linetype = 'dashed') +
    scale_x_continuous(breaks = seq(0.2, 2.2, by = 0.4)) +
    xlab(median.fr.label) + ylab(biomass.label)
equalizing.sads <- plot.sads(equalizing.examples.sp, color = oranges[3], max.x = 16)

# Panel (c): Productivity equalization
# ------------------------------------
prodeq.plot <- plot.process(prodeq.df, prodeq.examples,
                            aes(x = min.K, y = total.N),
                            color = greens[3], ref.color = NULL,
                            show.min = TRUE, ylim = c(1.95, 4.95), nudge = nudge.factor * -2.5) +
    xlab(K.label) + ylab(biomass.label)
prodeq.sads <- plot.sads(prodeq.examples.sp, color = greens[3], max.x = 16, max.y = 1) +
    scale_y_continuous(breaks = seq(0, 1, by = 1))

##### 3c. Export single panels #####
# All the panels have the same height, so we can keep them in proportion
# by exporting with the same dimensions
ggsave(paste0(outpath, '06_a_stabilizing.pdf'),
       stabilizing.plot + th, device = cairo_pdf, width = out.width, height = out.height)
ggsave(paste0(outpath, '06_a_stabilizing_sads.pdf'),
       stabilizing.sads + sad.extra + th, device = cairo_pdf, width = out.width, height = out.height * 2 / 3)

ggsave(paste0(outpath, '06_b_equalizing.pdf'),
       equalizing.plot + th, device = cairo_pdf, width = out.width, height = out.height)
ggsave(paste0(outpath, '06_b_equalizing_sads.pdf'),
       equalizing.sads + sad.extra + th, device = cairo_pdf, width = out.width, height = out.height * 2 / 3)

ggsave(paste0(outpath, '06_c_prodeq.pdf'),
       prodeq.plot + th, device = cairo_pdf, width = out.width, height = out.height)
ggsave(paste0(outpath, '06_c_prodeq_sads.pdf'),
       prodeq.sads + sad.extra + th, device = cairo_pdf, width = out.width, height = out.height * 2 / 3)


##### 3. Plot pairwise NFD for equalizing scenario (Figure S10) #####
# Calculate pairwise metrics
pairwise <- tibble()
for (i in 1:n.sp) {
    for (j in 1:n.sp) {
        if (Ks[i]>Ks[j]) {
            new.row <- tibble(
                K1=Ks[i], K2=Ks[j], 
                fr=nfds$frs[i,j], 
                rho=nfds$rhos[i,j])
            pairwise <- bind_rows(new.row, pairwise)
        }
    }
}

# Calculate fitness ratios for key events
eq.coex.fr.range <- equalizing.df %>% filter(richness>1) %>% 
    pull(median.ordered.fr) %>% range
eq.toy.fr.range <- equalizing.df %>% filter(total.N>max.K) %>% 
    pull(median.ordered.fr) %>% range
eq.ref.fr <- equalizing.df %>% filter(d==0) %>%
    pull(median.ordered.fr)
eq.opt.fr <- equalizing.df %>% arrange(desc(total.N)) %>% head(1) %>% 
    pull(median.ordered.fr)

# Calculate regions for (excess) NFD conditions
regions <- data.frame(rho=seq(0.5,1.1,by=0.001))
frs <- unique(c(eq.ref.fr, eq.opt.fr, eq.coex.fr.range, eq.toy.fr.range))
pairwise.examples <- data.frame(median.ordered.fr=frs) %>%
    arrange(median.ordered.fr) %>%
    mutate(id=row_number()) %>%
    left_join(equalizing.df, by="median.ordered.fr") %>%
    mutate(is.ref=d==0,
        facet.title=paste0("atop('", 
            paste(if_else(is.ref, ref.symbol, normal.symbol), id), 
            "', f[1]/f[2] == '",  sprintf("%.2f", median.ordered.fr), "')")) %>%
    mutate(facet.title=factor(facet.title,      levels=facet.title))
conditions <- c("transgressive overyielding", "coexistence")

# Plot faceted plot
pairwise.fct.plot <- pairwise %>%
    mutate(shift=sqrt(K1/K2)) %>%
    expand_grid(pairwise.examples) %>%
    mutate(factor=median.ordered.fr/eq.ref.fr) %>%
    expand_grid(condition=conditions) %>%
    # Coexistence compares -log rho and -log f1/f2,
    # while TOY divides these by sqrt(K1/K2)
    mutate(shift=if_else(
        condition=="transgressive overyielding", 
        shift, 1)) %>%
    mutate(x=log(1/rho/shift)) %>%
ggplot() +
    # Coexistence or TOY requirement zone
    geom_ribbon(aes(x=x, ymin=x, ymax=-x, 
            fill=condition), 
        data=expand_grid(x=seq(0,0.8,by=0.01), condition=conditions)) +
    geom_abline(intercept=0, slope=c(1,-1),             linetype='solid') +
    geom_hline(yintercept=0, linetype="dotted") +
    geom_point(aes(x=x, y=log(factor*fr/shift)),
        alpha=1, size=0.1) +
    coord_equal(ratio=1, 
        expand=F) +
    facet_grid(paste0("'", condition, "'")~facet.title, labeller=label_parsed) +   
    xlab(expression(paste(
        "(excess) niche difference, ", 
        - "log" ~ rho - Delta))) +
    ylab(expression(atop(
        "(excess) fitness difference, ", 
        "log" ~ f[1]/f[2] - Delta))) +
    scale_x_continuous(breaks=seq(0,0.5,by=0.5)) +
    scale_fill_manual(values=c(coexistence=coex.col, 
        `transgressive overyielding`=greens[2])) +
    guides(fill=guide_none()) +
    theme_cowplot() +
    theme(strip.background=element_blank(), strip.clip="off")
if (show.plots) print(pairwise.fct.plot)

# Re-plot the process plot with extra examples
pairwise.process.plot <- plot.process(equalizing.df, pairwise.examples, 
    aes(x=median.ordered.fr, y=total.N), color=oranges[3], ref.color=NULL,
    show.min=T, ylim=c(1.25, 4.25), nudge=nudge.factor*-1.5) +
    # Make sure 1 appears on the x axis
    scale_x_continuous(breaks=seq(0.2, 2.2, by=0.4)) +
    xlab(median.fr.label) + ylab(biomass.label) # fᵢ/fⱼ
if (show.plots) print(pairwise.process.plot)

# Save (Figure S10)
# -----------------
ggsave(paste0(outpath, 'S10_pairwise_NFD.pdf'),
    pairwise.fct.plot + th, device = cairo_pdf, width = out.width*1.5, height = out.height*1.5)
ggsave(paste0(outpath, 'S10_pairwise_process.pdf'),
    pairwise.process.plot + th, device = cairo_pdf, width = out.width*1.5, height = out.height*1.5)