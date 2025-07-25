#!/usr/bin/env Rscript
# 04_multispecies/02_supplemental_analyses.R
# ============================================
# Author: Joe Wan
# Generate supplemental figures for multispecies model analyses (Figures S08, S11).



##### 0. Helper functions, etc. #####
# Load utility functions and helpers
library(deSolve)
library(tidyr)
library(dplyr)
library(stringr)
library(tibble)
library(broom)

library(ggplot2)
library(cowplot)
library(cmocean)

source("species_pool_utils.R")
source("multispecies_params.R")
source("multispecies_helpers.R")
source("graph_params.R")
source("graph_helpers.R")

##### 1a. Define general parameters #####
# General parameters are in multispecies_params.R

##### 1b. Generate the system using a defined random seed #####
# Set the seed
set.seed(20240204)
# Generate the traits with the chosen seed
traits <- get.species.traits(metaparams, n.sp)
traits$R0 <- R0
params <- get.model.params(traits, R0)

# Use the original NFDs and Ks from the model defined above
nfds <- get.nfds(traits)
Ks <- get.Ks(traits)

# Output path
outpath <- "outputs/"
if (!str_ends(outpath, '/')) outpath <- paste0(outpath, '/')
dir.create(outpath, showWarnings = FALSE)

##### 2a. Alternate productivity equalization #####
# In the main figure, productivity scaling was around max(K). Here we
# scale around the geometric mean as a sensitivity analysis.
transform.prodeq.alt <- function(rhos, frs, Ks, d) {
  ref <- exp(mean(log(Ks)))
  new.Ks <- exp(log(ref) + (log(Ks) - log(ref)) * exp(d))
  return(list(rhos = rhos, frs = frs, Ks = new.Ks))
}

K.alt.ds <- seq(-0.25, 1.125, length = 76)
cat("Simulating alternate productivity equalization scenarios...\n")

prodeq.alt.df <- get.process.df(K.alt.ds, nfds, Ks, transform.prodeq.alt) %>%
  add.simulated.equilibria(time = 100, verbose = FALSE, refine = TRUE) %>%
  add.equilibrium.summary()

prodeq.alt.examples <- prodeq.alt.df %>%
  get.examples.df(rev = FALSE, num.sads = 3, pad = 5, pad.orig = 0) %>%
  mutate(facet.title.short = paste0(if_else(is.ref, ref.symbol, normal.symbol), " ", id),
         facet.title = paste0("\nKₘₐₓ - Kₘᵢₙ = ", sprintf("%.2f", min.K)))

prodeq.alt.examples.sp <- prodeq.alt.examples %>% get.example.sp.df()

# Plot alternate productivity equalization
prodeq.alt.plot <- plot.process(prodeq.alt.df, prodeq.alt.examples,
                                aes(x = max.K - min.K, y = total.N),
                                color = greens[3], show.min = TRUE, nudge = 0.3) +
  xlab(K.diff.label) + ylab(biomass.label)

prodeq.alt.sads <- plot.sads(prodeq.alt.examples.sp, color = greens[3], max.x = 16)

##### 2b. Histograms of niche, fitness, and productivity #####
# Plot histograms to summarize variation in key traits and pairwise metrics.
plot.histogram <- function(v, bins = 20, xlim = NULL, ylim = NULL) {
  data.frame(x = v) %>%
    ggplot() +
    geom_histogram(aes(x = x), bins = bins, alpha = 0.4) +
    coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
    ylab("# species pairs") +
    cowplot::theme_cowplot()
}

nd.histogram <- plot.histogram(-log(get.off.diags(nfds$rhos, triangular = TRUE)),
                               bins = 13, xlim = trad.to.log(c(0.385, 0.515)), ylim = c(0, 35)) +
  xlab(nd.label)

fd.histogram <- plot.histogram(get.ordered.frs(nfds$frs, Ks), bins = 13,
                               xlim = c(0.8, 1.4), ylim = c(0, 35)) +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  xlab(fd.label)

K.histogram <- plot.histogram(Ks, bins = 9, xlim = c(0, 3), ylim = c(0, 7)) +
  xlab("intrinsic yield, K")

##### 2c. Pairwise plot setup for equalizing scenario #####
# Placeholder for extended analysis of coexistence and overyielding
# (to match structure of original supplemental section).
pairwise <- tibble()
for (i in 1:n.sp) {
  for (j in 1:n.sp) {
    if (Ks[i] > Ks[j]) {
      pairwise <- bind_rows(pairwise,
                             tibble(K1 = Ks[i], K2 = Ks[j],
                                    fr = nfds$frs[i, j], rho = nfds$rhos[i, j]))
    }
  }
}

##### 3. Export plots #####
main.width <- out.width
main.height <- out.height
sads.width <- out.width
sads.height <- out.height * 2 / 3
sad.extra <- scale_x_continuous(breaks = seq(5, 15, by = 5))

# Export histograms
ggsave(paste0(outpath, 'S08_a_nd_histogram.pdf'),
       nd.histogram,
       device = cairo_pdf, width = main.width, height = main.height * 7 / 9)
ggsave(paste0(outpath, 'S08_b_fd_histogram.pdf'),
       fd.histogram,
       device = cairo_pdf, width = main.width, height = main.height * 7 / 9)
ggsave(paste0(outpath, 'S08_c_K_histogram.pdf'),
       K.histogram,
       device = cairo_pdf, width = main.width, height = main.height * 7 / 9)

# Export alternate productivity equalization
ggsave(paste0(outpath, 'S11_prodeq_alt.pdf'),
       prodeq.alt.plot + th,
       device = cairo_pdf, width = main.width, height = main.height)
ggsave(paste0(outpath, 'S11_prodeq_alt_sads.pdf'),
       prodeq.alt.sads + sad.extra + th,
       device = cairo_pdf, width = sads.width, height = sads.height)

