#!/usr/bin/env Rscript
# 04_multispecies/multispecies_helpers.R
# ========================
# Author: Joe Wan
# Description:
# Contains reusable helper functions for the multispecies model analysis.

library(dplyr)
library(tibble)
library(tidyr)
library(stringr)

# Iterates over a range of parameter perturbations (d values),
# applies a transform function to update rhos, frs, and Ks,
# recalculates alphas and B matrix, then stores results.
# -------------------------------------
get.process.df <- function(ds, nfds, Ks, transform.function) {
  df <- NULL
  for (i in seq_along(ds)) {
    factor <- exp(ds[i])
    result <- transform.function(nfds$rhos, nfds$frs, Ks, ds[i])
    new.rhos <- result$rhos
    new.frs <- result$frs
    new.Ks <- result$Ks
    
    # Perform the recalculation of alphas
    new.alphas <- backcalculate.alphas(new.rhos, new.frs, new.Ks)
    B.recalc <- get.B(new.alphas, traits$R0, traits$mus, traits$ups, traits$effs)
    
    # Create updated traits list
    new.traits <- traits
    new.traits$B <- B.recalc
    
    # Add new row to dataframe
    new.row <- tibble_row(
      index = i,
      d = ds[i],
      factor = factor,
      alphas = list(new.alphas),
      rhos = list(new.rhos),
      frs = list(new.frs),
      Ks = list(new.Ks),
      traits = list(new.traits)
    )
    
    df <- if (is.null(df)) new.row else bind_rows(df, new.row)
  }
  return(df)
}

# Species abundance rank helper
add.sad.ranks <- function(df) {
  return(
    df %>%
      group_by(index) %>%
      arrange(desc(N)) %>%
      mutate(rank = row_number())
  )
}

# Pick representative indices
choose.indices <- function(num.sads, ds, pad = 0, pad.orig = 0) {
  if (pad > 0) {
    return(choose.indices(num.sads, ds[(1 + pad):(length(ds) - pad)], pad = 0, pad.orig = pad.orig) + pad)
  }
  
  orig.index <- which(abs(ds) == min(abs(ds)))[1]
  max.index <- length(ds)
  
  orig.pos <- round(num.sads * (orig.index - 1) / (max.index - 1))
  orig.pos <- max(pad.orig + 1, min(num.sads - pad.orig, orig.pos))
  
  spots.before <- orig.pos - 1
  spots.after <- num.sads - orig.pos
  
  skip <- min(floor((max.index - orig.index) / spots.after),
              floor((orig.index - 1) / spots.before))
  
  return(seq(
    orig.index - spots.before * skip,
    orig.index + spots.after * skip,
    by = skip
  ))
}

# Pick example points
get.examples.df <- function(df, rev = FALSE, num.sads = 3, pad = 1, pad.orig = 0) {
  sad.indices <- choose.indices(num.sads, df$d, pad = pad, pad.orig = pad.orig)
  if (rev) sad.indices <- rev(sad.indices)
  
  ref.index <- which(abs(df$d) == min(abs(df$d)))[1]
  
  examples <- data.frame(id = 1:length(sad.indices), index = sad.indices) %>%
    left_join(df, by = "index") %>%
    mutate(is.ref = index == ref.index)
  
  return(examples)
}

# Convert example df into species-level rows
get.example.sp.df <- function(df) {
  return(
    df %>%
      get.species.df() %>%   # NOTE: Requires species_pool_utils.R
      add.sad.ranks() %>%
      ungroup() %>%
      as.data.frame() %>%
      arrange(id) %>%
      mutate(
        facet.title = factor.preserve.order(facet.title),
        facet.title.short = factor.preserve.order(facet.title.short)
      )
  )
}
