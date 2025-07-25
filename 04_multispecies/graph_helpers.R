#!/usr/bin/env Rscript
# 04_multispecies/graph_helpers.R
# ===============================
# Author: Joe Wan
# Helper functions for plotting multispecies analysis.

library(ggplot2)
library(cowplot)
library(dplyr)

# Create main process plots (biomass vs parameter)
# ----------------------------------------------------
plot.process <- function(df, example.df, aesthetic, color="black", ref.color=NULL,
                         show.min=TRUE, show.max=TRUE, show.ribbon=TRUE,
                         nudge=0.125, expand.x=NULL, expand.y=NULL, 
                         xlim=NULL, ylim=NULL) {
  if (is.null(ref.color)) ref.color <- color
  
  result <- ggplot(df, aesthetic) +
    geom_path(color=color, linewidth=0.75)
  
  # Add dotted lines for min and max yield
  if (show.min) result <- result + geom_line(aes(y=min.K), linetype='dotted')
  if (show.max) result <- result + geom_line(aes(y=max.K), linetype='dotted')
  
  # Add shaded ribbon for overyielding
  if (show.ribbon) {
    result <- result +
      geom_ribbon(aes(ymin=max.K, ymax=ifelse(total.N > max.K, total.N, NA)), 
                  fill=greens[2], alpha=0.1)
  }
  
  # Add labels and points for examples
  result <- result +
    geom_text(aes(label=id), data=example.df, nudge_y=nudge, hjust=0.5) +
    geom_point(data=filter(example.df, !is.ref), color=color) +
    geom_text(data=filter(example.df, is.ref), color=ref.color, label='★', size=6)
  
  # Axis and limits
  if (!is.null(expand.x)) result <- result + expand_limits(x=expand.x)
  if (!is.null(expand.y)) result <- result + expand_limits(y=expand.y)
  
  result <- result +
    coord_cartesian(expand=FALSE, xlim=xlim, ylim=ylim) +
    theme_cowplot()
  
  return(result)
}


# Create species abundance distribution bar plots
# -----------------------------------------------
plot.sads <- function(examples.sp, color="black", alpha=1, 
                      max.x=NULL, max.y=NULL, expand.y=NULL, shared.y=TRUE) {
  if (is.null(max.x)) max.x <- max(examples.sp$richness)
  if (is.null(max.y)) max.y <- max(examples.sp$N)
  
  yl <- c(0, max.y)
  if (!shared.y) yl <- NULL
  
  result <- ggplot(examples.sp) +
    geom_bar(aes(x=rank, y=N), stat="identity", fill=color, alpha=alpha) +
    facet_wrap(~facet.title, nrow=1, scales="free") +
    theme_cowplot() +
    coord_cartesian(xlim=c(0, max.x) + 0.5, ylim=yl, expand=FALSE) +
    xlab("species rank") + ylab("abundance") +
    theme(strip.background=element_blank())
  
  # Adjust y limits if necessary
  if (is.null(expand.y) && !shared.y) expand.y <- c(0)
  else if (!shared.y) expand.y <- c(expand.y, 0)
  
  if (!is.null(expand.y)) result <- result + expand_limits(y=expand.y)
  
  return(result)
}


# Generic histogram for distributions
# -----------------------------------
plot.histogram <- function(v, bins=20, xlim=NULL, ylim=NULL, xlab_text="Value") {
  plt <- data.frame(x=v) %>%
    ggplot() +
    geom_histogram(aes(x=x), bins=bins, alpha=0.4) +
    coord_cartesian(xlim=xlim, ylim=ylim, expand=FALSE) +
    xlab(xlab_text) +
    ylab("# species pairs") +
    theme_cowplot()
  
  return(plt)
}

# Mini-helpers for coordinate conversion, etc.
trad.to.log <- function(x) -log(1 - x)
d <- function(a) a[2] - a[1]