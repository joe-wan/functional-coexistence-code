#!/usr/bin/env Rscript
# 02_resource_model/graph_helpers.R
# ================================
# Author: Joe Wan

library(dplyr)
library(cowplot)
library(ggplot2)


# Helper functions for plotting

# Transform MCT's "niche difference" to the "linear" ND (= -log(rho))
to_nd_trans <- function() trans_new("to_nd", function(x) -log(1-x), function(x) 1-exp(-x))

# Custom label: format numbers according to sprintf format if they are in the
# set 'x'
custom_label <- function(format, x, digits=5) {
  return(function(breaks) {
    breaks.rounded <- round(breaks, digits)
    x.rounded <- round(x, digits)
    labels <- sprintf(format, breaks)
    labels <- if_else(breaks.rounded%in%x.rounded, labels, '')
    return(labels)
  })
}