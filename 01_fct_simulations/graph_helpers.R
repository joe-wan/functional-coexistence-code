# 01_fct_simulations/graph_helpers.R
# ==================================
# Author: Joe Wan
# Helper functions for plotting MCT and FCT figures.

library(dplyr)
library(cowplot)
library(ggplot2)


##### Helper functions for plotting #####
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

# Create label function to include r[1] and r[2], with resolution inc and
# print format format
get.label <- function(r, inc, format) {
  if (is.null(r) || is.null(inc) || is.null(format)) return(NULL)
  rounded.r <- inc*c(floor(r[1]/inc), ceiling(r[2]/inc))
  breaks <- seq(rounded.r[1], rounded.r[2], by=inc)
  return(custom_label(format, breaks))
}

# Return rectangle to fade out the plot
fade.box <- function(xlim, ylim, alpha=0.85) geom_rect(
  xmin=xlim[1], xmax=xlim[2], ymin=ylim[1], ymax=ylim[2],
  fill='white', alpha=alpha, data=data.frame(x=1))

##### MCT functions #####
# Get MCT data frame
get.mct.df <- function(nds) {
  result <- data.frame(ND=nds) %>%
  mutate(pers.1=-ND, pers.2=ND,
         coex.min=case_when(ND>=0 ~ pers.1), coex.max=case_when(ND>=0 ~ pers.2),
         ass.min=case_when(ND<=0 ~ pers.2), ass.max=case_when(ND<=0 ~ pers.1))
}

# Visualize MCT data frame
plot.mct.df <- function(df, simple.axes=T,
    xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL,
    zero.linetype=NULL) {
  result <- ggplot(df) +
    # Plot coexistence and ASS regions:
    geom_ribbon(aes(x=ND, ymin=coex.min, ymax=coex.max), fill=coex.col) +
    geom_ribbon(aes(x=ND, ymin=ass.min, ymax=ass.max), fill=ass.col) +
    # Plot boundaries (persistence conditions)
    geom_line(aes(x=ND, y=pers.1)) + geom_line(aes(x=ND, y=pers.2)) +
    # Add labels
    ggplot2::xlab(xlab) + ggplot2::ylab(ylab) +
    theme_cowplot()
  if (!is.null(zero.linetype)) result <- result +
    geom_hline(aes(yintercept=0), linetype=axes.lt) + 
    geom_vline(aes(xintercept=0), linetype=axes.lt)
  if (simple.axes) result <- result +
    # Only show '0' on axes, since this is a conceptual figure
    scale_x_continuous(breaks=breaks_pretty(1), expand=c(0,0)) +
    scale_y_continuous(breaks=breaks_pretty(1), expand=c(0,0))
  if (!is.null(xlim) || !is.null(ylim)) result <- result +
    coord_fixed(xlim=xlim, ylim=ylim, expand=F)
  return(result)
}


##### FCT functions #####
# Go from log back to normal NFD
undo.log.version <- function(df) df %>% 
  mutate(across(-ND, ~ exp(.x))) %>% 
  mutate(ND=1-exp(-ND))
  
# Go from normal NFD to log
get.log.version <- function(df) df %>% 
  mutate(across(-ND, ~ log(.x))) %>% 
  mutate(ND=-log(1-ND))

# Get FCT data frame
get.fct.df <- function(nds, Y1=NULL, Y2=NULL, log.version=F) {
  if (any(is.null(c(Y1, Y2)))) stop('Y1 and Y2 must be specified')
  if (is.null(log.version)) stop('log.version must be specified')
  result <- get.mct.df(nds)  %>%
    mutate(optimal=log(Y1/Y2)/2,
          toy.min=case_when(ND>=optimal~log(Y1/Y2)-ND), toy.max=coex.max,
          )
  # If needed, undo the log transformation
  if (!log.version) result <- undo.log.version(result)
  return(result)
}

# Visualize FCT data frame
plot.fct.df <- function(df, log.version=F,
    xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL, 
    nice.axes=T, zero.linetype=NULL, optimal.linetype=NULL,
    break.inc=0.1, label.inc=0.5, label.format="%1.1f") {
  result <- ggplot(df) +
    # Ribbons for coexistence outcomes
    geom_ribbon(aes(x=ND, ymin=coex.min, ymax=if_else(!is.na(toy.min), toy.min, coex.max)), fill=coex.col) +
    geom_ribbon(aes(x=ND, ymin=ass.min, ymax=ass.max), fill=ass.col) +
    # Ribbon for (transgressive) overyielding
    geom_ribbon(aes(x=ND, ymin=toy.min, ymax=coex.max), fill=greens[2]) +
    # Lines for coexistence boundaries
    geom_line(aes(x=ND, y=pers.1)) +
    geom_line(aes(x=ND, y=pers.2)) +
    # Line for overyielding boundary
    geom_line(aes(x=ND, y=toy.min)) +
    # Add labels
    ggplot2::xlab(xlab) + ggplot2::ylab(ylab)
  # If requested, add lines at x, y = 0
  if (!is.null(zero.linetype)) result <- result +
    geom_hline(aes(yintercept=1), linetype=zero.linetype) + 
    geom_vline(aes(xintercept=0), linetype=zero.linetype)
  if (!is.null(optimal.linetype)) result <- result +
    geom_line(aes(x=ND, y=optimal), linetype=optimal.linetype)
  if (nice.axes) result <- result +
    scale_x_continuous(breaks=breaks_width(break.inc), expand=c(0,0),
                     labels=get.label(xlim, label.inc, label.format)) +
    scale_y_continuous(breaks=breaks_width(break.inc), expand=c(0,0),
                     labels=get.label(ylim, label.inc, label.format))
  if (log.version) xtrans <- "log"
  else xtrans <- "to_nd"
  result <- result +
    coord_trans(x=xtrans, y="log", xlim=xlim, ylim=ylim, expand=F)
  return(result)
}

# Return a version of plot `p` with the layers in `classes` brought to the top.
# If `rev` is TRUE, the layers will be brought to the bottom.
bring.geoms.to.top <- function(p, classes, rev=F) {
  is.class <- function(layer, classes) any(class(layer$geom) %in% classes)
  x <- sapply(p$layers, function(x) is.class(x, classes))
  # print(x)
  if (rev) x <- !x
  p$layers <- p$layers[order(x)]
  return(p)
}

# Return a version of plot `p` with the layers in `classes` removed.
# If `keep` is TRUE, the layers will be kept and the others removed.
remove.geoms <- function(p, classes, keep=F) {
  is.class <- function(layer, classes) any(class(layer$geom) %in% classes)
  x <- sapply(p$layers, function(x) !is.class(x, classes))
  # print(x)
  if (keep) x <- !x
  p$layers <- p$layers[x]
  return(p)
}