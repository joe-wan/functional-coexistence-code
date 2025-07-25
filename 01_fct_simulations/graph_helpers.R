# 01_fct_simulations/graph_helpers.R
# ==================================
# Author: Joe Wan
# Helper functions for plotting MCT and FCT figures.

library(dplyr)
library(cowplot)
library(ggplot2)
library(stringr)

##### Helper functions for plotting #####
# Transform MCT's "niche difference" to the "linear" ND (= -log(rho))
to_nd_trans <- function() trans_new("to_nd", function(x) -log(1-x), function(x) 1-exp(-x))

# Custom label: format numbers according to sprintf format if they are in the set 'x'
custom_label <- function(format, x, digits=5) {
  return(function(breaks) {
    breaks.rounded <- round(breaks, digits)
    x.rounded <- round(x, digits)
    labels <- sprintf(format, breaks)
    labels <- if_else(breaks.rounded %in% x.rounded, labels, '')
    return(labels)
  })
}

# Create label function to include r[1] and r[2], with resolution inc and print format format
get.label <- function(r, inc, format) {
  if (is.null(r) || is.null(inc) || is.null(format)) return(NULL)
  rounded.r <- inc * c(floor(r[1]/inc), ceiling(r[2]/inc))
  breaks <- seq(rounded.r[1], rounded.r[2], by=inc)
  return(custom_label(format, breaks))
}

# Return rectangle to fade out the plot
fade.box <- function(xlim, ylim, alpha=0.85) geom_rect(
  xmin=xlim[1], xmax=xlim[2], ymin=ylim[1], ymax=ylim[2],
  fill='white', alpha=alpha, data=data.frame(x=1))

##### Raw transformation helpers #####
log.to.trad <- function(x) 1 - exp(-x)
trad.to.log <- function(x) -log(1 - x)

##### MCT functions #####
get.mct.df <- function(nds) {
  result <- data.frame(ND=nds) %>%
    mutate(pers.1=-ND, pers.2=ND,
           coex.min=case_when(ND>=0 ~ pers.1), coex.max=case_when(ND>=0 ~ pers.2),
           ass.min=case_when(ND<=0 ~ pers.2), ass.max=case_when(ND<=0 ~ pers.1))
}

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
undo.log.version <- function(df) df %>% 
  mutate(across(-ND, ~ exp(.x))) %>% 
  mutate(ND=1-exp(-ND))

get.log.version <- function(df) df %>% 
  mutate(across(-ND, ~ log(.x))) %>% 
  mutate(ND=-log(1-ND))

get.fct.df <- function(nds, Y1=NULL, Y2=NULL, log.version=F) {
  if (any(is.null(c(Y1, Y2)))) stop('Y1 and Y2 must be specified')
  if (is.null(log.version)) stop('log.version must be specified')
  result <- get.mct.df(nds) %>%
    mutate(optimal=log(Y1/Y2)/2,
           toy.min=case_when(ND>=optimal ~ log(Y1/Y2)-ND), toy.max=coex.max)
  if (!log.version) result <- undo.log.version(result)
  return(result)
}

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
  if (log.version) xtrans <- "identity" else xtrans <- "to_nd"
  result <- result +
    coord_trans(x=xtrans, y="log", xlim=xlim, ylim=ylim, expand=F)
  return(result)
}

##### Tools for compositing figures #####
get.ggtitle.string <- function(letter, title, spaces=2,
    parens=T, bold=T) {
    spaces <- paste(rep(" ", spaces), collapse="<i/>")
    if (!is.null(letter) && str_length(letter)>0) {
        if (letter=="c") letter <- paste0(letter, "<i/>")
        if (parens) letter <- paste0("(", letter, ")")
        if (bold) letter <- paste0("**", letter, "**")
    }
    return(paste0(letter, spaces, title))
}

format.title <- function(...) {
    return(annotate_figure(NULL, b=richtext_grob(
    get.ggtitle.string(...), 
    hjust=0, x=0, vjust=0, y=0, margin=margin(l=11, b=-10.5, unit="pt"),
    gp=gpar(fontfamily="Helvetica", fontsize=16))))
}

bring.geoms.to.top <- function(p, classes, rev=F) {
  is.class <- function(layer, classes) any(class(layer$geom) %in% classes)
  x <- sapply(p$layers, function(x) is.class(x, classes))
  if (rev) x <- !x
  p$layers <- p$layers[order(x)]
  return(p)
}

remove.geoms <- function(p, classes, keep=F) {
  is.class <- function(layer, classes) any(class(layer$geom) %in% classes)
  x <- sapply(p$layers, function(x) !is.class(x, classes))
  if (keep) x <- !x
  p$layers <- p$layers[x]
  return(p)
}
