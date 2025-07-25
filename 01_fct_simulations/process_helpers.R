#!/usr/bin/env Rscript
# 01_fct_simulations/process_helpers.R
# ====================================
# Author: Joe Wan
# Helper functions for calculating the effect of changing MCT/FCT components.

# Generate simulated data for a process
get.process.df <- function(nd, fr, Y1, Y2) {
  result <- expand_grid(nd=nd, fr=fr, Y1=Y1, Y2=Y2) %>%
  mutate(rho=1-nd) %>%
  as_tibble() %>% rowwise() %>% 
  mutate(result=list(with(
    get.equilibrium(get.alphas(rho, fr, Y1, Y2)), 
    data.frame(N1=N1, N2=N2, biomass=N1+N2)))) %>%
  unnest(result) %>% ungroup() %>% as.data.frame()
}

# Create appropriate color scale
get.scale <- function(colors=NULL, legend.name=NULL, 
    parse.legend.labels=T, is.fill=F) {
  # Set up color scale (a bit convoluted because we want to allow for 
  # color, for no legend name, and for parsing)
  if (!is.fill) {
    if (!is.null(colors)) myscale <- 
      function(...) scale_color_manual(values=colors, ...)
    else myscale <- 
      function(...) scale_color_discrete(...)
  } else {
    if (!is.null(colors)) myscale <- 
      function(...) scale_fill_manual(values=colors, ...)
    else myscale <- 
      function(...) scale_fill_discrete(...)
  }
  if (parse.legend.labels) result <- myscale(legend.name, labels=parse_format())
  else if (!is.null(legend.name)) result <- myscale(legend.name)
  else result <- myscale()
  return(result)
}

# Plot a process data frame
plot.process.df <- function(df, x.expr, facet.expr=NULL,
    show.Y1=T, show.Y2=F, yield.lt=axes.lt,
    highlight.condition=NULL, highlight.size=3,
    show.oy=T, oy.color=greens[2], oy.alpha=0.1,
    add.scale=T,
    colors=NULL, legend.name=NULL, parse.legend.labels=T,
    xlab=NULL, ylab=biomass.label) {
  # Add x and facet columns to the data frame (so we can refer to them using these names)
  if (is.null(facet.expr)) facet.expr <- T
  df <- mutate(df, x=eval(x.expr), facet=eval(facet.expr))
  # Add group info, if needed
  # if (!("group" %in% names(df))) df <- mutate(df, group=T)
  if (!("group.name" %in% names(df)) && "group" %in% names(df)) df <- mutate(df, group.name=group)
  # Build the plot
  result <- ggplot(df)
  # Plot intrinsic yields, if needed
  if (show.Y1) result <- result + geom_line(aes(x=x, y=Y1), linetype=yield.lt)
  if (show.Y2) result <- result + geom_line(aes(x=x, y=Y2), linetype=yield.lt)
  # Plot biomass as separate lines for each group
  if (is.null(highlight.condition)) {
    if ("group.name" %in% names(df)) result <- result +
      geom_line(aes(x=x, y=biomass, group=group, 
          color=as.factor(group.name)), linewidth=line.size)
    # else if ("group" %in% names(df)) result <- result +
      # geom_line(aes(x=x, y=biomass, group=group, color=as.factor(group)), linewidth=line.size)
    else result <- result + geom_line(aes(x=x, y=biomass), linewidth=line.size, color=colors[1])
  } else {
    result <- result +
      # The line does not need a color
      geom_line(aes(x=x, y=biomass), linewidth=line.size) +
      # Highlight points
      geom_point(aes(x=x, y=biomass, fill=as.factor(group.name)), 
                 size=highlight.size, 
                 data=filter(df, eval(highlight.condition)),
                 shape=21, color='white')
  }
  # If requested, add ribbon highlighting overyield
  if (show.oy) {
    # Facet is the same for all rows if we haven't ordered facet, so safe
    # to include it as a grouping variable.
    oy.df <- df %>% group_by(x, Y1, facet) %>%
      summarise(max.biomass=max(biomass)) %>% ungroup()
    # Plot using the newly calculated max.biomass
    result <- result + 
      geom_ribbon(aes(x=x, ymin=Y1, group=facet,
          ymax=ifelse(max.biomass>Y1, max.biomass, NA)), 
        fill=oy.color, alpha=oy.alpha, data=oy.df)
  }
  # Add scale if needed
  if (("group" %in% names(df) || !is.null(highlight.condition)) && add.scale) result <- result +
    get.scale(colors, legend.name, parse.legend.labels, is.fill=!is.null(highlight.condition))
  else result <- result + theme(legend.position="none")
  # Add labels
  result <- result + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)
  return(result)
}

# Plot arrows on top of already plotted base plot
plot.inset <- function(df, base.plot, x.expr, start.min=T,
    size=arrow.size, arrow=custom.arrow, avoid=arrow.avoid,
    colors=NULL) {
  if (start.min) myf <- min
  else myf <- max
  df <- mutate(df, x=eval(x.expr))
  result <- base.plot +
    geom_point(aes(x=nd, y=fr, group=group, color=as.factor(group)), 
              data=filter(df, x==myf(df$x))) +
    geom_line(aes(x=nd, y=fr, group=group, color=as.factor(group)), 
              linewidth=size, arrow=arrow, 
              data=filter(df, x>=myf(df$x)+avoid)) +
    get.scale(colors, NULL) +
    theme_void() + theme(legend.position="none", aspect.ratio=ratio)
  return(result)
}