#!/usr/bin/env Rscript
# 01_fct_simulations/05_multifunctionality.R
# ==========================================
# Author: Joe Wan
# Generates figures illustrating multifunctional overyielding (Figure 6a-b).

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

# New params for function
phi1 <- (1.02/1.2)^2
phi2 <- 1/Y2
Phi1 <- phi1*Y1
Phi2 <- phi2*Y2


##### 0c. Define visual language #####
source("graph_params.R")
source("graph_helpers.R")

# Alpha of ribbons for overyielding
ribbon.alpha <- 0.1

##### 0d. Specific parameters and helpers for varying MCT processes #####
source("process_params.R")
source("process_helpers.R")

# Override the default FR range to make things easier to see
eq.fr <- seq(0.65, 1.55, by=resolution/10)


##### 1. Plot FCT space, now with two functions #####
# Get data for the plot
multi.data <- get.fct.df(seq(min.nd, max.nd.wide, by=resolution), 
    Y1=Y1, Y2=Y2, log.version=T) %>%
  mutate(optimal.func=log(Phi1/Phi2)/2,
         toy.func.min=coex.min,
         toy.func.max=case_when(ND>=-optimal.func~log(Phi1/Phi2)+ND),
         optimal.multi=optimal + optimal.func) %>%
  undo.log.version()

# Plot the FCT space
multi.plot <- plot.fct.df(multi.data, 
      xlim=fct.x.wide, ylim=fct.y,
      zero.linetype=axes.lt, optimal.linetype=NULL,
      xlab=nd.label, ylab=fr.label) +
  # Add new ribbon for function overyield only
  geom_ribbon(aes(x=ND, ymin=coex.min,
      ymax=pmin(toy.func.max, if_else(!is.na(toy.min), toy.min, Inf))), 
    fill=blues[1]) +
  # ... and for multifunctional overyielding
  geom_ribbon(aes(x=ND, ymax=toy.func.max,
      ymin=if_else(toy.min<toy.func.max, toy.min, NA)),
    fill=blues[3]) +
  # Boundary for function overyield  
  geom_line(aes(x=ND, y=toy.func.max)) +
  theme_cowplot() + theme(aspect.ratio=ratio.wide)

# Tweak the plot:
# ---------------
# Lines need to be on top of the regions they delineate
multi.plot <- bring.geoms.to.top(multi.plot, c('GeomLine', 'GeomVline', 'GeomHline'))
# Save version before we add optimum lines
multi.plot.fade <- multi.plot
# Add lines for optima
multi.plot <- multi.plot +
  # Biomass
  geom_line(aes(x=ND, y=if_else(!is.na(toy.max-toy.min), optimal, NA)), 
    linetype=highlight.lt) +
  # Function
  geom_line(aes(x=ND, y=if_else(!is.na(toy.func.max-toy.func.min), optimal.func, NA)), 
    linetype=highlight.lt) +
  # Multifunctional overyielding
  geom_line(aes(x=ND, y=optimal.multi), 
    linetype=highlight.lt)
if (show.plots) print(multi.plot)

# Make faded version for insets
multi.plot.fade <- remove.geoms(multi.plot.fade, c('GeomVline', 'GeomHline')) +
  fade.box(xlim=fct.x.wide, ylim=fct.y) +
  coord_trans(x="to_nd", y="log",
              ylim=c(min.fr, max.fr), xlim=1-exp(-c(min.nd, max.nd)),
              expand=F)
if (show.plots) print(multi.plot.fade)


##### 3b. Equalization (increasing FR) #####
# From process_params.R:
# eq.colors, eq.nd
# Overwritten in this script:
# eq.fr

# Generate simulated data
eq.df <- get.process.df(eq.nd, eq.fr, Y1, Y2) %>%
  mutate(group=nd, group.name=paste(nd.prefix, nd)) %>%
  mutate(phi1=phi1, phi2=phi2, Phi1=Phi1, Phi2=Phi2) %>%
  mutate(func=phi1*N1+phi2*N2)

# Helper functions to reduce retyping
ribbon.min <- 0.675
ribbon.max <- 1.1
bg.highlight <- function(condition, color, alpha=ribbon.alpha) geom_ribbon(aes(
    x=fr, ymin=if_else(eval(condition), ribbon.min, NA), ymax=ribbon.max, group=group),
    fill=color, alpha=alpha)
above.highlight <- function(min.expr, max.expr, color, alpha=0.5) geom_ribbon(aes(
    x=fr, ymin=if_else(eval(max.expr)>eval(min.expr), eval(min.expr), NA), 
    ymax=eval(max.expr), group=group), fill=color, alpha=alpha)

# Plot biomass as we vary FR
eq.biomass <- plot.process.df(eq.df, expr(fr), facet.expr=expr(nd),
    show.Y1=T, show.Y2=F, show.oy=F, colors=eq.colors,
    xlab=fr.label, 
    ylab=expression(atop(paste('tot. biomass, ', Sigma * N), paste('tot. function, ', Sigma * Phi)))) +
  # # Add vertical line at optimal FR
  # geom_vline(xintercept=sqrt(Y1*Phi1/Y2/Phi2), linetype=highlight.lt) +
  geom_line(aes(x=fr, y=func, group=group, color=as.factor(group.name)), size=line.size, linetype='dashed') +
  # Highlight the background showing when different forms of overyielding occur
  bg.highlight(expr(biomass>1&func<1), greens[2]) +
  bg.highlight(expr(biomass<1&func>1), blues[2]) +
  bg.highlight(expr(biomass>1&func>1), blues[4]) +
  # Highlight the actual difference between the best monoculture and the community
  above.highlight(expr(pmax(func, 1)), expr(biomass), greens[2]) +
  above.highlight(expr(pmax(biomass, 1)), expr(func), blues[2]) +
  above.highlight(expr(1), expr(pmin(biomass, func)), blues[4]) + 
  # Make the x axis same as FCT plot's y axis
  scale_x_continuous(breaks=breaks_width(0.1), expand=c(0,0),
                     labels=custom_label('%1.1f', c(0.5, 1, 1.5, 2))) +
  scale_y_continuous(expand=c(0,0), breaks=breaks_width(0.2)) +
  facet_wrap(~facet, ncol=1) +
  # Make it the same aspect ratio as other panels
  cowplot::theme_cowplot() + 
  theme(strip.background = element_blank(),
    strip.text = element_blank(),
    legend.position="none")
if (show.plots) print(eq.biomass)

# Make inset plot
eq.inset <- eq.df %>%
  plot.inset(multi.plot.fade, 
    expr(fr), start.min=T,
    size=arrow.size, arrow=custom.arrow, avoid=arrow.avoid,
    colors=eq.colors) +
  theme_void() + theme(legend.position="none", aspect.ratio=ratio)
if (show.plots) print(eq.inset)


##### 3. Export plots #####
ggsave(paste0(outpath, '06_a_fct_space_multi_wide.pdf'), multi.plot + th, 
       device=cairo_pdf, width=out.width*1.5, height=out.height*1.5)
ggsave(paste0(outpath, '06_b_equalization_multi_biomass.pdf'), eq.biomass + th, 
       device=cairo_pdf, width=out.width*0.75, height=out.height*1.5)
ggsave(paste0(outpath, '06_b_equalization_multi_inset.pdf'), eq.inset + th, 
       device=cairo_pdf, width=out.width/2, height=out.height/2)
