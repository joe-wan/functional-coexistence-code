#!/usr/bin/env Rscript
# 01_fct_simulations/05_multifunctionality.R
# ==========================================
# Author: Joe Wan
# Generates figures illustrating multifunctional overyielding (Figure 9a-b).

# Data manipulation
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)

##### 0a. MCT helper functions #####
source("mct_helpers.R")

##### 0b. Define parameters of system and ranges for visualizations #####
source("mct_params.R")

# New parameters for function
phi1 <- (1.02/1.2)^2
phi2 <- 1/Y2
Phi1 <- phi1 * Y1
Phi2 <- phi2 * Y2

##### 0c. Define visual language #####
source("graph_params.R")
source("graph_helpers.R")

##### 0d. Specific parameters and helpers for varying MCT processes #####
source("process_params.R")
source("process_helpers.R")

# Override default FR range for visualization clarity
eq.fr <- seq(0.65, 1.55, by=resolution/10)

##### 1. Plot FCT space, now with two functions #####
# Get data for multifunctionality analysis
multi.data <- get.fct.df(seq(min.nd, max.nd.wide, by=resolution), Y1=Y1, Y2=Y2, log.version=T) %>%
  mutate(optimal.func = log(Phi1/Phi2)/2,
         toy.func.min = coex.min,
         toy.func.max = case_when(ND >= -optimal.func ~ log(Phi1/Phi2)+ND),
         optimal.multi = optimal + optimal.func) %>%
  undo.log.version() %>%
  mutate(ND = trad.to.log(ND))

# Create multifunctionality plot
multi.plot <- plot.fct.df(multi.data,
      log.version=T,
      xlim=trad.to.log(fct.x.wide), ylim=fct.y,
      zero.linetype=axes.lt, optimal.linetype=NULL,
      xlab=nd.label, ylab=fr.label) +
  # Ribbon for function overyielding only
  geom_ribbon(aes(x=ND, ymin=coex.min,
      ymax=pmin(toy.func.max, if_else(!is.na(toy.min), toy.min, Inf))), fill=blues[1]) +
  # Ribbon for multifunctional overyielding
  geom_ribbon(aes(x=ND, ymax=toy.func.max,
      ymin=if_else(toy.min<toy.func.max, toy.min, NA)), fill=blues[3]) +
  # Boundary line for function overyielding
  geom_line(aes(x=ND, y=toy.func.max)) +
  theme_cowplot() + theme(aspect.ratio=ratio.wide)

# Ensure lines appear on top
multi.plot <- bring.geoms.to.top(multi.plot, c('GeomLine', 'GeomVline', 'GeomHline'))

# Add lines for optima (biomass, function, multifunctional)
multi.plot <- multi.plot +
  geom_line(aes(x=ND, y=if_else(!is.na(toy.max-toy.min), optimal, NA)), linetype=highlight.lt) +
  geom_line(aes(x=ND, y=if_else(!is.na(toy.func.max-toy.func.min), optimal.func, NA)), linetype=highlight.lt) +
  geom_line(aes(x=ND, y=optimal.multi), linetype=highlight.lt)
if (show.plots) print(multi.plot)

##### 3b. Equalization (increasing FR) #####
# Generate simulated data with function-based metrics
eq.df <- get.process.df(eq.nd, eq.fr, Y1, Y2) %>%
  mutate(group=nd, group.name=paste(nd.prefix, nd),
         phi1=phi1, phi2=phi2, Phi1=Phi1, Phi2=Phi2,
         func=phi1*N1+phi2*N2)

# Helpers for ribbon highlights
ribbon.alpha <- 0.1
ribbon.min <- 0.675
ribbon.max <- 1.1
bg.highlight <- function(condition, color, alpha=ribbon.alpha) geom_ribbon(aes(
    x=fr, ymin=if_else(eval(condition), ribbon.min, NA), ymax=ribbon.max, group=group), fill=color, alpha=alpha)
above.highlight <- function(min.expr, max.expr, color, alpha=0.5) geom_ribbon(aes(
    x=fr, ymin=if_else(eval(max.expr)>eval(min.expr), eval(min.expr), NA), ymax=eval(max.expr), group=group), fill=color, alpha=alpha)

# Plot biomass and function under FR variation
eq.biomass <- plot.process.df(eq.df, expr(fr), facet.expr=expr(nd),
    show.Y1=T, show.Y2=F, show.oy=F, colors=eq.colors,
    xlab=fr.label,
    ylab=expression(atop(paste('tot. biomass, ', Sigma * N), paste('tot. function, ', Sigma * Phi)))) +
  geom_line(aes(x=fr, y=func, group=group, color=as.factor(group.name)), size=line.size, linetype='dashed') +
  bg.highlight(expr(biomass>1&func<1), greens[2]) +
  bg.highlight(expr(biomass<1&func>1), blues[2]) +
  bg.highlight(expr(biomass>1&func>1), blues[4]) +
  above.highlight(expr(pmax(func, 1)), expr(biomass), greens[2]) +
  above.highlight(expr(pmax(biomass, 1)), expr(func), blues[2]) +
  above.highlight(expr(1), expr(pmin(biomass, func)), blues[4]) +
  scale_x_continuous(breaks=breaks_width(0.1), expand=c(0,0), labels=custom_label('%1.1f', c(0.5, 1, 1.5, 2))) +
  scale_y_continuous(expand=c(0,0), breaks=breaks_width(0.2)) +
  facet_wrap(~facet, ncol=1) +
  theme_cowplot() + theme(strip.background=element_blank(), strip.text=element_blank(), legend.position="none")
if (show.plots) print(eq.biomass)

##### 3. Export plots #####
ggsave(paste0(outpath, '09_a_fct_space_multi_wide.pdf'), multi.plot + th,
       device=cairo_pdf, width=out.width*1.5, height=out.height*1.5)
ggsave(paste0(outpath, '09_b_equalization_multi_biomass.pdf'), eq.biomass + th,
       device=cairo_pdf, width=out.width*0.75, height=out.height*1.5)