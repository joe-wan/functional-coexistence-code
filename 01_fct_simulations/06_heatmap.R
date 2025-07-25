#!/usr/bin/env Rscript
# 01_fct_simulations/06_heatmap.R
# ==========================================
# Author: Joe Wan
# Generates figures illustrating heatmaps of productivity and additive partition
# components (Figure 3a-b).

library(dplyr)
library(ggplot2)
library(tidyr)
library(cowplot)
library(scales)


##### 0a. MCT helper functions #####
source("mct_helpers.R")

##### 0b. Define parameters of system and ranges for visualizations #####
source("mct_params.R")

##### 0c. Define visual language #####
source("graph_helpers.R")
source("graph_params.R")

source("heatmap_helpers.R")
source("heatmap_params.R")

##### 1. Calculations for heatmaps, color scales, etc. #####

# Calculate new aspect ratio
max.nd.heatmap <- log(max.fr)
heatmap.aspect <- (log(max.fr) - log(min.fr)) / max.nd.heatmap

# Calculate the heatmap resolution
delta <- (log(max.fr) - log(min.fr)) / heatmap.res
ND=seq(0, max.nd.heatmap, delta)
FD=seq(log(min.fr), log(max.fr), length.out=(heatmap.res+1))

# Set up the heatmap data frame
heatmap.df <- expand.grid(ND=ND, FD=FD) %>%
  mutate(rho=exp(-ND), FR=exp(FD), Y1=Y1, Y2=Y2) %>%
  rowwise() %>% 
  mutate(result=list(with(
    get.equilibrium(get.alphas(rho, FR, Y1, Y2)), 
    data.frame(N1=N1, N2=N2, biomass=N1+N2)))) %>%
  unnest(result) %>% ungroup() %>%
  # Fix the ND = 0 column
  mutate(biomass=case_when(
    ND != 0 ~ biomass,
    FD > 0 ~ Y1,
    FD < 0 ~ Y2)) %>%
  # Calculate partition components
    mutate(CE=(N1/Y1 + N2/Y2 - 1)*(Y1+Y2)/2) %>%
    mutate(CE=case_when(ND==0 ~ 0, TRUE ~ CE)) %>%
    mutate(SE=biomass-(Y1+Y2)/2-CE)

contour.delta <- (Y1 - Y2 - fudge) / contour.res
max.biomass <- max(heatmap.df$biomass)
contour.max <- ceiling((max.biomass - Y1) / contour.delta) * contour.delta + Y1
contour.min <- Y2+fudge-contour.delta
contour.breaks <- seq(contour.min, contour.max, by=contour.delta)

# Get colors for breaks
# Get position at which green[2] should appear
target.L <- col.to.lab(greens[2])[1]
target.pos <- contour.res + 2
delta.L <- (1 - target.L/100)/(target.pos-1)
breaks <- seq(0, by=delta.L, length.out=length(contour.breaks))
gray.cols <- gray.fun(breaks)
green.cols <- green.fun(breaks)
cols <- c(gray.cols[1:3], green.cols[4:length(green.cols)])

roy.col <- green.cols[4] # "#C4EB8C"

# Fake a stepped scale by repeating each color x times
fake.stepped <- function(cols, x=25) rep(cols, each=x)

# Set up a diverging scale for SE
max.SE <- max(heatmap.df$SE, na.rm=T) # 0.1527778
max.CE <- max(heatmap.df$CE, na.rm=T) # 0.3418616
component.step <- max.SE/3

# Calculate colors of the breaks
component.breaks <- seq(-q,q)*component.step
pos.cols <- blue.fun(seq(0,k,length.out=q+1))
neg.cols <- orange.fun(seq(0,k,length.out=q+1))
component.cols <- c(rev(neg.cols[-1]), pos.cols[-1])

# Calculate which brackets are used in each plot
start.SE <- first(which(component.breaks > -max.SE-1e-3)) - 1

##### 2. Create heatmaps #####
# Helpers to perform repetitive graphing tasks
add.roy <- function(x) x +
    # Plot an unfilled contour plot with the ROY boundary
    geom_contour(aes(x=xtran(ND), y=ytran(FD), z=biomass), 
        breaks=c((Y1+Y2)/2),
        linetype=roy.lt, color="black",
        data=heatmap.df)

add.toy <- function(x) x +
    # Plot an unfilled contour plot with the TOY boundary
    geom_line(aes(x=xtran(ND), y=ytran(TOY)),
        linetype=toy.lt)

# Shared base heatmap
base.hm <- data.frame(ND=seq(0, max.nd.heatmap, length.out=1000)) %>%
    mutate(TOY=case_when(ND>=log(Y1/Y2)/2~-ND+log(Y1/Y2))) %>%
    ggplot()

green.hm <- base.hm + 
    # Plot points with color to get our legend right
    geom_point(aes(x=xtran(ND), y=ytran(FD), color=if_else(biomass==Y2, biomass-contour.delta/2, biomass)),
        size=-1,
        data=heatmap.df) +
    scale_color_gradientn("biomass", colors=fake.stepped(cols), 
        limits=c(contour.min, contour.max+contour.delta),
        breaks=c(0.7,1,1.3),
        # Use plotmath in labels
        labels=c(expression(paste(0.7, ' ', 
        (K[2]))), expression(paste('1.0', ' ',  
        (K[1]))), 1.3),
        guide=guide_colorbar(frame.colour = "black"))

biomass.hm <- green.hm +
    geom_contour_filled(aes(x=xtran(ND), y=ytran(FD), z=if_else(FD>ND, biomass-1e-3, biomass)), 
            breaks=contour.breaks,
            data=heatmap.df) +
        scale_fill_manual(values=cols, guide="none")
biomass.hm <- biomass.hm %>% add.roy() %>% add.toy() %>% finish.heatmap()
if (show.plots) print(biomass.hm)

orbl.hm <- base.hm + 
    # Plot points with color to get our legend right
    geom_point(aes(x=xtran(ND), y=ytran(FD), color=SE),
        size=-1,
        data=heatmap.df) +
    scale_color_gradientn("value\n(SE or CE)", colors=fake.stepped(component.cols),  
        limits=c(-q+1,q-1)*component.step,
        guide=guide_colorbar(frame.colour = "black"))
# orbl.hm %>% finish.heatmap()

# Heat map of SE
selection.hm <- orbl.hm + 
    # Plot a niSE contour plot instead    
    geom_contour_filled(aes(x=xtran(ND), y=ytran(FD), z=SE+sign(SE)*1e-3), 
        breaks=component.breaks,
        data=heatmap.df) +
    # scale_fill_manual(values=c('white', component.cols[(q+1):(2*q)])) +
    scale_fill_manual(values=component.cols[start.SE:(2*q)], guide="none")
selection.hm <- selection.hm %>% finish.heatmap()
if (show.plots) print(selection.hm)

# Heat map of CE
complementarity.hm <- orbl.hm + 
    # Plot a niCE contour plot instead    
    geom_contour_filled(aes(x=xtran(ND), y=ytran(FD), z=CE-1e-3), 
        breaks=component.breaks,
        data=heatmap.df) +
    scale_fill_manual(values=c('white', component.cols[(q+1):(2*q)]), guide="none")
complementarity.hm <- complementarity.hm %>% finish.heatmap()
if (show.plots) print(complementarity.hm)


# Heat map of SE
ytran2 <- xtran
ggplot(filter(heatmap.df, ND>0)) + 
    # geom_contour_filled(aes(x=ND, y=FD, z=CE)) +
    # scale_fill_brewer(palette="RdBu") +
  # Diverging color scale
  geom_tile(aes(x=xtran(ND), y=ytran2(FD), fill=CE)) + 
  cmocean::scale_fill_cmocean(name="balance",
    direction=-1,
    limits=c(-0.35,0.35)) +
    # Make y axis in unlogged units
    coord_trans(x="identity", y="identity", 
        xlim=c(0, xtran(max.nd.heatmap)), ylim=c(ytran2(log(min.fr)), ytran2(log(max.fr))),
        expand=F) +
    theme_cowplot() + 
    xlab(nd.label) + ylab(fr.label) +
    # Hack to get ticks on one side outside
    theme(legend.ticks.length = unit(c(-3, 0), 'pt'), 
        legend.ticks = element_line(color='black'),
        aspect.ratio=heatmap.aspect)

# Save and composite the plots
ggsave(paste0(outpath, '03_a_heatmap_biomass.pdf'), biomass.hm + th, 
    device=cairo_pdf, width=out.width*1.5, height=out.height*1.5)
ggsave(paste0(outpath, '03_b_heatmap_selection.pdf'), selection.hm + th, 
    device=cairo_pdf, width=out.width*1.5, height=out.height*1.5)
ggsave(paste0(outpath, '03_b_heatmap_complementarity.pdf'), complementarity.hm + th, 
    device=cairo_pdf, width=out.width*1.5, height=out.height*1.5)
