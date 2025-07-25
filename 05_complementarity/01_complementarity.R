#!/usr/bin/env Rscript
# 05_complementarity/01_complementarity.R
# =======================================
# Author: Joe Wan
# Generates plot comparing niche and fitness metrics (Figura S3).

# Data manipulation
library(tidyr)
library(dplyr)
library(stringr)

# Plotting
library(ggplot2)
library(RColorBrewer)
library(cowplot)


##### 0a. Set up graphing options #####
show.plots <- interactive() # Will suppress plot output if not REPL
outpath <- 'outputs/'       # Directory to save plots
out.width <- 3.25           # Width of main plots
out.height <- 3.25          # Height of main plots
# Font, etc. will be added by this theme:
th <- theme(text=element_text(family="Helvetica"))

# Ensure the path ends with a slash
if (!str_ends(outpath, '/')) outpath <- paste0(outpath, '/')
# Ensure output directory exists
dir.create(outpath, showWarnings = FALSE)

# Color scheme
n.col <- 5
rb <- brewer.pal(n.col*2+1, 'RdBu')
# rb <- rev(cmocean('balance')(n.col*2+1))
re <- rb[1:n.col]
bu <- rev(rb)[1:n.col]

shapes.scale <- c(19, 1)
col.scale <- c(re[2], re[2], bu[3])
lin.scale <- c('solid', 'dashed', 'solid')

axes <- 'solid'
highlight <- 'dashed'
pointsize <- 2
pathsize <- 0.75
arrw <- arrow(length=unit(4,'pt'), type='closed')
ass.alpha <- 0.1
coe.alpha <- 0.2
line.colors <- c(bu[c(2,3)], 'gray50', re[c(4,3)])


##### 1. Generate data #####
############################
# Resolution: focal and non-focal parameter
maj <- 0.01
min <- 0.01

# Max/min value for geometric simulation
s.geo <- 5

# Generate processes for geometric definition
grid.geo <- bind_rows(
  expand_grid(var='FD', 
              ND.geo=seq(-s.geo, s.geo, by=min), 
              FD.geo=seq(-s.geo, s.geo, by=maj)),
  expand_grid(var='ND',
              ND.geo=seq(-s.geo, s.geo, by=maj), 
              FD.geo=seq(-s.geo, s.geo, by=min))) %>%
  mutate(grid='geo', 
         value=if_else(var=='ND', ND.geo, FD.geo),
         other=if_else(var!='ND', ND.geo, FD.geo),
         T1=1-exp(-ND.geo-FD.geo), T2=1-exp(-ND.geo+FD.geo))

# Max/min value for arithmetic simulation
s.ari <- 5

# Generate processes for arithmetic definition
grid.ari <- bind_rows(
  expand_grid(var='FD', 
              ND.ari=seq(-s.ari, s.ari, by=min), 
              FD.ari=seq(-s.ari, s.ari, by=maj)),
  expand_grid(var='ND',
              ND.ari=seq(-s.ari, s.ari, by=maj), 
              FD.ari=seq(-s.ari, s.ari, by=min))) %>%
  mutate(grid='ari', 
         value=if_else(var=='ND', ND.ari, FD.ari),
         other=if_else(var!='ND', ND.ari, FD.ari),
         T1=ND.ari+FD.ari, T2=ND.ari-FD.ari)

# Combine calculations
grid <- bind_rows(grid.geo, grid.ari) %>% #, grid.che, grid.tig) %>%
  mutate(ND.ari=(T1+T2)/2, FD.ari=T1-ND.ari,
         S1=1-T1, S2=1-T2,
         rho=sqrt(S1*S2), f1.f2=rho/S1,
         ND.geo=-log(rho), FD.geo=log(f1.f2),
         ND.che=1-rho, FD.che=f1.f2,
         scaling=2*abs(1/(1-rho^2) - 1/(T1+T2)),
         V1=T1*scaling, V2=T2*scaling,
         ND.com=(V1+V2)/2, FD.com=V1-ND.com) %>%
  mutate(grp=paste(var, paste0('(', grid, ')'), '=', sprintf('%.2f', value)))

# Generate regions for coexistence outcomes
regions.res <- 0.001
regions.geo <- data_frame(ND.geo=seq(-s.geo, s.geo, by=regions.res)) %>%
  mutate(inv1=ND.geo, inv2=-ND.geo,
         ass.min=case_when(ND.geo<0 ~ ND.geo),
         ass.max=case_when(ND.geo<0 ~ -ND.geo),
         coe.min=case_when(ND.geo>0 ~ -ND.geo),
         coe.max=case_when(ND.geo>0 ~ ND.geo)) %>%
  mutate(rho=exp(-ND.geo), ND.che=1-rho)
regions.ari <- data_frame(ND.ari=seq(-s.ari, s.ari, by=regions.res)) %>%
  mutate(inv1=ND.ari, inv2=-ND.ari,
         ass.min=case_when(ND.ari<0 ~ ND.ari),
         ass.max=case_when(ND.ari<0 ~ -ND.ari),
         coe.min=case_when(ND.ari>0 ~ -ND.ari),
         coe.max=case_when(ND.ari>0 ~ ND.ari))


##### 2. Illustrate complementarity as a measure of niche #####
# Tolerance for identifying points (accounts for floating point rounding error)
err <- 1e-3

# Focal processes for geometric definition
paths.com.geo <- data.frame(grp=c('ND (geo) = 0.40', 
                              'ND (geo) = 0.40', 
                              'FD (geo) = 0.10'),
                        index=c(1,     2,  3),
                        start=c(0.52, -0.00,  0.05),
                        end=c( -0.00, -0.52,  0.55)) %>%
  mutate(panel='a')
# Focal processes for arithmetic definition
paths.com.ari <- data.frame(grp=c('ND (ari) = 0.40', 
                              'ND (ari) = 0.40', 
                              'FD (ari) = 0.10'),
                        index=c(4,     5,  6),
                        start=c(0.52, -0.00,  0.05),
                        end=c( -0.00, -0.52,  0.55)) %>%
  mutate(panel='b')

# Combine paths for geometric and arithmetic definitions
paths.com <- bind_rows(paths.com.ari, paths.com.geo) %>%
  left_join(grid, 'grp') %>% 
  filter(other>=pmin(start, end)-err, other<=pmax(start, end)+err) %>%
  mutate(type=case_when(abs(other-start)<err ~ 'start', 
                        abs(other-end)<err ~ 'end', T ~ 'mid')) %>%
  arrange(-abs(other-end))
regions.com <- regions.ari %>% mutate(ND.com=ND.ari) %>% select(-ND.ari)

##### 2. Make panels #####
##########################
# Helper functions
# ----------------
# Function for plotting complementarity panels
plot.com <- function(x) {
  ggplot(regions.com, aes(x=ND.com)) +
  geom_hline(yintercept=0, linetype=axes) + 
  geom_vline(xintercept=0, linetype=axes) +
  geom_ribbon(aes(ymin=ass.min, ymax=ass.max), alpha=ass.alpha) +
  geom_ribbon(aes(ymin=coe.min, ymax=coe.max), alpha=coe.alpha) +
  geom_path(aes(x=ND.com, y=FD.com, group=index, 
                color=as.factor(index),
                linetype=as.factor(index)),
            arrow=arrw, size=pathsize,
            filter(x, abs(other-end)>0.02|(index==2&abs(other-end)>0.01))) +
  geom_point(aes(x=ND.com, y=FD.com), color='white',
             data=filter(x, type!='mid'&abs(other+0.22)>1e-4)) +
  geom_point(aes(x=ND.com, y=FD.com, color=as.factor(index), shape=type), size=pointsize,
             data=filter(x, type!='mid'&abs(other+0.22)>1e-4)) +
  geom_text(aes(x=ND.com, y=FD.com+0.04, label=index),
            filter(x, type=='start')) +
  coord_equal(xlim=c(-0.2,0.401), ylim=c(-0.301,0.301), expand=0) +
  scale_linetype_manual(values=lin.scale) +
  scale_color_manual(values=col.scale) +
  scale_shape_manual(values=shapes.scale, breaks=c('start', 'end')) +
  guides(shape='none', color='none', linetype='none') +
  xlab(expression(paste('complementarity, ', C * ( F[1] + F[2] ) / 2))) + 
  ylab(expression(paste('scaled fitness diff., ',
    C * ( F[1] - F[2] ) / 2))) +
  theme_cowplot()
}

# Function for plotting arithmetic panels
plot.ari <- function(x) {
  ggplot(regions.ari, aes(x=ND.ari)) +
  geom_hline(yintercept=0, linetype=axes) + 
  geom_vline(xintercept=0, linetype=axes) +
  geom_ribbon(aes(ymin=ass.min, ymax=ass.max), alpha=ass.alpha) +
  geom_ribbon(aes(ymin=coe.min, ymax=coe.max), alpha=coe.alpha) +
  geom_path(aes(x=ND.ari, y=FD.ari, group=index, 
                color=as.factor(index),
                linetype=as.factor(index)),
            arrow=arrw, size=pathsize,
            filter(x, abs(other-end)>0.03)) +
  geom_point(aes(x=ND.ari, y=FD.ari), color='white',
             data=filter(x, type!='mid'&abs(other+0.22)>1e-4)) +
  geom_point(aes(x=ND.ari, y=FD.ari, color=as.factor(index), shape=type), size=pointsize,
             data=filter(x, type!='mid'&abs(other+0.22)>1e-4)) +
  geom_text(aes(x=ND.ari-0.08, y=FD.ari, label=index),
            filter(x, type=='start')) +
  coord_equal(xlim=c(-0.5,0.75), ylim=c(-0.625,0.625), expand=0) +
  scale_linetype_manual(values=lin.scale) +
  scale_color_manual(values=col.scale) +
  scale_shape_manual(values=shapes.scale, breaks=c('start', 'end')) +
  guides(shape='none', color='none', linetype='none') +
  ylab(expression(paste('arith. fitness diff., ', 
    (F[1] - F[2])/2 ))) + 
  xlab(expression(paste('arith. niche diff., ', 
    (F[1] + F[2])/2))) +
  theme_cowplot()
}

# Function for plotting geometric panels
plot.geo <- function(x, factor=1) {
  # `factor` is used to adjust distances since things get distorted
  ggplot(regions.geo, aes(x=ND.geo)) +
  geom_hline(yintercept=0, linetype=axes) + 
  geom_vline(xintercept=0, linetype=axes) +
  geom_ribbon(aes(ymin=ass.min, ymax=ass.max), alpha=ass.alpha) +
  geom_ribbon(aes(ymin=coe.min, ymax=coe.max), alpha=coe.alpha) +
  geom_path(aes(x=ND.geo, y=FD.geo, group=index, 
                color=as.factor(index),
                linetype=as.factor(index)),
            arrow=arrw, size=pathsize,
            filter(x, abs(other-end)>0.03)) +
  geom_point(aes(x=ND.geo, y=FD.geo), color='white',
             data=filter(x, type!='mid'&abs(other+0.22)>1e-4)) +
  geom_point(aes(x=ND.geo, y=FD.geo, color=as.factor(index), shape=type), size=pointsize,
             data=filter(x, type!='mid'&abs(other+0.22)>1e-4)) +
  geom_text(aes(x=ND.geo-factor*0.08, y=FD.geo, label=index),
            filter(x, type=='start')) +
  coord_equal(xlim=factor*c(-0.5,0.75), ylim=factor*c(-0.625,0.625), expand=0) +
  scale_linetype_manual(values=lin.scale) +
  scale_color_manual(values=col.scale) +
  scale_shape_manual(values=shapes.scale, breaks=c('start', 'end')) +
  guides(shape='none', color='none', linetype='none') +
  ylab(expression(paste('geom. fitness diff., ', 
    log * ' ' * f[1]/f[2]))) + 
  xlab(expression(paste('geom. niche diff., ', 
    - log * ' ' * rho))) +
  theme_cowplot()
}

# Generate plots
# --------------
com.a <- plot.com(filter(paths.com, panel=='a'))
com.b <- plot.com(filter(paths.com, panel=='b'))

ari.a <- plot.ari(filter(paths.com, panel=='a'))
ari.b <- plot.ari(filter(paths.com, panel=='b'))

geo.a <- plot.geo(filter(paths.com, panel=='a'))
geo.b <- plot.geo(filter(paths.com, panel=='b'), factor=3.1)


# Save plots
# ----------
ggsave(paste0(outpath, 'S03_a_complementarity_tigr_geo.pdf'),  geo.a,
  device=cairo_pdf, width=out.width, height=out.height)
ggsave(paste0(outpath, 'S03_a_complementarity_tigr_ari.pdf'),  ari.a,
  device=cairo_pdf, width=out.width, height=out.height)
ggsave(paste0(outpath, 'S03_a_complementarity_tigr_com.pdf'),  com.a,
  device=cairo_pdf, width=out.width, height=out.height)
ggsave(paste0(outpath, 'S03_b_complementarity_tigr_geo.pdf'),  geo.b,
  device=cairo_pdf, width=out.width, height=out.height)
ggsave(paste0(outpath, 'S03_b_complementarity_tigr_ari.pdf'),  ari.b,
  device=cairo_pdf, width=out.width, height=out.height)
ggsave(paste0(outpath, 'S03_b_complementarity_tigr_com.pdf'),  com.b,
  device=cairo_pdf, width=out.width, height=out.height)