#!/usr/bin/env Rscript
# 03_cedar_creek/01_fit_cedar_creek.R
# ===================================
# Author: Joe Wan
# Fits the Cedar Creek data to the MCT model and outputs fit information.
# Produces no figures but is required for the next scripts in the analysis.

# Data manipulation libraries
library(dplyr)
library(tidyr)
library(tibble)
library(magrittr)

# Plotting libraries
library(ggplot2)
library(cowplot)


##### 0. Output options #####
source('analysis_params.R')


##### Part 1: Processing the monoculture data #####
convert.biomass <- function(biomass.raw) exp(log(404)*biomass.raw)-1

mono.info.raw <- read.csv('inputs/1_monocultures.csv', stringsAsFactors=F)
mono.info <- mono.info.raw %>% mutate(biomass=convert.biomass(biomass.raw))
mono.plot.explore <- ggplot(mono.info, aes(x=soil.nitrogen, y=biomass, color=species)) +
  geom_point(aes()) +
  geom_smooth(method='gam', formula=y~s(x, bs="cs", k=5),
              se=F, fullrange=T,
              data=group_by(mono.info, species)) +
  expand_limits(x=0, y=0) +
  theme_cowplot()
if (show.plots) print(mono.plot.explore)

# Estimate Rstar using a linear regression
mono.smallest <- mono.info %>% group_by(species) %>% arrange(soil.nitrogen) %>% slice(1:3) %>% ungroup
Rstar.estimates <- c()
for (sp in unique(mono.smallest$species)) {
  df <- filter(mono.smallest, species==sp)
  model <- lm(biomass~soil.nitrogen, df)
  Rstar.estimates[sp] <- -model$coefficients[[1]]/model$coefficients[[2]]
}

# Use linear regression to estimate coefficients
mono.params.rough <- data.frame(Rstar=Rstar.estimates, a0=NA, a1=NA)
for (sp in unique(mono.smallest$species)) {
  Rstar <- mono.params.rough[sp,]$Rstar
  df <- filter(mono.info, species==sp) %>%
    mutate(x=(soil.nitrogen-Rstar), y=(soil.nitrogen-Rstar)/biomass)
  # In this transformed model, biomass~(soil.nitrogen-Rstar)/(a1*(soil.nitrogen-Rstar)+a0) becomes
  # (soil.nitrogen-Rstar)/biomass ~ (a0 + a1*(soil.nitrogen-Rstar))
  model <- lm(y~x, df)
  mono.params.rough[sp,c('a0','a1')] <- model$coefficients
}
mono.params.rough <- rownames_to_column(mono.params.rough, 'species')
rownames(mono.params.rough) <- mono.params.rough$species

# Now use nls to get higher-quality fits
mono.params.refined <- data.frame(Rstar=c(), Rstar.se=c(), 
                         a0=c(), a0.se=c(), 
                         a1=c(), a1.se=c(),
                         nonlin.bic=c(), lin.bic=c())
param.names <- c('Rstar','a1','a0')
mono.info.filtered <- filter(mono.info, species != 'Agrostis' | !(floor(soil.nitrogen) %in% c(245,654)))
for (sp in unique(mono.smallest$species)) {
  df <- filter(mono.info.filtered, species==sp)
  model <- nls(biomass~(soil.nitrogen-Rstar)/(a1*soil.nitrogen+a0), df,
               start=mono.params.rough[sp,param.names])
  mono.params.refined[sp,param.names] <- coef(model)[param.names]
  mono.params.refined[sp,paste0(param.names,'.se')] <- summary(model)$coefficients[param.names,2]
  mono.params.refined[sp,'nonlin.bic'] <- BIC(model)
}

# Now compare a linear model
lin.param.names <- c('Rstar','a0')
for (sp in unique(mono.smallest$species)) {
  df <- filter(mono.info.filtered, species==sp)
  model <- nls(biomass~(soil.nitrogen-Rstar)/(a0), df,
               start=mono.params.rough[sp,lin.param.names])
  mono.params.refined[sp, 'lin.bic'] <- BIC(model)
  if (sp != 'Agrostis') { # BIC(model)<mono.params.refined[sp, 'nonlin.bic']) {
    mono.params.refined[sp, param.names] <- 0
    mono.params.refined[sp, paste0(param.names,'.se')] <- NA
    mono.params.refined[sp, lin.param.names] <- coef(model)[lin.param.names]
    mono.params.refined[sp, paste0(lin.param.names,'.se')] <- summary(model)$coefficients[lin.param.names,2]
  }
}
mono.params.refined <- rownames_to_column(mono.params.refined, 'species')


##### Part 2: Calculating interspecific coefficients #####
comp.files <- c('2a'='inputs/2a_Agropyron_Schizachyrium.csv',
           '4a'='inputs/4a_Poa_Schizachyrium.csv',
           '6a'='inputs/6a_Agropyron_Poa.csv',
           '8a'='inputs/8a_Agropyron_Agrostis.csv')

# We're only interested in fitting the Agropyron-Poa pair (the only one to show 
# evidence of coexisting in the original study)
cur.raw <- read.csv(comp.files['6a'], stringsAsFactors=F)
cur.clean <- cur.raw %>% mutate(biomass=convert.biomass(biomass.raw))

# Identify the species and assign them S1 or S2
cur.spp <- unique(cur.clean$species)
cur.pair <- paste(range(cur.spp), collapse=' vs. ')
# Flip these so species 1 is Poa (more productive)
cur.s1 <- max(cur.spp)
cur.s2 <- min(cur.spp)

# We want pairs of biomasses from the same plot; this can be reconstructed since the
# plots have the same x values for corresponding points
by.plot <- cur.clean %>% group_by(species) %>%
  mutate(plot=row_number()) %>% ungroup
cur.plot.info <- group_by(by.plot, plot) %>%
  summarize(mean.soil.nitrogen=mean(soil.nitrogen),
            total.biomass=sum(biomass))

# Clean the original data
cur.clean <- left_join(by.plot, select(cur.plot.info, plot, mean.soil.nitrogen), by='plot') %>%
  transmute(pair=cur.pair, plot, soil.nitrogen=mean.soil.nitrogen, species, biomass) %>%
  mutate(species.id=ifelse(species==cur.s1, 1, 2))

# Extract the intraspecies parameters:
# "f11" and "g11" are the intercept and slope components, respectively
par1 <- filter(mono.params.refined, species==cur.s1)
f11 <- par1$a0; g11 <- par1$a1; Rstar1 <- par1$Rstar
par2 <- filter(mono.params.refined, species==cur.s2)
f22 <- par2$a0; g22 <- par2$a1; Rstar2 <- par2$Rstar

input <- cur.clean %>%
  mutate(x=soil.nitrogen,xs1=soil.nitrogen-Rstar1,xs2=soil.nitrogen-Rstar2,
                alpha11=f11+g11*x,alpha22=f22+g22*x) 

get.biomass <- function(x, species.id, params) {
  with(params, {
    # Calculate terms that don't change during optimization
    xs1 <- x-Rstar1; xs2 <- x-Rstar2
    alpha11 <- f11+g11*x; alpha22 <- f22+g22*x
    # Calculate terms that do change
    alpha12 <- f12+g12*x; alpha21 <- f21+g21*x
    denom <- alpha11*alpha22 - alpha21*alpha12
    # These two terms determine coexistence or exclusion
    num1 <- xs1*alpha22 - xs2*alpha12
    num2 <- xs2*alpha11 - xs1*alpha21
    # These are the monoculture equilibria
    mono1 <- xs1/alpha11
    mono2 <- xs2/alpha22
    # Determine which equilibrium applies and return it
    result <- case_when(
      num1<0 & num2<0 ~ as.double(NA),
      num1>0 & num2<0 ~ if_else(species.id==1, mono1, 0),
      num1<0 & num2>0 ~ if_else(species.id==2, mono2, 0),
      TRUE ~ if_else(species.id==1, num1, num2)/denom
    )
    return(result)
  })
}

# First round: constrain f using g so that persistence boundary happens where we want
bound1 <- 0; bound2 <- 1250 # Flipped
g.guess <- c(g12=0.001, g21=0.001)
model.est <- nls(biomass ~ get.biomass(soil.nitrogen, species.id,
                                   list(g12=g12, g21=g21,
                                        f12=(bound1-Rstar1)*(f22+g22*bound1)/(bound1-Rstar2)-g12*bound1,
                                        f21=(bound2-Rstar2)*(f11+g11*bound2)/(bound2-Rstar1)-g21*bound2)),
             input, start=g.guess,
             control=nls.control(maxiter=1000, minFactor=1e-9, warnOnly=T), trace=T)
comp.params.rough <- list(g12=coef(model.est)[['g12']], g21=coef(model.est)[['g21']],
                   f12=(bound1-Rstar1)*(f22+g22*bound1)/(bound1-Rstar2)-coef(model.est)[['g12']]*bound1,
                   f21=(bound2-Rstar2)*(f11+g11*bound2)/(bound2-Rstar1)-coef(model.est)[['g21']]*bound2)

fit.df.rough <- expand.grid(soil.nitrogen=nitrogen.values, species.id=c(1,2)) %>%
  mutate(x=soil.nitrogen,xs1=soil.nitrogen-Rstar1,xs2=soil.nitrogen-Rstar2,
         alpha11=f11+g11*x,alpha22=f22+g22*x,
         alpha12.est=comp.params.rough$f12+comp.params.rough$g12*x,alpha21.est=comp.params.rough$f21+comp.params.rough$g21*x) %>% select(-x) %>%
  mutate(biomass.est=get.biomass(soil.nitrogen, species.id, comp.params.rough))
comp.plot.rough <- ggplot(fit.df.rough, aes(x=soil.nitrogen, color=as.factor(species.id))) +
  geom_line(aes(y=biomass.est)) + geom_point(aes(y=biomass), data=input) +
  ylim(-0.1,200) +
  theme_cowplot()
if (show.plots) print(comp.plot.rough)

# Second round: starting from the previously-found parameters, do an unconstrained optimization
model.refined <- nls(biomass ~ get.biomass(soil.nitrogen, species.id,
                                       list(g12=g12, g21=g21,
                                            f12=(bound1-Rstar1)*(f22+g22*bound1)/(bound1-Rstar2)-g12*bound1,
                                            f21=(bound2-Rstar2)*(f11+g11*bound2)/(bound2-Rstar1)-g21*bound2)),
                 input, start=c(g12=comp.params.rough[['g12']], g21=comp.params.rough[['g21']],
                                bound1=bound1, bound2=bound2),
                 control=nls.control(maxiter=1000, minFactor=1e-9, warnOnly=T), trace=T)
comp.params.refined <- with(as.list(coef(model.refined)),
                   list(g12=g12, g21=g21,
                        f12=(bound1-Rstar1)*(f22+g22*bound1)/(bound1-Rstar2)-g12*bound1,
                        f21=(bound2-Rstar2)*(f11+g11*bound2)/(bound2-Rstar1)-g21*bound2))

# Calculate fits
fit.df.refined <- as.data.frame(fit.df.rough) %>% 
  select(soil.nitrogen, species.id, alpha11, alpha22) %>%
  mutate(x=soil.nitrogen,xs1=soil.nitrogen-Rstar1,xs2=soil.nitrogen-Rstar2,
         alpha12=comp.params.refined$f12+comp.params.refined$g12*x,alpha21=comp.params.refined$f21+comp.params.refined$g21*x) %>%
  mutate(biomass.est=get.biomass(soil.nitrogen, species.id, comp.params.refined)) %>%
  mutate(biomass.mono=ifelse(species.id==1, xs1/alpha11, xs2/alpha22)) %>%
  mutate(species=ifelse(species.id==1,cur.s1,cur.s2))
fit.df.compare <- bind_rows(
  mutate(fit.df.rough, pass='rough'),
  mutate(fit.df.refined, pass='refined'))

comp.plot.compare <- ggplot(fit.df.compare, aes(x=soil.nitrogen, color=as.factor(species.id))) +
  geom_line(aes(y=biomass.est, linetype=pass)) + 
  geom_point(aes(y=biomass), data=input) +
  ylim(-0.1,200) +
  theme_cowplot()
if (show.plots) print(comp.plot.compare)

# Export cleaned data
# -------------------
saveRDS(mono.info, paste0(outpath, '08_monoculture_data.rds'))
write.csv(mono.info, paste0(outpath, '08_monoculture_data.csv'), row.names=F)
saveRDS(cur.clean, paste0(outpath, '08_competition_data.rds'))
write.csv(cur.clean, paste0(outpath, '08_competition_data.csv'), row.names=F)


# Export parameters and fits
# --------------------------
# The raw params (Rstar and f, g)
all.params <- c(
  list(Rstar1=Rstar1, Rstar2=Rstar2,
    f11=f11, g11=g11, f22=f22, g22=g22), 
  comp.params.refined)
print(as.data.frame(all.params))
saveRDS(all.params, paste0(outpath, '08_raw_params.rds'))
# Dump them as a bare-bones CSV
as.data.frame(all.params) %>% t() %>%
  set_colnames(c('value')) %>% as.data.frame %>% 
  rownames_to_column('parameter') %>%
  write.csv(paste0(outpath, '08_raw_params.csv'), row.names=F, quote=F)

# Will be useful to work directly with the more familiar parameters
with(all.params, {
  param.df <- data.frame(soil.nitrogen=nitrogen.values)
  param.df <- mutate(param.df,
    Rstar1=Rstar1, Rstar2=Rstar2,
    a11=f11 + g11*soil.nitrogen, a22=f22 + g22*soil.nitrogen,
    a12=f12 + g12*soil.nitrogen, a21=f21 + g21*soil.nitrogen,
    alpha11=a11/(soil.nitrogen-Rstar1), alpha22=a22/(soil.nitrogen-Rstar2),
    alpha12=a12/(soil.nitrogen-Rstar1), alpha21=a21/(soil.nitrogen-Rstar2))
  saveRDS(param.df, paste0(outpath, '08_param_df.rds'))
  write.csv(param.df, paste0(outpath, '08_param_df.csv'),
    row.names=F, quote=F)
})
