#!/usr/bin/env Rscript
# 04_multispecies/01_multispecies.R
# ==================================
# Author: Joe Wan
# Generates all figures for the multispecies model (Figures 7a-c, S9.1-5).


# Data manipulation libraries
library(tidyr)
library(dplyr)
library(stringr)
library(tibble)
library(broom)

# deSolve for ODE simulation
library(deSolve)

# Plotting libraries
library(ggplot2)
library(cowplot)


##### 0. Helper functions, etc. #####
# Load functions from utility file
source("species_pool_utils.R")


##### 1a. Define general parameters ######
# General parameters
n.sp <- 20        # Number of species to simulate
R0 <- 100         # Initial resource concentration
# Metaparameters for random generation of traits
base.cv <- 0.1    # Coefficient of variation for all traits
mu.mean <- 0.5    # Mean mortality rate
up.mean <- 2      # Mean resource uptake ability
    # Note: this was increased somewhat from previous versions to make beta a bit more important
eff.mean <- 0.05  # Mean resource use efficiency
B.inter.mean <- 1 # Mean interspecific interference
B.intra.mean <- 5 # Mean intraspecific interference
metaparams <- list(base.cv=base.cv, mu.mean=mu.mean, up.mean=up.mean, 
    eff.mean=eff.mean, B.inter.mean=B.inter.mean, B.intra.mean=B.intra.mean)

##### 1b. Generate the system using a defined random seed #####
# Set the seed
set.seed(20240204)
# Generate the traits with the chosen seed
traits <- get.species.traits(metaparams, n.sp)
# Get parameters from traits
params <- get.model.params(traits, R0)
# Append R0 to the traits
traits$R0 <- R0


##### 2a. Set up back-calculation visualization #####
# We'll vary NFD by changing the pairwise metrics, then back-calculate model parameters. 
# Here are the basic options for this analysis:

# Use the original NFDs and Ks from the model defined above
nfds <- get.nfds(traits)
Ks <- get.Ks(traits)
# Label options
ref.symbol <- "★"   # Symbol denoting the reference parameters
normal.symbol <- "" # "⦁" # Symbol denoting modified, back-calculated parameters
# Ranges for each of our analyses

# The following are helpers for this analysis:
# --------------------------------------------
get.process.df <- function(ds, nfds, Ks, transform.function) {
    df <- NULL
    for (i in 1:length(ds)) {
        factor <- exp(ds[i])
        result <- transform.function(nfds$rhos, nfds$frs, Ks, ds[i])
        new.rhos <- result$rhos; new.frs <- result$frs; new.Ks <- result$Ks
        # Perform the recalculation
        new.alphas <- backcalculate.alphas(new.rhos, new.frs, new.Ks)
        B.recalc <- get.B(new.alphas, traits$R0, traits$mus, traits$ups, traits$effs)
        # Stick it into a new traits list
        new.traits <- traits; new.traits$B <- B.recalc
        # Add to the data frame
        new.row <- tibble_row(
            index=i, d=ds[i], factor=factor, 
            alphas=list(new.alphas),
            rhos=list(new.rhos), 
            frs=list(new.frs),
            Ks=list(new.Ks), 
            traits=list(new.traits))
        if (is.null(df)) df <- new.row
        else df <- bind_rows(df, new.row) 
    }
    return(df)
}

# Sort and add species ranks for SAD plotting
add.sad.ranks <- function(df) {
    return(df %>% 
        group_by(index) %>%
        arrange(desc(N)) %>%
        mutate(rank=row_number())
    )
}

# Pick out the representative indices for the SADs
choose.indices <- function(num.sads, ds, pad=0, pad.orig=0) {
    if (pad>0) return(choose.indices(num.sads, ds[(1+pad):(length(ds)-pad)], pad=0, pad.orig=pad.orig) + pad)
    orig.index <- which(abs(ds)==min(abs(ds)))[1]
    max.index <- length(ds)
    # Original index should be included; calculate which-th of the selected panels it will be
    orig.pos <- round(num.sads*(orig.index-1)/(max.index-1))
    # If the original index is too close to the edge, we'll need to pad it
    orig.pos <- max(pad.orig+1, min(num.sads-pad.orig, orig.pos))
    # First and last species will be (orig.pos - 1) and (num.sads - orig.pos)
    # spots before and after the original index.
    spots.before <- orig.pos - 1
    spots.after <- num.sads - orig.pos
    skip <- min(floor((max.index-orig.index)/spots.after),
                floor((orig.index-1)/spots.before))
    return(seq(
        orig.index - spots.before*skip, 
        orig.index + spots.after*skip, 
        by=skip))
}

# Pick out example points
get.examples.df <- function(df, rev=F, num.sads=3, pad=1, pad.orig=0) {
    sad.indices <- choose.indices(num.sads, df$d, pad=pad, pad.orig=pad.orig)
    if (rev) sad.indices <- rev(sad.indices)
    ref.index <- which(abs(df$d)==min(abs(df$d)))[1]
    examples <- data.frame(
            id=1:length(sad.indices), index=sad.indices) %>%
        left_join(df, by="index") %>%
        mutate(is.ref=index==ref.index) 
    return(examples)
}

# Convert the example data frame to have one row per species
get.example.sp.df <- function(df) {
    return(df %>%
        get.species.df() %>% add.sad.ranks() %>%
        ungroup %>% as.data.frame %>%
        arrange(id) %>%
        mutate(facet.title=factor.preserve.order(facet.title),
            facet.title.short=factor.preserve.order(facet.title.short))
    )
}


##### 2b. Simulate scenarios (Figure 7) #####
# Panel (a): Stabilizing 
# ----------------------
# Stabilizing i.e. decreasing niche overlap corresponds to increasing rho. 
# To do so, we will add some factor d to log(rho_ij), equivalent to 
# multiplying by exp(d).

# Define the ds used in the analysis
rho.ds <- seq(-0.15, 0.6, by=0.01)

# Simulate the scenarios
cat("Simulating stabilizing scenarios...\n")
stabilizing.df <- get.process.df(rho.ds, nfds, Ks, 
        function(rhos, frs, Ks, d) {
            # New rhos (equivalent to multiplying by exp(d))
            new.rhos <- exp(log(rhos) + d)
            # Reset diagonals to one
            for (j in 1:length(Ks)) new.rhos[j,j] <- 1
            return(list(rhos=new.rhos, frs=frs, Ks=Ks))
    }) %>%
    # Simulate to equilibrium and add summary stats
    add.simulated.equilibria(time=100, verbose=F, refine=T) %>% 
    add.equilibrium.summary()

# Get example points for the main figure
stabilizing.examples <- stabilizing.df %>%
    get.examples.df(rev=T, num.sads=3, pad=1, pad.orig=0) %>%
    mutate(facet.title.short=paste0(
        if_else(is.ref, ref.symbol, normal.symbol), " ", id),
        facet.title=paste0(facet.title.short, "\n",
        "1−ρ = ", sprintf("%.2f", 1-median.rho)))
stabilizing.examples.sp <- stabilizing.examples %>%
    get.example.sp.df()

# Panel (b): (De)equalization 
# ---------------------------
# We want to increase fitness in favor of the more productive species
# To do so, we will add some factor sgn(K_i-K_j)*d to log(f_i/f_j),  
# equivalent to multiplying/dividing by exp(d)

# Define the ds used in the analysis
fr.ds <- seq(-0.9, 0.6, by=0.02)

# Simulate the scenarios
cat("Simulating equalizing scenarios...\n")
equalizing.df <- get.process.df(fr.ds, nfds, Ks, 
        function(rhos, frs, Ks, d) {
            n.sp <- length(Ks)
            signs <-  sign(diag(Ks) %*% ones(n.sp) - ones(n.sp) %*% diag(Ks))
            new.frs <- exp(log(frs) + signs * d)
            # Signs will be 0 along the diagonals, so new diagonals should remain one
            return(list(rhos=rhos, frs=new.frs, Ks=Ks))
    }) %>%
    # Simulate to equilibrium and add summary stats
    add.simulated.equilibria(time=100, verbose=F, refine=T) %>% 
    add.equilibrium.summary() # %>%
    # add.simulated.timeseries(seq(0, 100, by=0.1))

# Get example points for the main figure
equalizing.examples <- equalizing.df %>%
    get.examples.df(rev=F, num.sads=3, pad=1, pad.orig=1) %>%
    mutate(facet.title.short=paste0(
        if_else(is.ref, ref.symbol, normal.symbol), " ", id),
        facet.title=paste0(facet.title.short, "\n",
        "f₁/f₂ = ", sprintf("%.2f", median.ordered.fr)))
equalizing.examples.sp <- equalizing.examples %>%
    get.example.sp.df()

# Panel (c): Functional (productivity) equalization
# -------------------------------------------------
# Let's increase or decrease variation in productivity. To do 
# so, we will scale the distribution of Ks around some reference value
# by exp(d). There are multiple possibilities; let's use the most productive 
# species for the main figure (since it matches the two-species case) and 
# also examine the effect of scaling around the mean for the supplement.

# Define the ds used in the analysis
K.ds <- seq(-1.5, 1.5, length=76)

# Simulate the scenarios
cat("Simulating productivity equalization scenarios...\n")
prodeq.df <- get.process.df(K.ds, nfds, Ks, 
        function(rhos, frs, Ks, d) {
            ref <- max(Ks)
            new.Ks <- exp(log(ref) + (log(Ks) - log(ref)) * exp(d))
            return(list(rhos=rhos, frs=frs, Ks=new.Ks))
    }) %>%
    # Simulate to equilibrium and add summary stats
    add.simulated.equilibria(time=100, verbose=F, refine=T) %>% 
    add.equilibrium.summary() # %>%
    # add.simulated.timeseries(seq(0, 100, by=0.1))

# Get example points for the main figure
prodeq.examples <- prodeq.df %>%
    get.examples.df(rev=T, num.sads=3, pad=5, pad.orig=1) %>%
    mutate(
        facet.title.short=paste0(
            if_else(is.ref, ref.symbol, normal.symbol), " ", id),
        facet.title=paste0(facet.title.short, "\n",
            "Kₘᵢₙ = ", sprintf("%.2f", min.K) #,
            # "(cv = ", sprintf("%.2f", cv.K), ")"
        ))
prodeq.examples.sp <- prodeq.examples %>%
    get.example.sp.df()

##### 3a. Set up plotting #####
# Define visuals
oranges <- c('#ffc27d', '#ffa154', '#f37329', '#cc3b02', '#a62100')
greens <- c('#d1ff82', '#9bdb4d', '#68b723', '#3a9104', '#206b00')
blues <- c('#8cd5ff', '#64baff', '#3689e6', '#0d52bf', '#002e99')
coex.col <- 'gray65'
colored.line.width <- 0.75

# Labels with math
nd.label <- expression(paste('niche difference, ', 1 - rho))
median.nd.label <- expression(paste('median niche difference, ', 1 - rho))
# nd.label <- "niche difference, 1 − ρ"
fd.label <- expression(paste('fitness difference, ', log(f[1] / f[2])))
# fd.label <- "fitness difference, log(f₁ / f₂)"
fr.label <- expression(paste('fitness ratio, ', f[1] / f[2]))
median.fr.label <- expression(paste('median fitness ratio, ', f[1] / f[2]))
# fr.label <- "fitness ratio, f₁ / f₂"
K.label <- expression(paste('minimum intrinsic yield, ', K[min]))
# K.label <- "minimum yield, Kₘᵢₙ"
K.diff.label <- expression(paste('yield range, ', K[max] - K[min]))
# K.diff.label <- "yield range, Kₘₐₓ - Kₘᵢₙ"
resource.label <- "resource lvl., R"
abundance.label <- expression(paste("abund., ", N[i])) # "Nᵢ"

# Output options
show.plots <- interactive() # Will suppress plot output if not REPL
outpath <- 'outputs/'       # Directory to save plots
out.width <- 4              # Width of main plots
out.height <- 3             # Height of main plots
# Font, etc. will be added by this theme:
th <- theme(text=element_text(family="Helvetica"))

# Make outpath if it doesn't exists
if (!str_ends(outpath, '/')) outpath <- paste0(outpath, '/')
dir.create(outpath, showWarnings=F)

# Nudge factor for text labels
nudge.factor <- 0.15

# Axis titles for re-use
biomass.label <- expression(paste("total biomass, ", Sigma * N))
biomass.label.short <- expression(paste("tot. biomass, ", Sigma * N))
stabilization.label <- "Stabilization"
equalization.label <- "Equalization (fitness--function)"
prodeq.label <- "Functional equalization"
# Unicode combining character (poor alignment): n*N̂

# Helper functions for plotting:
# ------------------------------
# Main plot
plot.process <- function(df, example.df, aesthetic, color="black", ref.color=NULL,
        show.min=T, show.max=T, show.ribbon=T,
        nudge=0.125, expand.x=NULL, expand.y=NULL, xlim=NULL, ylim=NULL) {
    if (is.null(ref.color)) ref.color <- color
    result <- ggplot(df, aesthetic) +
        geom_path(color=color, linewidth=colored.line.width)
    # Dotted lines for max (and optionally min) yield
    if (show.min) result <- result + geom_line(aes(y=min.K), linetype='dotted')
    if (show.max) result <- result + geom_line(aes(y=max.K), linetype='dotted')
    if (show.ribbon) result <- result +
        # Add ribbon highlighting overyield
        geom_ribbon(aes(ymin=max.K, ymax=ifelse(total.N>max.K, total.N, NA)), 
            fill=greens[2], alpha=0.1)
    result <- result +
        # Add text labels below and horizontally centered
        geom_text(aes(label=id), 
            nudge_y=nudge, hjust=0.5, 
            data=example.df) +
        geom_point(data=filter(example.df, !is.ref), 
            color=color) +
        # Add ref as five pointed star
        geom_text(data=filter(example.df, is.ref), 
            color=ref.color, 
            label='★', size=6)
    if (!is.null(expand.x)) result <- result + expand_limits(x=expand.x)
    if (!is.null(expand.y)) result <- result + expand_limits(y=expand.y)
    result <- result + 
        coord_cartesian(expand=F, xlim=xlim, ylim=ylim) +
        theme_cowplot()
    return(result)
}

# SAD plot
plot.sads <- function(examples.sp, color="black", alpha=1, 
        max.x=NULL, max.y=NULL, expand.y=NULL, shared.y=T) {
    if (is.null(max.x)) max.x <- max(examples.sp$richness)
    if (is.null(max.y)) max.y <- max(examples.sp$N)
    yl <- c(0, max.y)
    if (!shared.y) yl <- NULL
    result <- ggplot(examples.sp) +
        geom_bar(aes(x=rank, y=N), stat="identity",
            fill=color, alpha=alpha) +
        facet_wrap(~facet.title, nrow=1, scales="free") +
        theme_cowplot() +
        coord_cartesian(
            xlim=c(0, max.x) + 0.5, 
            ylim=yl,
            expand=FALSE) +
        xlab("species rank") + ylab("abund.") +
        # Remove gray bar behind facet labels
        theme(strip.background=element_blank())
    if (is.null(expand.y) && !shared.y) expand.y <- c(0)
    else if (!shared.y) expand.y <- c(expand.y, 0)
    if (!is.null(expand.y)) result <- result + expand_limits(y=expand.y)
    return(result)
}


##### 3b. Plot the scenarios (Figure 7) #####
# Panel (a): Stabilizing
# ----------------------
# Main plot
stabilizing.y <- c(2.85, 4.15) # Need this later
stabilizing.df.trimmed <- stabilizing.df %>% 
    # Chop off some of the right and top, so the undershoot is better visible
    filter(median.rho > 0.52)
    # 1 - rho < 0.49 ==> rhos > 0.52
stabilizing.plot <- stabilizing.df.trimmed %>%
    plot.process(stabilizing.examples, 
        aes(x=1-median.rho, y=total.N), color=blues[2], ref.color=NULL, # Used "black" but it looked awkward
        show.min=F, expand.y=stabilizing.y, nudge=nudge.factor*1) +
        xlab(median.nd.label) + ylab(biomass.label)
if (show.plots) print(stabilizing.plot)
# SAD plot
stabilizing.sads <- plot.sads(stabilizing.examples.sp, color=blues[2], max.x=16)
if (show.plots) print(stabilizing.sads)

# Panel (b): Equalizing
# ---------------------
# Main plot
optimal.fitness <- equalizing.df %>% arrange(desc(total.N)) %>% head(1) %>% pull(median.ordered.fr)
equalizing.plot <- plot.process(equalizing.df, equalizing.examples, 
    aes(x=median.ordered.fr, y=total.N), color=oranges[3], ref.color=NULL, # Used "black" but it looked awkward
    show.min=T, ylim=c(1.25, 4.25), nudge=nudge.factor*-2.25) +
    # Add dashed line for optimal fitness
    geom_vline(xintercept=optimal.fitness, linetype='dashed') +
    # Make sure 1 appears on the x axis
    scale_x_continuous(breaks=seq(0.2, 2.2, by=0.4)) +
    xlab(median.fr.label) + ylab(biomass.label) # fᵢ/fⱼ
if (show.plots) print(equalizing.plot)
# SAD plot
equalizing.sads <- plot.sads(equalizing.examples.sp, color=oranges[3], max.x=16)
if (show.plots) print(equalizing.sads)

# Panel (c): Productivity equalization
# ------------------------------------
# Main plot
prodeq.plot <- plot.process(prodeq.df, prodeq.examples, 
    aes(x=min.K, y=total.N), color=greens[3], ref.color=NULL, # Used "black" but it looked awkward
    show.min=T, ylim=c(1.95, 4.95), nudge=nudge.factor*-2.5) +
    xlab(K.label) + ylab(biomass.label)
if (show.plots) print(prodeq.plot)
# SAD plot
prodeq.sads <- plot.sads(prodeq.examples.sp, color=greens[3], max.x=16, max.y=1) +
        scale_y_continuous(breaks=seq(0,1,by=1))
if (show.plots) print(prodeq.sads)


##### 3. Export plots #####
if (!str_ends(outpath, '/')) outpath <- paste0(outpath, '/')
dir.create(outpath)

main.width <- out.width
main.height <- out.height
sads.width <- out.width
sads.height <- out.height*2/3

sad.extra <- scale_x_continuous(breaks=seq(5,15,by=5))

# All the panels have the same height, so we can keep them in proportion
# by exporting with the same dimensions
ggsave(paste0(outpath, '07_a_stabilizing.pdf'), 
    stabilizing.plot + th,
    device=cairo_pdf, width=main.width, height=main.height)
ggsave(paste0(outpath, '07_a_stabilizing_sads.pdf'), 
    stabilizing.sads + sad.extra + th,
    device=cairo_pdf, width=sads.width, height=sads.height)
ggsave(paste0(outpath, '07_b_equalizing.pdf'), 
    equalizing.plot + th,
    device=cairo_pdf, width=main.width, height=main.height)
ggsave(paste0(outpath, '07_b_equalizing_sads.pdf'), 
    equalizing.sads + sad.extra + th,
    device=cairo_pdf, width=sads.width, height=sads.height)
ggsave(paste0(outpath, '07_c_prodeq.pdf'), 
    prodeq.plot + th,
    device=cairo_pdf, width=main.width, height=main.height)
ggsave(paste0(outpath, '07_c_prodeq_sads.pdf'), 
    prodeq.sads + sad.extra + th,
    device=cairo_pdf, width=sads.width, height=sads.height)


##### 4. Plots for appendix #####
# Export a zoomed-in version of the stabilizing plot (Figure S9.2)
# ----------------------------------------------------------------
# Define desired aspect ratio of the new plot and its x range
new.aspect <- 1/3
zoom.x <- c(0.09, 0.19)
zoom.y.min <- 2.9


# Aspect and x range of the old plot
old.aspect <- main.height/main.width
stabilizing.x <- range(1 - stabilizing.df.trimmed$median.rho)
# Shortcut to calculate the difference of a range
d <- function(a) a[2]-a[1]
# Calculate the ratio of coordinate units for coord_fixed
heights <- rep(c(7/9, 3/9), 3)
fix.ratio <- heights[1]*old.aspect*d(stabilizing.x)/d(stabilizing.y)
# Calculate the new plot's y range according to specifications above
zoom.y <- zoom.y.min  + c(0, new.aspect*d(zoom.x)/fix.ratio)

# Create the zoomed-in plot
stabilizing.plot.zoom <- stabilizing.plot +
    coord_fixed(ratio=fix.ratio, xlim=zoom.x, ylim=zoom.y, expand=F) +
    scale_y_continuous(breaks=seq(2.5, 3.5, by=0.1)) +
    scale_x_continuous(breaks=seq(0, 1, by=0.025))

# Save the zoomed-in plot
ggsave(paste0(outpath, 'S9_02_stabilizing_zoom.pdf'),
    stabilizing.plot.zoom,
    device=cairo_pdf, width=main.width, height=main.width*new.aspect)


# Bonus analysis for productivity equalization (Figure S9.4)
# ----------------------------------------------------------
# In the main figure, we scaled the productivity distribution around its
# maximum (i.e. keeping the most productive species constant). This raises
# the question of whether the results are simply because most species increased
# in productivity. To address this, we'll also examine the effect of scaling
# the distribution around its center.

# Define the ds used in the analysis
K.alt.ds <- seq(-0.25, 1.125, length=76)
# Note: values higher than this cause numerical issues

# Simulate the scenarios
cat("Simulating alternative productivity equalization scenarios...\n")
prodeq.alt.df <- get.process.df(K.alt.ds, nfds, Ks, 
        function(rhos, frs, Ks, d) {
            ref <- exp(mean(log(Ks)))
            new.Ks <- exp(log(ref) + (log(Ks) - log(ref)) * exp(d))
            return(list(rhos=rhos, frs=frs, Ks=new.Ks))
    }) %>%
    # Simulate to equilibrium and add summary stats
    add.simulated.equilibria(time=100, verbose=F, refine=T) %>% 
    add.equilibrium.summary()

# Get example points for the main figure
prodeq.alt.examples <- prodeq.alt.df %>%
    get.examples.df(rev=F, num.sads=3, pad=5, pad.orig=0) %>%
    mutate(facet.title.short=paste0(
        if_else(is.ref, ref.symbol, normal.symbol), " ", id),
        facet.title=paste0("\n",
        "Kₘₐₓ - Kₘᵢₙ = ", sprintf("%.2f", min.K) #,
        # "(cv = ", sprintf("%.2f", cv.K), ")"
        ))
prodeq.alt.examples.sp <- prodeq.alt.examples %>%
    get.example.sp.df()

# Plot, as above:
prodeq.alt.plot <- plot.process(prodeq.alt.df, prodeq.alt.examples, 
    aes(x=max.K-min.K, y=total.N), color=greens[3], 
    show.min=T, 
    nudge=0.3
    # ylim=c(1.95, 4.95), nudge=0.225
    ) +
    xlab(K.diff.label) + ylab(biomass.label)
if (show.plots) print(prodeq.alt.plot)
prodeq.alt.sads <- plot.sads(prodeq.alt.examples.sp, color=greens[3], max.x=16)
if (show.plots) print(prodeq.alt.sads)

# Save, as above
ggsave(paste0(outpath, 'S9_04_prodeq_alt.pdf'), 
    prodeq.alt.plot + th,
    device=cairo_pdf, width=main.width, height=main.height)
ggsave(paste0(outpath, 'S9_04_prodeq_alt_sads.pdf'), 
    prodeq.alt.sads + sad.extra + 
        scale_y_continuous(breaks=seq(0,1,by=1)) + 
        th,
    device=cairo_pdf, width=sads.width, height=sads.height)


# Plot histograms of niche, fitness, and productivity (Figure S9.1)
# -----------------------------------------------------------------
# Helper functions to plot histograms
plot.histogram <- function(v, bins=20, xlim=NULL, ylim=NULL) {
    plt <- data.frame(x=v) %>%
        ggplot() +
            geom_histogram(aes(x=x), bins=bins, alpha=0.4) +
            coord_cartesian(xlim=xlim, ylim=ylim, expand=FALSE) +
            ylab("# species pairs") +
            theme_cowplot()
    return(plt)
}

nd.histogram <- plot.histogram(
        1-get.off.diags(nfds$rhos, triangular=T), bins=13,
        xlim=c(0.385, 0.515), ylim=c(0, 35)) +
    xlab(nd.label)
if (show.plots) print(nd.histogram)

fd.histogram <- plot.histogram(
        get.ordered.frs(nfds$frs, Ks), bins=13,
        xlim=c(0.8, 1.4+1e-8), ylim=c(0, 35)) +
    geom_vline(xintercept=1, linetype='dashed') +
    xlab(fd.label)
if (show.plots) print(fd.histogram)

K.histogram <- plot.histogram(Ks, bins=9,
    xlim=c(0, 3), ylim=c(0, 7)) +
    xlab("intrinsic yield, K")
if (show.plots) print(K.histogram)

# Save as above
ggsave(paste0(outpath, 'S9_01_a_nd_histogram.pdf'), 
    nd.histogram, 
    device=cairo_pdf, width=main.width, height=main.height*7/9)
ggsave(paste0(outpath, 'S9_01_b_fd_histogram.pdf'),
    fd.histogram, 
    device=cairo_pdf, width=main.width, height=main.height*7/9)
ggsave(paste0(outpath, 'S9_01_c_K_histogram.pdf'),
    K.histogram, 
    device=cairo_pdf, width=main.width, height=main.height*7/9)


# Plot pairwise NFD for equalizing scenario (Figure S9.3)
# -------------------------------------------------------
# Calculate pairwise metrics
pairwise <- tibble()
for (i in 1:n.sp) {
    for (j in 1:n.sp) {
        if (Ks[i]>Ks[j]) {
            new.row <- tibble(
                K1=Ks[i], K2=Ks[j], 
                fr=nfds$frs[i,j], 
                rho=nfds$rhos[i,j])
            pairwise <- bind_rows(new.row, pairwise)
        }
    }
}

# Calculate fitness ratios for key events
eq.coex.fr.range <- equalizing.df %>% filter(richness>1) %>% 
    pull(median.ordered.fr) %>% range
eq.toy.fr.range <- equalizing.df %>% filter(total.N>max.K) %>% 
    pull(median.ordered.fr) %>% range
eq.ref.fr <- equalizing.df %>% filter(d==0) %>%
    pull(median.ordered.fr)
eq.opt.fr <- equalizing.df %>% arrange(desc(total.N)) %>% head(1) %>% 
    pull(median.ordered.fr)

# Calculate regions for (excess) NFD conditions
regions <- data.frame(rho=seq(0.5,1.1,by=0.001))
frs <- unique(c(eq.ref.fr, eq.opt.fr, eq.coex.fr.range, eq.toy.fr.range))
pairwise.examples <- data.frame(median.ordered.fr=frs) %>%
    arrange(median.ordered.fr) %>%
    mutate(id=row_number()) %>%
    left_join(equalizing.df, by="median.ordered.fr") %>%
    mutate(is.ref=d==0,
        facet.title=paste0("atop('", 
            paste(if_else(is.ref, ref.symbol, normal.symbol), id), 
            "', f[1]/f[2] == '",  sprintf("%.2f", median.ordered.fr), "')")) %>%
    mutate(facet.title=factor(facet.title,      levels=facet.title))
conditions <- c("transgressive overyielding", "coexistence")

# Plot faceted plot
pairwise.fct.plot <- pairwise %>%
    mutate(shift=sqrt(K1/K2)) %>%
    expand_grid(pairwise.examples) %>%
    mutate(factor=median.ordered.fr/eq.ref.fr) %>%
    expand_grid(condition=conditions) %>%
    # Coexistence compares -log rho and -log f1/f2,
    # while TOY divides these by sqrt(K1/K2)
    mutate(shift=if_else(
        condition=="transgressive overyielding", 
        shift, 1)) %>%
    mutate(x=log(1/rho/shift)) %>%
ggplot() +
    # Coexistence or TOY requirement zone
    geom_ribbon(aes(x=x, ymin=x, ymax=-x, 
            fill=condition), 
        data=expand_grid(x=seq(0,0.8,by=0.01), condition=conditions)) +
    geom_abline(intercept=0, slope=c(1,-1),             linetype='solid') +
    geom_hline(yintercept=0, linetype="dotted") +
    geom_point(aes(x=x, y=log(factor*fr/shift)),
        alpha=1, size=0.1) +
    coord_equal(ratio=1, 
        expand=F) +
    facet_grid(paste0("'", condition, "'")~facet.title, labeller=label_parsed) +   
    xlab(expression(paste(
        "(excess) niche difference, ", 
        - "log" ~ rho - Delta))) +
    ylab(expression(atop(
        "(excess) fitness difference, ", 
        "log" ~ f[1]/f[2] - Delta))) +
    scale_x_continuous(breaks=seq(0,0.5,by=0.5)) +
    scale_fill_manual(values=c(coexistence=coex.col, 
        `transgressive overyielding`=greens[2])) +
    guides(fill=guide_none()) +
    theme_cowplot() +
    theme(strip.background=element_blank(), strip.clip="off")
if (show.plots) print(pairwise.fct.plot)

# Re-plot the process plot with extra examples
pairwise.process.plot <- plot.process(equalizing.df, pairwise.examples, 
    aes(x=median.ordered.fr, y=total.N), color=oranges[3], ref.color=NULL,
    show.min=T, ylim=c(1.25, 4.25), nudge=nudge.factor*-1.5) +
    # Make sure 1 appears on the x axis
    scale_x_continuous(breaks=seq(0.2, 2.2, by=0.4)) +
    xlab(median.fr.label) + ylab(biomass.label) # fᵢ/fⱼ
if (show.plots) print(pairwise.process.plot)

# Save, making things a bit wider
pairwise.fac <- 1.5
pairwise.heights <- c(2/3, 1.1)
ggsave(paste0(outpath, 'S9_03_equalizing_pairwise_process.pdf'), pairwise.process.plot, 
    device=cairo_pdf, width=main.width*pairwise.fac, height=main.height*pairwise.fac*sum(pairwise.heights))
ggsave(paste0(outpath, 'S9_03_equalizing_pairwise_nfd.pdf'), pairwise.fct.plot, 
    device=cairo_pdf, width=main.width*pairwise.fac, height=main.height*pairwise.fac*sum(pairwise.heights))


# Pick out time series for the reference model and plot them separately
# ---------------------------------------------------------------------
# Numbers of species to simulate
example.ns <- c(1,2,5,20)
# Time series length
example.timeseries.length <- 100

# Simulate the time series
example.ts <- stabilizing.examples %>% 
    # Grab the reference system only
    filter(d==0) %>%
    # Use join on dummy variable to expand to all numbers of species
    transmute(traits=traits, dummy=T) %>%
    right_join(data.frame(n.sp=example.ns, dummy=T), by="dummy") %>% 
    select(-dummy) %>%
    # Add traits and simulate
    rowwise %>% mutate(traits=list(subset.system(traits, 1:n.sp))) %>% 
    add.simulated.timeseries(seq(0, example.timeseries.length, by=0.1)) %>%
    # Get one row for each species at each time point
    unnest(timeseries) %>% 
    # Add the resource use efficiency, for plotting
    left_join(data.frame(species=1:n.sp, eff=traits$effs), by="species")

# Calculate resource level, etc.
example.ts.R <- example.ts %>% 
    group_by(time, n.sp) %>%
    summarize(consumed=sum(N / eff), Ntot=sum(N)) %>%
    mutate(R=traits$R0-consumed) %>%
    ungroup

# Grapical options
ts.line.width <- 0.5
 # Theme for all panels
tts <- th +
    # Increase margin on the right
    theme(plot.margin=margin(t=7, b=7, l=7, r=14))

# Plot the time series by looping through number of species
example.ts.plots <- list()
example.ts.R.plots <- list()
example.ts.composites <- list()
for (n in unique(example.ts$n.sp)) {
    # Subset needed data
    cur.ts <- filter(example.ts, n.sp==n)
    cur.R <- filter(example.ts.R, n.sp==n)
    max.n <- max(cur.ts$N)

    # Plot the time series 
    example.ts.plots[[n]] <- cur.ts %>%
        ggplot(aes(x=time)) +
        geom_line(aes(y=N, color=as.factor(species)), linewidth=ts.line.width) +
        expand_limits(y=c(0, 1, max.n + 0.01)) +
        scale_y_continuous(breaks=seq(0,10,by=0.5)) +
        coord_cartesian(expand=F) +
        ylab(abundance.label) +
        theme_cowplot() + tts +
        theme(legend.position="none")

    # Plot resource level
    example.ts.R.plots[[n]] <- cur.R %>%
        ggplot(aes(x=time)) +
        geom_line(aes(y=R), linewidth=ts.line.width) +
        expand_limits(y=c(0, traits$R0)) +
        coord_cartesian(expand=F) +
        ylab(resource.label) +
        theme_cowplot() + tts +
        theme(legend.position="none")
    
    # Add more lines to plots depending on number of species
    if (n > 5) example.ts.R.plots[[n]] <- example.ts.R.plots[[n]] +
        geom_line(aes(y=consumed), 
        linewidth=ts.line.width, linetype="dashed")
    # For fewer, we can show each species' consumption
    else example.ts.R.plots[[n]] <- example.ts.R.plots[[n]] +
        geom_line(aes(y=N/eff, color=as.factor(species)), 
        linewidth=ts.line.width, linetype="dashed",
            data=cur.ts)
    
    # Save
    ggsave(paste0(outpath, 'S6_01_example_timeseries_n=', n, '.pdf'), 
        example.ts.plots[[n]],
        device=cairo_pdf, width=main.width, height=main.height*0.75)
    ggsave(paste0(outpath, 'S6_01_example_timeseries_n=', n, '_R.pdf'), 
        example.ts.R.plots[[n]],
        device=cairo_pdf, width=main.width, height=main.height*0.75)
}


# Analyze invasion growth rate (Figure S9.5)
# ------------------------------------------
# Helper function: get the initial abundances for the invasion
get.n0s <- function(n.sp, i, n0=1e-5, focal=0) {
    n0s <- rep(n0, n.sp)
    n0s[i] <- 0
    return(n0s)
}
# Helper function: get the equilibrium abundance of species i
del.sp <- function(equilibrium, i) {
    result <- equilibrium
    # Remove the species
    result[i] <- 0
    return(result)
}

# Calculate equilibrium for the reference model
ref.df <- data.frame(x=NA) %>%
    rowwise %>%
    mutate(traits=list(traits)) %>%
    mutate(equilibrium=list(
        get.equilibrium(resource.model, get.model.params(traits), resource.growth)
    )) %>%
    select(-x)
ref.eq <- ref.df$equilibrium[[1]]

# Helper function: remove species i from the community
get.subset <- function(n.sp, i) {
    result <- 1:n.sp
    result <- result[-i]
    return(result)
}
# Helper function: restore species i into the equilibrium
expand.equilibrium <- function(equilibrium, i) {
    # Expand the equilibrium to include the invader
    if (i == 1) head <- c()
    else head <- equilibrium[1:(i-1)]
    if (i == length(equilibrium)+1) tail <- c()
    else tail <- equilibrium[i:length(equilibrium)]
    return(c(head, 0, tail))
}

# Perform the invasion analysis
invasion.df <- data.frame(i=1:n.sp) %>%
    rowwise %>%
    mutate(sub.traits=list(
        subset.system(traits, get.subset(n.sp, i))
    )) %>% mutate(sub.equilibrium=list(
        get.equilibrium(resource.model, get.model.params(sub.traits), resource.growth)
    )) %>% mutate(equilibrium=list(
        expand.equilibrium(sub.equilibrium, i)
    )) %>% mutate(
        igr=resource.growth(equilibrium, get.model.params(traits))[i],
        r=resource.growth(rep(0, n.sp), get.model.params(traits))[i],
        K=get.Ks(traits)[i]
    ) %>% ungroup

# Calculate the equilibrium abundance of the invader and add its IGR
full.df <- ref.df %>%
    crossing(i=1:n.sp) %>%
    rowwise %>%
    mutate(N=equilibrium[i]) %>% 
    ungroup %>%
    left_join(invasion.df, by="i")

# Estimate the slope of the N ~ FK relationship
FK.slope <- full.df %>% filter(igr>0) %>%
    transmute(FK=K*igr/r, N) %>%
    lm(N ~ 0 + FK, data=.) %>%
    tidy() %>%
    filter(term=='FK') %>% pull(estimate)
cat("Slope of N ~ FK relationship: ", FK.slope, "\n")

# Plot the N ~ FK relationship
FK.plot <- ggplot(full.df) +
    geom_point(aes(x=igr/r * K, y=N)) +
    # Use pmax to get the "dog-leg" or "hockey stick" shape
    geom_line(aes(x=x, y=pmax(0, FK.slope*x)), 
        data=data.frame(x=seq(-0.05,0.1, len=1000)),
        linetype='dashed') +
    xlab(expression(paste('fitness' %*% 'intrinsic yield, ', F[i]~K[i]))) + 
    ylab(expression(paste('equilbrium abundance, ', N[i]))) +
    scale_x_continuous(expand=c(0,0)) +
    theme_cowplot() +
    theme(legend.position="none")
if (show.plots) print(FK.plot)

# Save the plot
ggsave(paste0(outpath, 'S9_05_invasion_growth.pdf'), 
    FK.plot + th, 
    device=cairo_pdf, width=main.width, height=main.height)
