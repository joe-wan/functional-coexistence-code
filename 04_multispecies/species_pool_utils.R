#!/usr/bin/env Rscript
# 04_multispecies/species_pool_utils.R
# ====================================
# Author: Joe Wan
# Helper functions for working with the multispecies model.

# Data manipulation libraries
library(tidyr)
library(dplyr)
library(stringr)
library(tibble)

# deSolve for ODE simulation
library(deSolve)


##### 0. General helpers for calculations #####
# "Clip" abundances to zero
clip <- function(N) pmax(N, 0)

# Generate a matrix of ones
ones <- function(n.sp) matrix(rep(1, n.sp^2), nrow=n.sp, ncol=n.sp) 

# Useful function for reading off the terms of interest (i.e. off-diagonals)
get.off.diags <- function(mat, triangular=F, upper=T) {
    n.sp <- dim(mat)[1]
    is <- rep(seq(1, n.sp), n.sp)
    js <- rep(seq(1, n.sp), each=n.sp)
    if (triangular && upper) return(mat[is<js])
    if (triangular && !upper) return(mat[is>js])
    return(mat[is!=js])
}

# Get the fitness ratio with more productive species in the numerator
get.ordered.frs <- function(frs, Ks) {
    result <- list()
    for (i in 1:length(Ks)) {
        for (j in 1:length(Ks)) {
            if (i!=j && Ks[i] > Ks[j]) {
                result = append(result, frs[i,j])
            }
        }
    }
    return(as.numeric(result))
}

# Convert some vector to a vector of factors, preserving order
factor.preserve.order <- function(x) factor(x, levels=unique(x))


##### 1. Defining the ODE model and its traits #####
# Function implementing the derivatives for deSolve, ported from the Julia model:
# # dn[:] = n .* ((R_0 - dot(c, n))*r ./ (B*n .+ 1) .- mu)
resource.model <- function(t, N, params) {
    # Parameters:
    #     r: vector of utilization abilities (uptake * efficiency)
    #     c: vector of resource requirements
    #     mu: vector of mortalities
    #     B: matrix of interaction coefficients, where B[i,j] is
    #        beta_ij (effect of species j on i)
    #     R_0: the total resource in the system (scalar)
    dn <- clip(N) * ((params$R0 - sum(params$c * N)) * params$r / 
        (params$B %*% N + 1) - params$mu)
    # # Originally this was:
    # dn[:] = n .* ((R_0 - dot(c, n))*r ./ (B*n .+ 1) .- mu)
    # # but clipping was added to ensure no negative population blowup
    return(list(dn))
}

resource.growth <- function(N, params) {
    # Parameters:
    #     r: vector of utilization abilities
    #     c: vector of resource requirements
    #     mu: vector of mortalities
    #     B: matrix of interaction coefficients, where B[i,j] is
    #        beta_ij (effect of species j on i)
    #     R_0: the total resource in the system (scalar)
    return(((params$R0 - sum(params$c * N)) * params$r / 
        (params$B %*% N + 1) - params$mu))
}

# Helper functions to calculate parameters from desired moments
get.lognormal.params <- function(mean, sd=NULL, cv=NULL) {
    if (!is.null(sd) & !is.null(cv)) stop("Cannot specify both sd and cv")
    if (!is.null(cv)) return(get.lognormal.params(mean, cv*mean))
    return(list(mu=log(mean^2/sqrt(mean^2+sd^2)), 
        sigma=sqrt(log(1+sd^2/mean^2))))
}
get.beta.params <- function(mean, sd=NULL, cv=NULL) {
    if (!is.null(sd) & !is.null(cv)) stop("Cannot specify both sd and cv")
    if (!is.null(cv)) return(get.beta.params(mean, cv*mean))
    # The standard deviation cannot be too large
    if (sd >= sqrt(mean * (1-mean))) {
        stop("Beta distribution must satisfy sd < sqrt(mean * (1-mean))")
    }

    # Knowing that parameters are good, return them
    v <- mean*(1-mean)/sd^2 - 1
    return(list(alpha=mean*v, beta=(1-mean)*v))
}

# Generate species traits from metaparameters describing their distributions
get.species.traits <- function(metaparams, n.sp) {
    # Convert mean/cv and generate the vector-valued parameters
    mu.params <- get.lognormal.params(metaparams$mu.mean, cv=metaparams$base.cv)
    mus <- rlnorm(n.sp, meanlog=mu.params$mu, sdlog=mu.params$sigma)
    eff.params <- get.beta.params(metaparams$eff.mean, cv=metaparams$base.cv)
    effs <- rbeta(n.sp, shape1=eff.params$alpha, shape2=eff.params$beta)
    up.params <- get.lognormal.params(metaparams$up.mean, cv=metaparams$base.cv)
    ups <- rlnorm(n.sp, meanlog=up.params$mu, sdlog=up.params$sigma)

    # Generate matrices of mean and sds for B
    B.intra.params <- get.lognormal.params(metaparams$B.intra.mean, cv=base.cv)
    B.inter.params <- get.lognormal.params(metaparams$B.inter.mean, cv=base.cv)
    B.mu <- (1-diag(n.sp))*B.inter.params$mu + diag(n.sp)*B.intra.params$mu
    B.sigma <- (1-diag(n.sp))*B.inter.params$sigma + diag(n.sp)*B.intra.params$sigma
    B <- exp(B.mu + matrix(rnorm(n.sp^2), nrow=n.sp, ncol=n.sp)*B.sigma)
    # # Double check means
    # sum(B*diag(n.sp))/n.sp
    # sum(B*(1-diag(n.sp)))/(n.sp^2-n.sp)

    return(list(mus=mus, effs=effs, ups=ups, B=B))
}

get.model.params <- function(traits, R0=NULL) {
    # If R0 is not specified, try to get it from traits
    if (is.null(R0)) {
        if (!"R0" %in% names(traits)) stop("R0 must be in trait list if not specified")
        R0 <- traits$R0
    }
    # Calculate actual model parameters
    rs <- traits$up * traits$eff
    cs <- traits$eff^-1
    
    return(list(r=rs, cs=cs, mu=traits$mu, B=traits$B, R0=R0))
}

subset.system <- function(system, species) {
    # Iterate through the parameters in the named list
    result <- list()
    for (k in names(system)) {
        v <- system[[k]]
        if (is.matrix(v) && length(dim(v)) == 2) {
            new.v <- v[species, species]
            if (length(species)==1) new.v <- matrix(new.v, 1, 1)
        } else if (is.vector(v) && length(v) > 1) {
            new.v <- v[species]
        } else {
            new.v <- v
        }
        new.list <- list(new.v)
        names(new.list) <- k
        result <- c(result, new.list)
    }
    return(result)
}

##### 2a. Functions to get pairwise MCT metrics #####
summarize.traits <- function(params) {
    # Vector of Rstars: straightforwardly apply the formula
    #   Rstar = mu / (eff * up)
    Rstars <- params$mu / (params$eff * params$up)
    # Matrix of as: apply the formula
    #   a_ij = beta_ij * Rstar_i + eff_j^-1
    n.sp <- dim(params$B)[1]
    A <- diag(Rstars) %*% params$B + ones(n.sp) %*% diag(params$eff^-1)
    return(c(list(Rstars=Rstars, A=A), params))
}

get.pairwise.rho <- function(summary, i, j) {
    # Check that the input is a summary
    if(any(!(c('Rstars', 'A') %in% names(summary))))
        stop("Input must be the result of summarize.traits")
    
    # Implement the formula:
    #   rho_ij = sqrt(a_ij * a_ji / (a_ii * a_jj))
    A <- summary$A
    return(sqrt(A[i,j]*A[j,i]/(A[i,i]*A[j,j])))
}

get.rhos <- function(summary) {
    # If the input is the raw traits, recurse on the summary
    if(any(!(c('Rstars', 'A') %in% names(summary))))
        get.rho(summarize.traits(summary))
    # Else we have a summary, so can proceed.

    # Implement the formula:
    #   rho_ij = sqrt(a_ij * a_ji / (a_ii * a_jj))
    A <- summary$A
    self.as <- diag(A)
    return(sqrt(A * t(A) /
        (self.as %*% t(self.as))))
}

get.pairwise.fr <- function(summary, i, j) {
    # Check that the input is a summary
   if(any(!(c('Rstars', 'A') %in% names(summary))))
        stop("Input must be the result of summarize.traits")
    
    # Implement the formula:
    #                        _______________
    #             R - R*_i   | a_ji * a_jj 
    #   f_i/f_j = ________   | ____________
    #             R - R*_j  \|  a_ij * a_ii
    A <- summary$A
    R0 <- summary$R0
    Rstars <- summary$Rstars
    return((R0-Rstars[i])/(R0-Rstars[j]) * sqrt(A[j,i]*A[j,j]/(A[i,j]*A[i,i])))
}

get.frs <- function(summary) {
    # # If the input is the raw traits, recurse on the summary
    if(any(!(c('Rstars', 'A') %in% names(summary))))
        return(get.pairwise.rho(summarize.traits(summary)))
    # # Else we have a summary, so can proceed.
        
    # Implement the formula:
    #                        _______________
    #             R - R*_i   | a_ji * a_jj 
    #   f_i/f_j = ________   | ____________
    #             R - R*_j  \|  a_ij * a_ii
    A <- summary$A
    R0 <- summary$R0
    Rstars <- summary$Rstars
    n.sp <- dim(A)[1]
    extra.R <- R0-Rstars
    resource.part <- extra.R %*% t(extra.R^-1)
    # We need to match each a_ij with its corresponding a_ii
    products <- A * (diag(diag(A)) %*% ones(n.sp))
    root.part <- sqrt(t(products) * products^-1)
    return(resource.part * root.part)
}

get.single.K <- function(summary, i) {
    # Check that the input is a summary
    if(any(!(c('Rstars', 'A') %in% names(summary))))
        stop("Input must be the result of summarize.traits")
    
    # Implement the formula:
    #   K_i = (R_0 - R*_i) / a_ii
    return((summary$R0 - summary$Rstars[i])/summary$A[i,i])
}

get.Ks <- function(summary) {
    # # If the input is the raw traits, recurse on the summary
    if(any(!(c('Rstars', 'A') %in% names(summary))))
        return(get.Ks(summarize.traits(summary)))
    # # Else we have a summary, so can proceed.
    
    # Implement the formula:
    #   K_i = (R_0 - R*_i) / a_ii
    return((summary$R0 - summary$Rstars)/diag(summary$A))
}

get.nfds <- function(summary) {
    # # If the input is the raw traits, recurse on the summary
    if(any(!(c('Rstars', 'A') %in% names(summary))))
        return(get.nfds(summarize.traits(summary)))
    # # Else we have a summary, so can proceed.
    
    rhos <- get.rhos(summary)
    frs <- get.frs(summary)

    return(list(rhos=rhos, frs=frs))
}


##### 2b. Functions to back-calculate parameters from MCT metrics #####
backcalculate.alphas <- function(rhos, frs, Ks) {
    self.alphas <- Ks^-1
    # Implement the formula
    #   alpha_ij = rho / (f_i/f_j) * alpha_jj
    # We will need the alpha_jj for each alpha_ij i.e. a matrix with element [i,j] = alpha_jj.
    # This is obtained by multiplying a matrix of all ones by a diagonal matrix with the self-alphas.
    n.sp <- dim(rhos)[1]
    return(rhos * (frs^-1) * (ones(n.sp) %*% diag(self.alphas)))
    # Note this gives the correct result for alpha_ii, since rho = fr = 1 in those cases.
}

get.B <- function(alphas, R0, mus, ups, effs) {
    # Implement the formula
    #   beta_ij = (alpha_ij*(R0 - R*_i) - eff_j^-1) / (R0 - Rstar_i)
    # We will need the Rstar_i for each beta_ij i.e. a matrix with element [i,j] = Rstar_i,
    # and the eff_j for each beta_ij i.e. a matrix with element [i,j] = eff_j.
    # These are obtained by multiplying a matrix of all ones by a diagonal matrix 
    # in the correct order: diag(Rstars) %*% ones and ones %*% diag(effs), respectively:
    n.sp <- dim(alphas)[1]
    Rstars <- mus / (effs * ups)
    row.Rstars <- diag(Rstars) %*% ones(n.sp)
    col.effs <- ones(n.sp) %*% diag(effs)
    # Now the formula can straightforwardly be implemented as:
    return((alphas*(R0 - row.Rstars) - col.effs^-1) / row.Rstars)
}


##### 3. Functions for automatically detecting model equilbria #####
get.n.sp <- function(params) {
    # Find the first vector in params
    v <- NULL
    for (i in 1:length(params)) {
        if (is.vector(params[[i]])) {
            v <- params[[i]]
            break
        }
    }
    # Ensure that a vector was found
    if (is.null(v)) stop("No vector found in params")
    # Number of species is the number of elements in the vector
    return(length(v))
}

detect.zero.growth <- function(final.n, growth.function, params, threshold=1e-7)
    abs(growth.function(final.n, params)) < abs(threshold)
detect.is.present <- function(final.n, threshold=1e-7) final.n > abs(threshold)

at.equilibrium <- function(final.n, params, growth.function, growth.threshold=1e-7, abund.threshold=1e-7) {
    # Consider a species present if final abundance is above abs(abund.threshold)
    # Then the system is at equilibrium when all species that are present have
    # growths sufficiently near zero: i.e. magnitude less than abs(growth.threshold).

    is.present <- detect.is.present(final.n, threshold=abund.threshold)
    zero.growth <- detect.zero.growth(final.n, growth.function, params, threshold=growth.threshold)
    return(all(!is.present | zero.growth))
}

get.equilibrium <- function(fn, params, growth.function, ..., time=1e6, n0=1e-5, maxtime=1e12, refine=F, verbose=F) {
    n.sp <- get.n.sp(params)
    if (length(n0) > 1) n0s <- n0 
    else n0s <- rep(n0, n.sp)
    # ode()'s default solver is LSODA, which uses adaptive time steps, so no need for intermediate times
    times <- seq(0, time, len=1000)
    
    # Solve the ODE with deSolve
    if (verbose) print(paste("Solving ODE with time", time))
    out <- ode(..., y=n0s, times=times, func=fn, parms=params)
    # Grab the last row of the output
    final.n <- as.vector(out[nrow(out),-1]) # Remove the time column

    # Check if the system is at equilibrium
    if (at.equilibrium(final.n, params, growth.function)) {
        if (refine) {
            if (verbose) print("Refining results")
            new.start <- ifelse(detect.is.present(final.n) & detect.zero.growth(final.n, growth.function, params), final.n, 0)
            return(get.equilibrium(fn, params, growth.function, 
                time=time, n0=new.start, maxtime=maxtime, refine=F, verbose=verbose))
        }
        return(final.n)
    } else if (time*10 < maxtime) {
        if (verbose) print(paste("Increasing time to", 10*time))
        return(get.equilibrium(fn, params, growth.function, 
            time=10*time, n0=n0, maxtime=maxtime, refine=refine, verbose=verbose))
    } else {
        if (verbose) print("Reached maximum simulation time; returning NULL")
        return(NULL)
    }
}


##### 3b. Functions for organizing results #####
# Make data frame from deSolve output
format.df <- function(out) { 
    return(data.frame(out) %>%
        gather(key="species", value="N", -time) %>%
        mutate(species=str_replace(species, "X", "")) %>%
        mutate(species=as.numeric(species)))
}

add.simulated.equilibria <- function(df, ...) {
    result <- df %>% 
        rowwise() %>%
        mutate(eq=list(
            get.equilibrium(resource.model, get.model.params(traits), resource.growth, ...)
            )) %>%
        # Turn off rowwise
        ungroup()
    return(result)
}

add.equilibrium.summary <- function(df, na.rm=T, ...) {
    result <- df %>% 
        rowwise() %>%
        mutate(
            ...,
            median.rho=median(get.off.diags(rhos), na.rm=na.rm),
            median.ordered.fr=median(get.ordered.frs(frs, Ks), na.rm=na.rm),
            min.K=min(Ks, na.rm=na.rm), max.K=max(Ks, na.rm=na.rm),
            median.K=median(Ks, na.rm=na.rm),
            sd.K=sd(Ks, na.rm=na.rm), cv.K=sd.K/mean(Ks, na.rm=na.rm),
            mean.log.rho=mean(log(get.off.diags(rhos)), na.rm=na.rm),
            mean.ordered.log.fr=mean(log(get.ordered.frs(frs, Ks)), na.rm=na.rm),
            mean.log.K=mean(log(Ks), na.rm=na.rm),
            sd.log.K=sd(log(Ks), na.rm=na.rm), cv.log.K=sd.log.K/mean(log(Ks), na.rm=na.rm),
            total.N=sum(eq),
            richness=sum(eq>1e-7)) %>%
        # Turn off rowwise
        ungroup()
    return(result)
}

add.simulated.timeseries <- function(df, times, only.present=T, 
        threshold=1e-1, n0=1e-5, ...) {
    result <- df %>% 
        rowwise() %>%
        mutate(timeseries=list(format.df(
            ode(..., 
                y=rep(n0, dim(traits$B)[1]),
                times=times, 
                func=resource.model, 
                parms=get.model.params(traits))
            ))) %>%
        # Turn off rowwise
        ungroup()
    return(result)
}

get.species.df <- function(df) {
    result <- df %>% 
        rowwise() %>%
        mutate(eq=list(data.frame(species=1:n.sp, N=eq))) %>% 
        ungroup() %>%
        unnest(eq) %>% 
        select(-c(rhos, frs, Ks, traits, alphas))
    if ("timeseries" %in% names(df)) results <- select(result, -timeseries)
    return(result)
}