# 01_fct_simulations/mct_params.R
# ===============================
# Author: Joe Wan
# Shared parameters for MCT and FCT simulations

# Bounds for ND (= -log(rho))
max.nd.lo <- 0.55 # Upper and lower limit for (origin-centered) MCT figure
max.nd <- -log(1-0.55) # 0.95 # Upper bound for FCT figure
min.nd <- max.nd-(2*max.nd.lo) # Calculate so figures have same aspect ratio

# Founds for FR (f1/f2 = exp(FD)), shared between MCT/FCT panels
min.fr <- 0.425
max.fr <- 1/0.425

# Calculate necessary aspect ratio
ratio <- (log(max.fr)-log(min.fr))/(max.nd-min.nd)
fct.x <- 1-exp(-c(min.nd, max.nd))
fct.y <- c(min.fr, max.fr)

# Calculate necessary aspect ratio for wider plot
max.nd.wide <- min.nd+log(max.fr)-log(min.fr)
ratio.wide <- (log(max.fr)-log(min.fr))/(max.nd.wide-min.nd)
fct.x.wide <- 1-exp(-c(min.nd, max.nd.wide))

# Increment for continuous axes
resolution <- 0.001

# Intrinsic yields
Y1 <- 1
Y2 <- 1/(1.2^2)