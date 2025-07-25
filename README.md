# Code for: "Functional coexistence theory: a mechanistic framework linking biodiversity to ecosystem function"

This repository contains the original analysis code for the following 
manuscript:

> Functional coexistence theory: a mechanistic framework linking biodiversity 
> to ecosystem function  
> Joe Wan*, Po-Ju Ke*, Iris Hordijk, Lalasia Bialic-Murphy, and Thomas W. 
> Crowther  
> Preprint (bioRxiv). doi: https://doi.org/10.1101/2024.05.05.591902


## File Structure
Each set of analyses is in its own folder, `0X_analysis_name/`, containing 
all scripts, files with parameters or 
helper functions (e.g., `.../graph_helpers.R`), and inputs (for the Cedar Creek 
analysis in `03_cedar_creek/`). In each directory the script 
`.../00_run_analysis_name.sh` handles the entire analysis. Each analysis is meant
to be run within its own working directory.

The contents of this repository are:
- `00_run_all.sh` runs analyses 1–5 in their respective working directories.
- `01_fct_simulations/` performs simulations for the general (i.e. 
   Lotka–Volterra) results, producing main text Figures 1-4, 9, and Supplemental 
   Figure S1.
- `02_resource_model/` handles the theoretical simulations of the resource 
   competition model for two species, producing main text Figure 7.
- `03_cedar_creek/` contains digitized results and analyses for Wedin and Tilman
   (1993)'s experiment at Cedar Creek using our resource competition model, 
   producing main text Figure 8 as well as Supplemental Figures S5-S7.
- `04_multispecies/` simulates the resource competition from above for more than
   two species, producing main text Figure 6 and Supplemental Figures S4, 
   S8-12.
- `05_tigr_complementarity/` compares geometric, arithmetic, and 
   complementarity-based metrics for niche and fitness difference, and produces 
   Supplemental Figure S3.


## Software
All analyses were conducted in [R](https://www.r-project.org/) (version 4.3.0) 
and [GNU bash](https://www.gnu.org/software/bash/) (version 5.0.17(1)-release). 
The following package versions were used for the R analyses:
> [broom](https://broom.tidymodels.org/) (v1.0.5), 
>[cowplot](https://wilkelab.org/cowplot/) (v1.1.1), 
>[deSolve](http://desolve.r-forge.r-project.org/) (v1.40), 
>[dplyr](https://dplyr.tidyverse.org) (v1.1.2), 
>[ggplot2](https://ggplot2.tidyverse.org) (v3.4.4), 
>[magrittr](https://magrittr.tidyverse.org) (v2.0.3), 
>[RColorBrewer](NA) (v1.1-3), 
>[scales](https://scales.r-lib.org) (v1.2.1), 
>[stringr](https://stringr.tidyverse.org) (v1.5.0), 
>[tibble](https://tibble.tidyverse.org/) (v3.2.1), 
>[tidyr](https://tidyr.tidyverse.org) (v1.3.0)


## Author Contributions
JW and IH designed the resource competition model, and JW and PJK designed the 
remaining analyses with input from all authors. JW wrote analysis code and 
derived the theoretical results.