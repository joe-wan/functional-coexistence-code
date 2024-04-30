#!/usr/bin/env bash
# 00_run_all.sh
# =============
# Author: Joe Wan
# Runs all scripts in the functional-coexistence-code/ analysis repository.

# Prompt user to see if we should clean up previous outputs
# =========================================================
# If trash not installed, define it as rm -rf
read -p "Do you want to clean up previous outputs? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
  if ! command -v trash &> /dev/null; then
    read -p "Outputs will be permanently deleted. Continue? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
      (rm -rf outputs/ logs/ */outputs/ */logs/)
    fi
  else
    (trash outputs/ logs/ */outputs/ */logs/)
  fi
fi

# Setup for running all scripts
# ==============================
# Create log directory
log_dir="logs"
mkdir -p "$log_dir"

# Define a function to run scripts
rule1="==========================="
rule2="---------------------------"
function run_script {
  dir=$(dirname $1)
  script=$(basename $1)
  echo "Running: $script"
  echo "Working directory: $dir"
  echo $rule1
  set -o pipefail # Capture exit status of Rscript
  (cd "$dir" && bash "$script") 2>&1 | tee "$log_dir"/$(basename $script .R).log && \
    echo $rule2$'\nSTATUS: Script '$script$' completed successfully'$'\n'$rule2 || \
    echo $rule2$'\nSTATUS: Failed'$'\n'$rule2
  echo $'\n\n\n\n'
}


# Run scripts using our function
# ==============================
run_script 01_fct_simulations/00_run_fct_simulations.sh
run_script 02_resource_model/00_run_resource_model.sh
run_script 03_cedar_creek/00_run_cedar_creek.sh
run_script 04_multispecies/00_run_multispecies.sh
run_script 05_complementarity/00_run_complementarity.sh


# Put all outputs in a single directory
# =====================================
read -p "Do you want to collect all outputs in one directory? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
  mkdir -p outputs
  # Copy all outputs to the outputs directory
  for i in $(ls */outputs/*); do
    cp $i outputs/$(basename $i)
  done
fi
