#!/usr/bin/env bash
# 02_resource_model/00_run_resource_model.sh
# ==========================================
# Author: Joe Wan
# Runs all scripts in the `02_resource_model` directory.

# Setup for running all scripts
# ==============================
# Create log directory
log_dir="logs"
mkdir -p "$log_dir"
mkdir -p "outputs/"

# Define a function to run scripts
rule1="==========================="
rule2="---------------------------"
function run_script {
  echo "Running: $1"
  echo $rule1
  set -o pipefail # Capture exit status of Rscript
  Rscript $1 2>&1 | tee "$log_dir"/$(basename $1 .R).log && \
    echo $rule2$'\nSTATUS: Completed successfully'$'\n'$rule2 || \
    echo $rule2$'\nSTATUS: Failed'$'\n'$rule2
  echo $'\n\n\n\n'
}


# Run scripts using our function
# ==============================
run_script 01_resource_model.R