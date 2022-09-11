# Dominance GWAS demo

## To simulate data

The script `R/simulate_toy_population.R` creates fake data
and outputs it into the `simulations` folder (where there
already is an example simulated dataset).


## To run GWAS

The `R/toy_dominance_gwas.R` runs genome scans on the
simulated data and makes some comparisons between the
additive only and dominance+additive model, and
with the simulated true variants.


## Functions for GWAS

The helper functions to run GWAS with hglm are in 
`R/hglm_functions.R`.
