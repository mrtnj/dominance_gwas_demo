
## Run a dominance GWAS on toy data

library(dplyr)
library(hglm)
library(readr)
library(patchwork)


source("R/hglm_functions.R")



## Read in simulated data and true positions

geno <- read_tsv("simulations/geno.txt",
                 col_names = FALSE)

pheno <- read_tsv("simulations/pheno.txt",
                  col_names = FALSE)

marker_pos <- read_tsv("simulations/marker_positions.txt")

marker_pos$marker_id <- paste("X", marker_pos$id, sep = "")

qtl_pos <- read_tsv("simulations/qtl_positions.txt")


## Remove monomorphic SNPs (the simulated data has no missingness, but
## SNPs are occasionally monomorphic in the mapping population)

geno_pruned <- prune_snp_matrix(geno)



## Construct GRM; in a real analysis this would probably need to be replaced with
## an externally computed GRM made with some more efficient methods

G <- construct_grm(geno_pruned)


## Decomposition of the GRM to create an incidence matrix for the additive
## genetic effect

Z_grm <- decompose_grm(G, FALSE)


## Create a design matrix for fixed effects (intercept only)

X <- matrix(1, ncol = 1, nrow = dim(geno))


## Number of colums for each random effect; here only one

RandC <- ncol(Z_grm)


## GWAS with and without dominance term for SNP effects

gwasA <- run_gwasA(pheno = pheno$X1,
                   X = X,
                   Z = Z_grm,
                   RandC = RandC,
                   geno_pruned)


gwasAD <- run_gwasAD(pheno = pheno$X1,
                     X = X,
                     Z = Z_grm,
                     RandC = RandC,
                     geno_pruned)


## Add positions to results

gwasA_pos <- inner_join(gwasA, marker_pos)

gwasAD_pos <- inner_join(gwasAD, marker_pos)


## Plot the estimates, with the position of the true causative variant as
## red triangles

plot_true_estimated_posA <- qplot(x = pos, y = -log10(p), data = gwasA_pos) +
  facet_wrap(~chr) +
  geom_point(aes(x = pos, y = 15), colour = "red", data = qtl_pos, shape = 2)

plot_true_estimated_posAD <- qplot(x = pos, y = -log10(p), data = gwasAD_pos) +
  facet_wrap(~chr) +
  geom_point(aes(x = pos, y = 15), colour = "red", data = qtl_pos, shape = 2)

plot_true_estimated <- plot_true_estimated_posA / plot_true_estimated_posAD


## Look at significant positions in the two scans

filter(gwasA_pos, p < 1e-8)

filter(gwasAD_pos, p < 1e-8)
