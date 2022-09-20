
## Remove missing and monomorphic SNPs

prune_snp_matrix <- function(snp_matrix) {
  
  missing <- lapply(snp_matrix, function(x) sum(is.na(x)))
  
  genotype_values <- lapply(snp_matrix, function(x) length(unique(x)))
  
  snp_matrix[, missing == 0 & genotype_values > 1]
  
}


## Only remove monomorphic -- for data with missing genotypes

remove_monomorphic <- function(snp_matrix) {

  genotype_values <- lapply(snp_matrix, function(x) length(unique(na.exclude(x))))
  
  snp_matrix[, genotype_values > 1]
  
}



## Decompose the gnomic relationship matrix to find the modified incidence
## matrix

decompose_grm <- function(grm, missing) {
  svd <- svd(grm[!missing, !missing])
  svd$u %*% diag(sqrt(svd$d))
}



## Construct genomic relationship matrix from genotypes,
## centring and scaling by the sum of allele frequency variances

construct_grm <- function(geno) {
  
  f <- colSums(geno)/2/nrow(geno)
  
  P <- matrix(2 * (f - 0.5), 
              byrow = TRUE,
              ncol = ncol(geno),
              nrow = nrow(geno))
  
  Z <- as.matrix(geno - P)
  
  scaling <- 2 * sum(f *  (1 - f))
  
  G <- Z %*% t(Z) / scaling
  
  G
}


## Create the V matrix for fitting SNP effects efficiently

make_V_matrix <- function(Z, RandC, ratio) {
  n_obs <- nrow(Z)
  n_components <- length(ratio)
  V <- diag(n_obs)
  
  start <- c(1, cumsum(RandC)[-length(RandC)] + 1)
  end <- cumsum(RandC)
  for (ix in 1:n_components) {
    Z_component <- Z[, start[ix]:end[ix]]
    V <- V + tcrossprod(Z_component) * ratio[ix]
  }
  V
}



## Run GWAS with additive and dominance coding

run_gwasAD <- function(pheno,
                       X,
                       Z,
                       RandC,
                       snp_matrix) {
  
  n_ind <- length(pheno)
  n_snp <- ncol(snp_matrix)
  
  ## Fit baseline model
  
  model_baseline <- hglm(y = pheno,
                         X = X,
                         Z = Z,
                         RandC = RandC)
  
  
  ratio <- model_baseline$varRanef/model_baseline$varFix
  
  V <- make_V_matrix(Z,
                     RandC,
                     ratio)
  
  eigV <- eigen(V)
  
  transformation_matrix <- diag(1/sqrt(eigV$values)) %*% t(eigV$vectors)
  
  transformed_y <- transformation_matrix %*% pheno
  transformed_X <- transformation_matrix %*% X
  
  
  ## Null model
  
  qr0 <- qr(transformed_X)
  
  est0 <- qr.coef(qr0, transformed_y)
  
  null_residual <- transformed_y - transformed_X %*% est0
  
  RSS_null <- sum(null_residual^2)/n_ind
  
  
  ## SNP models
  
  estimate_a <- numeric(n_snp)
  estimate_d <- numeric(n_snp)
  LRT <- numeric(n_snp)
  p <- numeric(n_snp)
  
  for (snp_ix in 1:n_snp) {
    
    missing_genotype <- is.na(snp_matrix[, snp_ix])
    
    transformed_snp <- transformation_matrix[!missing_genotype, !missing_genotype] %*%
      as.matrix(snp_matrix[!missing_genotype, snp_ix] - 1)
    transformed_d <- transformation_matrix[!missing_genotype, !missing_genotype] %*%
      as.matrix(as.numeric(snp_matrix[!missing_genotype, snp_ix] == 1))
    
    X1 <- cbind(transformed_snp, transformed_d, transformed_X[!missing_genotype,])
    qr1 <- qr(X1)
    est1 <- qr.coef(qr1, transformed_y[!missing_genotype])
    residual1 <- transformed_y[!missing_genotype] - X1 %*% est1
    RSS1 <- sum(residual1^2)/(n_ind - sum(missing_genotype))
    
    estimate_a[snp_ix] <- est1[1]
    estimate_d[snp_ix] <- est1[2]
    LRT[snp_ix] <- -(n_ind - sum(missing_genotype)) * (log(RSS1) - log(RSS_null))
    p[snp_ix] <- 1 - pchisq(LRT[snp_ix],
                            df = 1)
  }
  
  data.frame(marker_id = colnames(snp_matrix),
             estimate_a,
             estimate_d,
             LRT,
             p)
}


## Run GWAS with additive only coding

run_gwasA <- function(pheno,
                      X,
                      Z,
                      RandC,
                      snp_matrix) {
  
  n_ind <- length(pheno)
  n_snp <- ncol(snp_matrix)
  
  ## Fit baseline model
  
  model_baseline <- hglm(y = pheno,
                         X = X,
                         Z = Z,
                         RandC = RandC)
  
  
  ratio <- model_baseline$varRanef/model_baseline$varFix
  
  V <- make_V_matrix(Z,
                     RandC,
                     ratio)
  
  eigV <- eigen(V)
  
  transformation_matrix <- diag(1/sqrt(eigV$values)) %*% t(eigV$vectors)
  
  transformed_y <- transformation_matrix %*% pheno
  transformed_X <- transformation_matrix %*% X
  
  
  ## Null model
  
  qr0 <- qr(transformed_X)
  
  est0 <- qr.coef(qr0, transformed_y)
  
  null_residual <- transformed_y - transformed_X %*% est0
  
  RSS_null <- sum(null_residual^2)/n_ind
  
  
  ## SNP models
  
  estimate <- numeric(n_snp)
  LRT <- numeric(n_snp)
  p <- numeric(n_snp)
  
  for (snp_ix in 1:n_snp) {
    
    missing_genotype <- is.na(snp_matrix[, snp_ix])
    
    transformed_snp <- transformation_matrix[!missing_genotype,!missing_genotype] %*%
      as.matrix(snp_matrix[!missing_genotype, snp_ix] - 1)
    
    X1 <- cbind(transformed_snp, transformed_X[!missing_genotype,])
    qr1 <- qr(X1)
    est1 <- qr.coef(qr1, transformed_y[!missing_genotype])
    residual1 <- transformed_y[!missing_genotype] - X1 %*% est1
    RSS1 <- sum(residual1^2)/(n_ind - sum(missing_genotype))
    
    estimate[snp_ix] <- est1[1]
    LRT[snp_ix] <- -(n_ind - sum(missing_genotype)) * (log(RSS1) - log(RSS_null))
    p[snp_ix] <- 1 - pchisq(LRT[snp_ix],
                            df = 1)
  }
  
  data.frame(marker_id = colnames(snp_matrix),
             estimate,
             LRT,
             p)
}


## Get heritability estimate from a hgml model

get_h2 <- function(model) {
  model$varRanef[1] / (model$varFix + sum(model$varRanef))
}
