
## Simulate a toy population to prototype estimation code on

library(AlphaSimR)


n_chr <- 5


## Create founder genomes with Macs

founders <- runMacs2(nInd = 200,
                       nChr = n_chr,
                       histGen = NULL,
                       histNe = NULL,
                       Ne = 5000)

simparam <- SimParam$new(founders)



## Pick common QTL. Do this per chromsome in order to find
## the within-chromsome index (lociLoc) rather than
## the global index

threshold <- 0.2

qtl_ix <- numeric(n_chr)
f_qtl <- numeric(n_chr)


for (chr in 1:n_chr) {
  
  geno <- pullSegSiteGeno(founders, chr = chr, simParam = simparam)
  
  f <- colSums(geno)/2/nrow(geno)
  
  candidate_qtl <- which(f > threshold & f < 1 - threshold)
  
  qtl_ix[chr] <- sample(candidate_qtl, 1)
  f_qtl[chr] <- f[qtl_ix[chr]]

}


## Set true a and d effects

a <- c(1, 1, 1, 0, 0)
d <- rep(1, 5)


## Add this trait to the simulation, setting a high heritability

trait <- new("TraitAD",
             intercept = 0,
             nLoci = as.integer(n_chr),
             lociPerChr = rep(1L, n_chr),
             lociLoc = as.integer(qtl_ix),
             addEff = a,
             domEff = d)

provisional_genetic_parameters <- AlphaSimR:::calcGenParam(trait, founders, 1)

resulting_additive_genetic_variance <- popVar(provisional_genetic_parameters$bv)[1]

trait@intercept <- - mean(provisional_genetic_parameters$gv)

simparam$manAddTrait(lociMap = trait, varE = resulting_additive_genetic_variance/3)


simparam$setSexes("yes_sys")

simparam$addSnpChip(nSnpPerChr = 5000)


## Create parent population

pop <- newPop(founders, simParam = simparam)


## Cross them to create mapping population

mapping_pop <- randCross(pop,
                         nProgeny = 10,
                         nCrosses = 100,
                         simParam = simparam)


## Pull out genotypes, phenotypes and metadata

geno <- pullSnpGeno(mapping_pop,
                    simParam = simparam)


marker_pos <- getSnpMap(simParam = simparam)

qtl_pos <- getQtlMap(simParam = simparam)

qtl_pos$a <- a
qtl_pos$d <- d
qtl_pos$f <- f_qtl

qtl_pos$varA <- 2 * qtl_pos$f * (1 - qtl_pos$f) * (qtl_pos$a + qtl_pos$d * (1 - 2 * qtl_pos$f))^2

qtl_pos$phen_var_expl <- qtl_pos$varA/varP(mapping_pop)[1]



## For each chromosome, find the correlations between the SNPs and the
## causative variant on that chromosome

snp_correlations <- vector(mode = "list", length = n_chr)

for (chr in 1:n_chr) {
  
  snps_on_chr <- pullSnpGeno(mapping_pop, simParam = simparam, chr = chr)
  
  qtl_on_chr <- pullQtlGeno(mapping_pop, simParam = simparam, chr = chr)
  
  snp_correlations[[chr]] <- as.vector(cor(snps_on_chr, qtl_on_chr))
  
}

marker_pos$qtn_cor <- unlist(snp_correlations)



## Save data

dir.create("simulations")

write.table(geno,
            file = "simulations/geno.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)


write.table(mapping_pop@pheno,
            file = "simulations/pheno.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)


write.table(marker_pos,
            file = "simulations/marker_positions.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)


write.table(qtl_pos,
            file = "simulations/qtl_positions.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
