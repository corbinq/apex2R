# -------------
# in bash
# -------------
## git clone https://github.com/corbinq/apex2R.git

# -------------
# in R
# -------------

setwd("apex2R")
source("Apex2R.r")
setwd("..")

ssm <- metaSumStats(
    G1000 ="G1000_DATA/OUTPUT/G1000_chr21",
    HapMap ="HAPMAP_DATA/OUTPUT/hapmap_chr21"
)

# problematic gene
gene <- "ENSG00000154645"

## Get XtY and XtX objects for specified gene
gene_sf <- ssm$getSuffStats(gene)

# SNPs that were selected based on 'apex meta --stepwise'
snps <- c("21_17025359_G_A", "21_17915436_G_A")
idx <- match(snps, colnames(gene_sf$XtX))

# get summary statistics for the SNPs of interest
xtx <- gene_sf$XtX[idx,idx]
xty <- gene_sf$Xty[idx]

# covariate-adjusted SNP covariance matrix
xtx

# covariate-adjusted SNP corr matrix
cov2cor(xtx)

# OLS/WLS slope estimates
beta <- solve(xtx, xty)

sigma_hat <- sqrt((gene_sf$yty - sum(xty * beta))/( gene_sf$n_samples - gene_sf$n_covariates - length(idx) ) )

beta_se <- sigma_hat * sqrt(diag( solve( xtx ) ))

# Check hand-calculated estimates & SEs ... 
data.frame(SNP = snps, beta = beta, SE = beta_se)
