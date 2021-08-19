
source("Apex2R.r")

# single study

prefix <- "chr20"

ss <- sumStats(prefix)

# meta-analysis

system.time(ssm <- metaSumStats(
	studyA ="chr20.A",
	studyB ="chr20.B",
	studyC = "chr20.C"
))

# genes present in all studies	
i_genes <- Reduce(intersect, lapply(ssm$study_data, function(x) x$genes))

gene <- "ENSG00000101004"

## Get XtY and XtX objects for specified gene
gene_sf <- ssm$getSuffStats(gene)
## Get single-variant summary statistics for specified gene
gene_sm <- ssm$getSumStats(gene)

require(ggplot2)

# Find the top SNP for this gene
top_snp <- which.min(gene_sm$pval)

# Construct nice SNP labels, and calculate Rsq with top SNP
gene_sm[,`:=`(
	snp = paste(chr,pos,ref,alt,sep='_'),
	Rsq = (gene_sf$XtX[top_snp,]^2)/(diag(gene_sf$XtX) * gene_sf$XtX[top_snp,top_snp])
),]

# Make a little manhattan plot
ggplot(gene_sm, aes(x = pos/(1000000), y= log10(1/pval), colour = discr_rsq(Rsq))) + geom_point(size=0.75) + theme_minimal() + coord_cartesian(expand=FALSE,clip='off') + ggsci::scale_colour_locuszoom() + ylab(NULL) + xlab("Chromosome position") + theme(strip.placement = "outside", panel.spacing = unit(2, "lines")) + guides(colour = guide_legend(title = expression(italic(r)^2~"with lead SNP"), override.aes = list(shape = 15, size = 4))) + ggtitle(levels(factor(tmp_l$gene))[ii]) + theme(panel.spacing = unit(0.65, "lines"), axis.text=element_text(size=7, color = "black"), title = element_text(size = 9))

# function to convert correlation matrix to long table 
getCorr <- function(x){
	require(data.table)
	data.table
	out <- melt(data.table(as.data.frame(cov2cor(x)), keep.rownames = "snp1"), id.var = "snp1", variable.name = "snp2", value.name = "corr")
	out$snp1 <- factor(out$snp1, levels = rownames(x), order = TRUE)
	out$snp2 <- factor(out$snp2, levels = rownames(x), order = TRUE)
	out
}

# Nearest 200 SNPs surrounding the top SNP 
wk <- ( which.min(ss$pval) - 100 ):( which.min(ss$pval) + 100 )

# Plot the LD for SNPs in the top-hit region
p1 <- ggplot(getCorr(test[[1]]$XtX[wk,wk]), aes(x = snp1, y = snp2, fill = corr^2)) + geom_tile() +theme_void() +  scale_fill_distiller(palette = "RdYlBu", limits = c(0,1)) + theme(legend.position = "none") +  theme(plot.margin = unit(c(1,0.5,1,0.5)/5, "cm")) + ggtitle(ssm$studies[1])


