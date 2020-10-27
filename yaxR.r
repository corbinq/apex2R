
require(data.table)
require(Matrix)
require(Rcpp)

Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
Rcpp::sourceCpp('xzReader.cpp') 

# xzReader$methods(finalize = function(){})

suffStats <- methods::setRefClass("suffStats", 
	fields = c(
		"Xty" = "numeric",
		"XtX" = "matrix",
		"yty" = "numeric",
		"n_samples" = "numeric",
		"n_covariates" = "numeric",
		"variants" = "character", 
		"trait" = "character"
	)
)

suffStats$methods(
	show = function(){
		cat("\n")
		cat(paste0("  SuffStats for trait '",trait,"'.\n"))
		cat(paste0("    Number of variants = ",nrow(XtX),"\n"))
		cat(paste0("    Number of samples = ",n_samples,"\n"))
		cat(paste0("    Number of technical covariates = ",n_covariates,"\n\n"))
	}
)
suffStats$lock(
"Xty", "XtX", "yty", "n_samples", "n_covariates", "variants", "trait"
)

sumStats <- methods::setRefClass("sumStats",
	fields = c(
		"prefix" = "character",
		"dscore" = "list",
		"dreg" = "data.frame",
		"ds_gene" = "data.table",
		"vb" = "Rcpp_xzReader",
		"vx" = "data.table",
		"kp" = "numeric",
		"genes" = "character",
		"n_snps" = "numeric",
		"n_genes" = "numeric",
		"n_samples" = "numeric",
		"n_covar" = "numeric"
	)
)

sumStats$lock(
	"prefix", "dscore", "dreg", "ds_gene", "vb", "vx", "genes", "n_snps", "n_genes", "n_samples", "n_covar"
)

sumStats$methods(
	initialize = function(file_prefix){
	
		if(missing(file_prefix)) stop("File prefix not specified.")
	
		ss_f <- paste0(file_prefix, ".cis_sumstats.txt.gz")
		sg_f <- paste0(file_prefix, ".cis_gene_table.txt.gz")
		vb_f <- paste0(file_prefix, ".vcov.bin")
		vx_f <- paste0(file_prefix, ".vcov.idx.gz")
		
		for(x in c(ss_f, sg_f, vb_f, vx_f)){
			if(!file.exists(x)) stop(paste0("Required input file '", x, "' does not exist.\n"))
		}

		ds <- lapply(readLines(ss_f), function(x) scan(text = x, quiet=TRUE))
		
		dreg_t <- as.data.frame(do.call(rbind,lapply(ds, function(x) head(x, 4))))
		colnames(dreg_t) <- c("chr", "start", "end", "sv")

		.self$initFields(
			prefix = file_prefix,
			dscore = sapply(ds, function(x) x[-c(1:4)]),
			dreg = dreg_t,
			ds_gene = data.table::fread(sg_f),
			vb = xzReader$new(vb_f),
			vx = fread(vx_f, skip = 1)
		)
		.self$ds_gene[,gene:=gsub("\\..*", "",gene),]
		data.table::setnames(.self$vx, names(.self$vx)[1:10], c("chr", "pos", "ref", "alt", "flipped", "mac", "dV", "rid", "b_s", "b_n"))
		.self$initFields(
			genes = .self$ds_gene$gene,
			n_snps = nrow(vx),
			n_genes =  nrow(.self$dreg),
			n_samples = head(ds_gene$n_samples,1),
			n_covar =  head(ds_gene$n_covar,1)
		)
		.self$kp <- integer(0)
	}
)

sumStats$methods(
	setKeepList = function(keep_list){
		.self$kp <- keep_list
	},
	finalize = function(){},
	show = function(){
		cat("\n")
		cat(paste0("  sumStats querying object.\n\n"))
		cat(paste0("    File prefix = '",.self$prefix,"'\n\n"))
		cat(paste0("    Sample information:\n"))
		cat(paste0("      Sample size = ", .self$n_samples, "\n"))
		cat(paste0("      Number of covariates = ", .self$n_covar, "\n"))
		cat(paste0("    sumStat data:\n"))
		cat(paste0("      Number of genes = ", .self$n_genes, "\n"))
		cat(paste0("      Number of variants = ", .self$n_snps, "\n"))
		cat(paste0("      Total number of tests = ", sum(sapply(.self$dscore, length)), "\n"))
		cat("\n")
	}
)

sumStats$methods(
	getCisRegion = function(gene){
		if( all(is.numeric(gene)) ){
			dreg[gene,]
		}else{
			dreg[match(gene, ds_gene$gene),]
		}
	},
	getSNPs = function(gene, chrom, start, end){
		if(missing(chrom)){
			idx <- getCisRegion(gene)
			chrom <- idx$chr
			start <- idx$start
			end <- idx$end
		}
		inp <- vx[,1:10]
		inp[,rid := seq_along(pos),]
		inp <- subset(inp, chr == chrom & pos >= start & pos <= end)
		inp
	}
)

sumStats$methods(
	getSumStats = function(gene){
		k <- gene
		if( !all(is.numeric(gene)) ){
			k <- match(gene, ds_gene$gene)
		}

		dsnp <- getSNPs(gene)
		
		if( length(kp) == 0 ){
			kp_gene <- TRUE
		}else{
			kp_gene <- which(dsnp$rid %in% kp)
		}
		
		s_gene <- ds_gene[k,]

		U <- dscore[[k]] * s_gene$resid_sd
		dV <- dsnp$dV
		b <- U/dV
		SSR <- (s_gene$n_samples - 1) * (s_gene$resid_sd^2) 
		DFR <- s_gene$n_samples - s_gene$n_covar - 1
		se <- sqrt((SSR/dV - b*b)/DFR)
		zscores <- b/se
		pvals <- pf((zscores^2), 1, s_gene$n_samples - s_gene$n_covar - 1, low = F)

		data.table(dsnp[,1:6], beta = b, se = se, pval = pvals)[kp_gene,]
	},
	getV = function(gene, chrom, start, end, snp_names = TRUE){
		
		if(missing(chrom)){
			idx <- getCisRegion(gene)
			chrom <- idx$chr
			start <- idx$start
			end <- idx$end
		}
		
		ss_idx <- with(vx, which(chr == chrom & pos >= start & pos <= end))
		
		dsnp <- vx[ss_idx,1:10]

		GtU <- as.matrix(vx[ss_idx,-c(1:10)])
		
		x_s <- as.numeric(2 * min(dsnp$b_s))
		x_e <- as.numeric(sum(dsnp$b_n))

		b_start <- as.numeric(dsnp$b_s - dsnp$b_s[1])

		C <- buildMatrixC(
			s = b_start, 
			n = dsnp$b_n, 
			m = dsnp$mac, 
			x = sm$vb$getData(x_s, x_e)
		)

		V <- flipMatrix(
			C - tcrossprod(GtU), 
			which(dsnp$flipped==1), which(dsnp$flipped==0)
		)
		rm(C, GtU)
		
		if( snp_names ){
			colnames(V) <- dsnp[,paste(chr,pos,ref,alt,sep='_'),]
			rownames(V) <- colnames(V)
		}
		if( length(kp) > 0 ){
			kp_g <- which(ss_idx %in% kp)
			V <- V[kp_g, kp_g]
		}
		V
	}
)

sumStats$methods(
	getSuffStats = function(gene){
		
		k <- gene
		if( !all(is.numeric(gene)) ){
			k <- match(gene, ds_gene$gene)
		}
		
		dsnp <- getSNPs(gene)
		s_gene <- ds_gene[k,]
		
		if( length(kp) == 0 ){
			kp_gene <- TRUE
		}else{
			kp_gene <- which(dsnp$rid %in% kp)
		}
		
		U <- dscore[[k]] * s_gene$resid_sd
		dV <- dsnp$dV
		b <- U/dV
		SSR <- (s_gene$n_samples - 1) * (s_gene$resid_sd^2)
		ss <- list()
		tm <- Sys.time()
		ss$XtX <- .self$getV(gene)
		diag(ss$XtX) <- dV[kp_gene]
		tm <- round(Sys.time() - tm, 2)
		cat(paste0("Processed XtX for ", nrow(ss$XtX), " variants in ", tm, " seconds.\n"))

		ss$Xty <- U[kp_gene]
		ss$yty <- SSR
		ss$n_samples <- mean(s_gene$n_samples)
		ss$n_covariates <- mean(s_gene$n_covar)
		
		invisible(ss)
	}
)



metaSumStats <- methods::setRefClass("metaSumStats",
	fields = c(
		"study_data" = "list",
		"prefixes" = "character",
		"studies" = "character",
		"variant_list" = "data.table",
		"genes" = "character"
	)
)

metaSumStats$lock(
	"study_data", "prefixes", "studies", "variant_list"
)

metaSumStats$methods(
	initialize = function(...){
		pfx <- c(...)
		pfx_names <- names(pfx)
		if(is.null(pfx_names)) pfx_names <- pfx
		ni <- length(pfx)
		read_study <- function(ii){
			cat(paste0("Processing study ", ii , " out of ", ni))
			if(!is.null(names(pfx))){
				cat(paste0(": ", names(pfx)[ii], " ... \n"))
			}else{
				cat(" ...\n")
			}
			sumStats$new(pfx[ii])
		}
		.self$initFields(
			prefixes = pfx, 
			studies = pfx_names,
			study_data = lapply(seq_along(pfx), read_study)
		)
		cat("Processed single-study summary statistics, ")
		chr_list <- sapply(.self$study_data, function(x){
			gsub("^chr", "", unique(x$dreg$chr))
		})
		if( length(unique(chr_list)) > 1 ){
			stop(paste("Multiple chromosomes present: ", chr_list))
		}else{
			cat(paste0("chromosome ", unique(chr_list), ".\n"))
		}
		pos_list <- lapply(.self$study_data, function(x){
			x$vx$pos
		})
		ref_list <- lapply(.self$study_data, function(x){
			x$vx$ref
		})
		alt_list <- lapply(.self$study_data, function(x){
			x$vx$alt
		})
		idx_list <- mergeIntersect(pos_list, ref_list, alt_list)
		for(i in 1:length(.self$study_data)){
			.self$study_data[[i]]$setKeepList(idx_list[[i]])
		}
		.self$initFields(
			variant_list = .self$study_data[[1]]$vx[idx_list[[1]],1:4]
		)
	},
	show = function(){
		cat("\n")
		cat(paste0("  metaSumStats querying object.\n\n"))
		cat(paste0("    File prefixes : \n"))
		for(i in 1:length(prefixes)){
		cat(paste0("      ",names(prefixes)[i], ": ",prefixes[i],"\n"))
		}
		cat("\n")
		cat(paste0("    Sample information:\n"))
		cat(paste0("      Total sample size = ", sum(sapply(.self$study_data, function(x) x$n_samples)), "\n"))
		cat(paste0("      Total of covariates = ", sum(sapply(.self$study_data, function(x) x$n_covar)), "\n"))
		cat(paste0("    sumStat data:\n"))
		cat(paste0("      Number of genes = ", length(unique(unlist(lapply(.self$study_data, function(x) x$genes)))), "\n"))
		cat(paste0("      Number of variants = ", nrow(.self$variant_list), "\n"))
		cat("\n")
	},
	finalize = function(){}
)

metaSumStats$methods(
	getSuffStats = function(gene, include_single_studies = FALSE){
		use_studies <- which(sapply(.self$study_data, function(x) gene %in% x$genes))
		sstats <- lapply(.self$study_data[use_studies], function(x){
			x$getSuffStats(gene)
		})
		n_s <- length(use_studies)
		denom <- 0
		for( i in 1:n_s ){
			sstats[[i]]$w <- with(sstats[[i]], 
				(n_samples - n_covariates)/yty
			)
			denom <- denom + sstats[[i]]$w
		}
		out <- list()
		for(x in c("XtX", "Xty", "yty")){
			out[[x]] <- Reduce(`+`, 
				lapply(sstats, function(i){
					i[[x]] * (i$w)
				})
			)
		}
		for(x in c("n_samples", "n_covariates")){
			out[[x]] <- Reduce(`+`, 
				lapply(sstats, function(i){
					i[[x]]
				})
			)
		}
		
		if( include_single_studies ){
			names(sstats) <- .self$studies[use_studies]
			out$studies <- sstats
		}
		out
	}
)

metaSumStats$methods(
	getSumStats = function(gene, intersect = TRUE){
		use_studies <- which(sapply(.self$study_data, function(x) gene %in% x$genes))
		sstats <- data.table::rbindlist(lapply(.self$study_data[use_studies], function(x){
			out <- x$getSumStats(gene)
			out[,`:=`(n = x$ds_gene$n_samples[1], m = x$ds_gene$n_covar[1]),]
			out
		}))
		sstats[,w := 1/(se^2),]
		sstats <- sstats[,list(
			beta = sum(w*beta)/sum(w),
			se = sqrt( sum((w^2)*(se^2))/(sum(w)^2) ),
			n = sum(n), m = sum(m), mac = sum(mac), 
			n_studies = .N
		),by=list(chr,pos,ref,alt)]
		sstats[,pval := pf( (beta/se)^2, 1, n-m-1, low = FALSE),]
		sstats[,n_allele := .N,by=list(chr,pos)]
		if( intersect ){
			sstats <- subset(sstats, n_studies == length(.self$study_data) & n_allele == 1)
		}
		sstats[,list(chr,pos,ref,alt,n_studies,mac,beta,se,pval),][order(chr,pos)]
	}
)

pruneRsq <- function(MAT, thresh = 1){

	kp <- 1:nrow(MAT)
	
	dm <- diag(MAT)
	diag(MAT) <- 0
	
	ctn <- TRUE
	
	nc <- length(kp)
	
	while(ctn){
		for(i in kp){
			wr <- which((MAT[i,kp]^2) >= thresh*(dm[i] * dm[kp]))
			if( length(wr) > 0 ){
				kp <- kp[-wr]
				break;
			}
		}
		if( length(kp) == nc ){
			ctn <- FALSE
		}
		nc <- length(kp)
	}
	kp
}

finemapGene <- function(gene, object, gene_sf = NULL, gene_sm = NULL, L = 10, verbose = TRUE, track_fit = TRUE, ...){

	require(susieR)

	if( is.null(gene_sf) ){
		gene_sf <- object$getSuffStats(gene)
	}
	if( is.null(gene_sm) ){
		gene_sm <- object$getSumStats(gene)
	}

	if( "flipped"  %in% names(gene_sm) ){
		setnames(gene_sm, "flipped", "n_studies")
		gene_sm[,n_studies := as.numeric(NA),]
	}
	gene_sm[,n_studies := as.numeric(n_studies),]

	# diag(gene_sf$XtX) <- abs(diag(gene_sf$XtX))

	print(system.time(
		fm_fit <- susie_suff_stat(
			XtX = gene_sf$XtX, 
			Xty = gene_sf$Xty,
			yty = gene_sf$yty,
			n = gene_sf$n_samples, 
			L = L, verbose = verbose, track_fit = track_fit, ...
		)
	))

	top_snp <- which.min(gene_sm$pval)

	gene_sm[,`:=`(
		"-log10(p-value)" = (-1)*log10(pval),
		snp = paste(chr,pos,ref,alt,sep='_'),
		PIP = fm_fit$pip, signal = apply(fm_fit$alpha, 2, which.max),
		Rsq = (gene_sf$XtX[top_snp,]^2)/(diag(gene_sf$XtX) * gene_sf$XtX[top_snp,top_snp])
	),]
	
	gene_sm$gene <- gene
	
	gene_sm
}
