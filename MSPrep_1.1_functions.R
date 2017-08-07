### Load required packages
### source("http://bioconductor.org/biocLite.R")
### biocLite("pcaMethods")

library(crmn)
library(preprocessCore)
library(sva)
library(pcaMethods)


### Function readdata reads in all the datafiles and summarizes the replicates
readdata_sj <- function (directory, clinicalfile, quantificationfile, linkfile, 
						cvmax = 0.5, missing = 1, linktxt, sean_extra_test_col = F, sean_extra_test_row = F, sean_exclude = NULL, 
						name_list = NULL,
						id_list = NULL,
						id_replace_string = NULL)
{
	# sean_extra_test takes care of a situation in which there are duplicates in the "quant" file, but not replicates.
	# sean_exclude is the name of a metabolite to exclude.
	# obviously, sean_extra_test and sean_exclude functionality were added by Sean Jacobson
	# name_list is where we have a list of replicates - each element in the list is a vector of replicates. If it is NULL, the program 
	# 	tries to figure out what the list is manually.
	# id_list is a list of ids that the program finds in the column labels and creates the name_list using ids (which are presumably found 
	# 	in the column names)
	# id_replace_string is two elements, the string that should be replaced in the IDs, and what it should be replaced by.
	# 	For example, in the Lipids, the ID is listed as 21732E and should be 21731E, so you would have
	# 	id_replace_string = c("21732E", "21731E")
	# 	it can also be a list of sets of 2 if more than one string needs to be fixed
	if(!is.null(name_list) & !is.null(id_list)) stop("only one of name_list and id_list can be specified")
	if(!is.null(id_replace_string) & typeof(id_replace_string) == "character" & length(id_replace_string) != 2) 
		stop("id_replace_string should have a length of 2, the string to be replaced, and the string with which to replace it")
	if(!is.null(id_replace_string) & typeof(id_replace_string) == "list" & any(sapply(id_replace_string, length) != 2))
		stop("elements in id_replace_string should have a length of 2, the strings to be replaced, and the strings with which to replace it")


	quant <- read.table(paste(directory, quantificationfile, 
		sep = ""), header = TRUE, sep = ",")
	if(!is.null(id_replace_string) & typeof(id_replace_string) == "character")
		names(quant) <- gsub(id_replace_string[1], id_replace_string[2], names(quant))
	if(!is.null(id_replace_string) & typeof(id_replace_string) == "list"){
		for(anel in id_replace_string)
			names(quant) <- gsub(anel[1], anel[2], names(quant))
	}
	if(length(sean_exclude) > 0){
		quant <- quant[,-which(names(quant) %in% sean_exclude)]
	}
	if(sean_extra_test_col){
		dupnames <- sapply(names(quant), function(x){unlist(strsplit(x, split = "_"))[1]})[
						duplicated(sapply(names(quant), function(x){unlist(strsplit(x, split = "_"))[1]}))]
		dups <- which(sapply(names(quant), function(x){unlist(strsplit(x, split = "_"))[1]}) %in% dupnames)
		for(adupname in dupnames){
			aa <- quant[,grep(adupname, names(quant))]
			aa$comb <- apply(aa, 1, function(x){if(any(x == 1)) return(max(x)) else return(mean(x))})
			eval(parse(text = sprintf("quant$%s <- aa$comb", gsub("_B_Aq", "_AB_Aq", names(aa)[2]))))
		}
		quant <- quant[, -dups]
	}
	if(sean_extra_test_row){
		#rows
		mzdup_test <- function(){
			mzdup <- quant$mz[duplicated(quant$mz)]
			mzdupnum <- which(quant$mz %in% mzdup)
			quant_mzdup <- quant[mzdupnum, ]
			quant_mzdup <- quant_mzdup[order(quant_mzdup$mz), ]
		}
		mzrt <- apply(quant[,1:2], 1, function(x){paste0(as.character(x), collapse = "_")})
		duprows <- mzrt[duplicated(mzrt)]
		duprs <- which(mzrt %in% duprows)
		for(aduprow in duprows){
			aa <- quant[mzrt == aduprow, ]
			quant <- rbind(quant, apply(aa, 2, mean))
		}
		quant <- quant[-duprs, ]
	}
	# end of initial stuff added by Sean Jacobson
	link <- read.table(paste(directory, linkfile, sep = ""), 
		header = TRUE, sep = ",")
	if(!("COPDGeneID" %in% names(link))) stop("expected to have the string 'COPDGeneID' in the link file")
	
	if(!is.null(id_replace_string) & typeof(id_replace_string) == "character"){
		link$LCMS_Run_Id <- gsub(id_replace_string[1], id_replace_string[2], link$LCMS_Run_Id)
		link$COPDGeneID <- gsub(id_replace_string[1], id_replace_string[2], link$COPDGeneID)
	}
	if(!is.null(id_replace_string) & typeof(id_replace_string) == "list"){
		for(anel in id_replace_string){
			link$LCMS_Run_Id <- gsub(anel[1], anel[2], link$LCMS_Run_Id)
			link$COPDGeneID <- gsub(anel[1], anel[2], link$COPDGeneID)
		}
	}

	clinical <- read.table(paste(directory, clinicalfile, sep = ""), 
		header = TRUE, sep = ",")
	if(length(sean_exclude) > 0){
		excl_clin <- unique(sapply(sean_exclude, function(x){unlist(strsplit(unlist(strsplit(x, split = "X"))[2], split = "_"))[1]}))
		clinical <- clinical[-which(clinical$SubjectID %in% excl_clin), ]
	}
	compounds <- nrow(quant)
	subjects <- nrow(clinical)
	metaf <- matrix(NA, ncol = compounds, nrow = subjects)
	# replicates <- (ncol(quant) - 2)/subjects
	# if (replicates%%1 != 0) {
	# 	 stop("Check datafile for missing replicates")
	# }
	# temp3 <- matrix(NA, ncol = replicates, nrow = subjects * compounds)
	kk <- 1
	rownames(quant) <- paste(quant$mz, "_", quant$rt, sep = "")
	#metab <- subset(quant, select = -c(mz, rt))
	library(dplyr)
	metab <- dplyr::select(quant, -c(mz, rt))
	## More stuff added: here's where we intelligently select replicates
	if(is.null(name_list)){
		if(is.null(id_list))
			unique_ids <- unique(sapply(names(metab), function(x) paste0(unlist(strsplit(unlist(strsplit(x, split = "_"))[1], split = ""))[-1], collapse = "")))
		if(!is.null(id_list))
			unique_ids <- unique(id_list)
		name_list <- lapply(unique_ids, function(x) names(metab)[grep(x, names(metab))])
	}

	temp3 <- matrix(NA, ncol = max(sapply(name_list, length)), nrow = subjects * compounds)
	temp4 <- matrix(NA, ncol = 2, nrow = subjects * compounds)

	for(ii in 1:length(name_list)){#1:subjects) {
		# bcall <- ""
		# for (m in 1:replicates) {
		# 	bcall <- paste(bcall, "metab[,i*replicates-", eval(m - 
		# 		1), "],", sep = "")
		# }
		# bcall2 <- paste("matrix(c(", substr(bcall, 1, nchar(bcall) - 
		# 	1), "),ncol=", replicates, ")", sep = "")
		# temp <- eval(parse(text = bcall2))
		temp <- metab[,name_list[[ii]]]
		replicates <- length(name_list[[ii]])
		temp2 <- matrix(NA, ncol = 4, nrow = compounds)
		# thething <- numeric(compounds)
		for(jj in 1:compounds){
			if(replicates == 1){
				temp2[jj, 1] <- temp[jj]#, 1]
			}else{
				for(m in 1:replicates) {
					if(temp[jj, m] == missing) {
						temp[jj, m] <- NA
					}
				}
				temp2[jj, 1] <- mean(unlist(temp[jj, ]), na.rm = TRUE)
				temp2[jj, 2] <- sd(unlist(temp[jj, ]), na.rm = TRUE)
				temp2[jj, 3] <- temp2[jj, 2]/temp2[jj, 1]
				if (is.na(temp2[jj, 3])) {
					temp2[jj, 3] <- cvmax + 1
				}
				temp2[jj, 4] <- length(na.omit(unlist(temp[jj, ])))
				if (temp2[jj, 1] == "NaN") {
				  temp2[jj, 1] <- 0
				}
				if(temp2[jj, 4] < ceiling(replicates/2) || temp2[jj, 3] > cvmax){
					# thething[jj] <- 1
					temp2[jj, 1] <- 0
					if(temp2[jj, 4] == replicates){
						temp2[jj, 1] <- median(unlist(temp[jj, ]), na.rm = TRUE)
						temp3[kk, ] <- c(unlist(temp[jj, ]), rep(NA, 3 - length(unlist(temp[jj, ]))))
						temp4[kk, 1] <- rownames(metab)[jj]
						temp4[kk, 2] <- name_list[[ii]][1]#colnames(metab)[i]
						kk <- kk + 1
					}
				}
			}
		}
		metaf[ii, ] <- t(temp2[, 1])
	}
	if(max(sapply(name_list, length)) != 1){
		med_comp <- temp3[1:eval((compounds * subjects) - sum(is.na(temp3[, 1]))), ]
		med_comp2 <- temp4[1:eval((compounds * subjects) - sum(is.na(temp3[, 1]))), ]
		med_comp <- cbind(med_comp, med_comp2)
	} else {
		med_comp <- 0
	}
	colnames(metaf) <- rownames(metab)
	link2 <- matrix(colnames(metab), ncol = 1)
	link2 <- matrix(sapply(link2, function(x){gsub("AB_A", "A_A", x)}))
	colnames(link2) <- "LCMS_Run_ID"
	if(all(sapply(link2[,1], function(x){unlist(strsplit(x, split = ""))[1] == "X"})))
		link2[,1] <- sapply(link2[,1], function(x){paste0(unlist(strsplit(x, split = ""))[-1], collapse = "")})
	# added by me:
	if(nrow(link) > 0 & ncol(link) > 0){
		link3 <- unique(merge(link2, link, by.x = "LCMS_Run_ID", 
			by.y = as.character(linktxt))[2])
	}else link3 <- link2

	#link3 <- unique(merge(, link, by.x = "LCMS_Run_ID", 
	#	by.y = as.character(linktxt))[2])
	rownames(metaf) <- link3[, 1]
	list(sum_data1 = metaf, clinical = clinical, medians = med_comp)
}

# this is just a copied and pasted version of the "pca" function from the pcaMethods package.
# I have no idea at all why it works when I do this, but it does, so, I guess I will use it.
pca_sj <- function (object, method, nPcs = 2, scale = c("none", "pareto", 
	"vector", "uv"), center = TRUE, completeObs = TRUE, subset = NULL, 
	cv = c("none", "q2"), ...) 
{
	if (inherits(object, "data.frame")) {
		num <- vapply(object, is.numeric, logical(1))
		if (sum(num) < 2) 
			stop("no numeric data in supplied data.frame")
		Matrix <- as.matrix(object[, num])
	}
	else if (inherits(object, "ExpressionSet")) {
		Matrix <- t(exprs(object))
	}
	else Matrix <- as.matrix(object, rownames.force = TRUE)
	if (!is.null(subset)) 
		Matrix <- Matrix[, subset]
	cv <- match.arg(cv)
	scale <- match.arg(scale)
	if (nPcs > ncol(Matrix)) {
		warning("more components than matrix columns requested")
		nPcs <- min(dim(Matrix))
	}
	if (nPcs > nrow(Matrix)) {
		warning("more components than matrix rows requested")
		nPcs <- min(dim(Matrix))
	}
	if (!checkData(Matrix, verbose = interactive())) 
		stop("Invalid data format.", "Run checkData(data, verbose=TRUE) for details")
	missing <- is.na(Matrix)
	if (missing(method)) {
		if (any(missing)) 
			method <- "nipals"
		else method <- "svd"
	}
	if (any(missing) & method == "svd") {
		warning("data has missing values using nipals instead of user requested svd")
		method <- "nipals"
	}
	method <- match.arg(method, choices = listPcaMethods())
	prepres <- prep(Matrix, scale = scale, center = center, simple = FALSE, 
		...)
	switch(method, svd = {
		res <- svdPca(prepres$data, nPcs = nPcs)
	}, nipals = {
		res <- nipalsPca(prepres$data, nPcs = nPcs)
	}, rnipals = {
		res <- RnipalsPca(prepres$data, nPcs = nPcs)
	}, bpca = {
		res <- bpca(prepres$data, nPcs = nPcs)
	}, ppca = {
		res <- ppca(prepres$data, nPcs = nPcs)
	}, svdImpute = {
		res <- svdImpute(prepres$data, nPcs = nPcs)
	}, robustPca = {
		res <- robustPca(prepres$data, nPcs = nPcs)
	}, nlpca = {
		res <- nlpca(prepres$data, nPcs = nPcs)
	})
	nPcs <- ncol(res@scores)
	if (is.null(scores(res)) | is.null(loadings(res)) | is.null(R2cum(res)) | 
		is.null(method(res))) 
		stop(paste("bad result from pca method", method))
	colnames(res@scores) <- paste("PC", 1:nPcs, sep = "")
	rownames(res@scores) <- rownames(Matrix)
	if (all(dim(loadings(res)) == c(ncol(Matrix), nPcs))) {
		colnames(res@loadings) <- paste("PC", 1:nPcs, sep = "")
		rownames(res@loadings) <- colnames(Matrix)
	}
	if (!is.null(subset)) 
		res@subset <- subset
	res@missing <- missing
	res@nPcs <- nPcs
	res@nObs <- nrow(Matrix)
	res@nVar <- ncol(Matrix)
	res@sDev <- apply(scores(res), 2, sd)
	res@center <- prepres$center
	res@centered <- center
	res@scale <- prepres$scale
	res@scaled <- scale
	res@R2 <- res@R2cum[1]
	if (length(res@R2cum) > 1) 
		res@R2 <- c(res@R2, diff(res@R2cum))
	if (completeObs) {
		cObs <- Matrix
		if (method %in% listPcaMethods("nonlinear")) 
			cObs[missing] <- fitted(res, Matrix, pre = TRUE, 
				post = TRUE)[missing]
		else cObs[missing] <- fitted(res, post = TRUE)[missing]
		res@completeObs <- cObs
	}
	if (cv == "q2") 
		res@cvstat <- Q2(res, Matrix, nruncv = 1, ...)
	return(res)
}



### Filterft to only include compounds that are found in specified percentage of subjects and perform imputation of missing data


filterft_sj <- function (metaf, filterpercent = 0.5) 
{
	count <- matrix(NA, nrow = ncol(metaf), ncol = 1)
	rownames(count) <- colnames(metaf)
	toss <- matrix(1, nrow = eval(ncol(metaf)), ncol = 1)
	for (j in 1:eval(ncol(metaf))) {
		k <- 0
		for (i in 1:dim(metaf)[1]) {
			# sean jacobson added this so that missing data would be ignored
			if(is.na(metaf[i, j])) next
			if (metaf[i, j] == 0) {
				k <- k + 1
			}
		}
		count[j, 1] <- 1 - (k/nrow(metaf))
		if (k >= eval(dim(metaf)[1] * (1 - filterpercent))) {
			toss[j, 1] <- 0
		}
	}
	colnames(toss) <- "toss"
	tests <- cbind(toss, t(metaf))
	metafin <- t(subset(tests, toss != 0)[, -1])
	metaimp <- matrix(NA, nrow = nrow(metafin), ncol = ncol(metafin))
	for (i in 1:nrow(metafin)) {
		for (j in 1:ncol(metafin)) {
			# sean jacobson added this so that missing data would be ignored
			# if(is.na(metafin[i, j])) next
			if (metafin[i, j] == 0) {
				metaimp[i, j] <- NA
			}
			else {
				metaimp[i, j] <- metafin[i, j]
			}
		}
	}
	present <- matrix(1, nrow = nrow(metaf), ncol = ncol(metaf))
	colnames(present) <- colnames(metaf)
	rownames(present) <- rownames(metaf)
	for (j in 1:eval(ncol(metaf))) {
		for (i in 1:dim(metaf)[1]) {
			# sean jacobson added this for missing data
			# if(is.na(metaf[i, j])){present[i,j] <- NA; next}
			if (metaf[i, j] == 0) {
				present[i, j] <- 0
			}
		}
	}
	minval <- metafin
	for (i in 1:nrow(metaf)) {
		for (j in 1:ncol(metafin)) {
			# sean jacobson added this for missing data
			# if(is.na(metafin[i, j])){minval[i,j] <- NA; next}
			if (metafin[i, j] == 0) {
				minval[i, j] <- summary(metaimp[, j])[1]/2
			}
			else {
				minval[i, j] <- metafin[i, j]
			}
		}
	}
	metabpca <- pca_sj(metaimp, nPcs = 3, method = "bpca")
	bpca <- completeObs(metabpca)
	colnames(bpca) <- colnames(metafin)
	rownames(bpca) <- rownames(metafin)
	coldrop <- 1
	for (i in 1:nrow(bpca)) {
		for (j in 1:ncol(bpca)) {
			if (bpca[i, j] < 0) {
				bpca[i, j] <- minval[i, j]
			}
		}
	}
	list(minval = as.data.frame(minval), withzero = as.data.frame(metafin), 
		bpca = as.data.frame(bpca), count = as.data.frame(count))
}

### Create pdf of histogram of distribution of values

graphimputations <- function (directory, minval, withzero, bpca, meanval1 = 0, meanval2 = 1e+06, 
	xmax1 = 5e+05, ymax1 = 25000, xmax2 = 20, ymax2 = 1000, xmax3 = 20, 
	ymax3 = 1500, nbreaks = 500) 
{
	compounds <- ncol(bpca)
	subjects <- nrow(bpca)
	stackdata <- function(dset) {
		temp <- matrix(NA, ncol = 5, nrow = ncol(dset) * nrow(dset))
		comps <- ncol(dset)
		subjects <- nrow(dset)
		for (j in 0:eval(comps - 1)) {
			for (i in 1:subjects) {
				temp[(subjects * j + i), 1] <- as.character(dset[i, 
				  j + 1])
			}
		}
		return(temp)
	}
	keeper <- matrix(0, nrow = compounds, ncol = 1)
	log_zer <- matrix(NA, nrow = subjects, ncol = compounds)
	log_min <- matrix(NA, nrow = subjects, ncol = compounds)
	log_bpca <- matrix(NA, nrow = subjects, ncol = compounds)
	for (i in 1:compounds) {
		if (mean(withzero[, i]) >= meanval1 && mean(withzero[, 
			i]) <= meanval2) {
			keeper[i, 1] <- i
		}
		for (j in 1:subjects) {
			if (withzero[j, i] == 0) {
				log_zer[j, i] <- log(1, 2)
			}
			else {
				log_zer[j, i] <- log(withzero[j, i], 2)
			}
			log_min[j, i] <- log(minval[j, i], 2)
			log_bpca[j, i] <- log(bpca[j, i], 2)
		}
	}
	keep <- subset(keeper, keeper[, 1] > 0)
	minst <- stackdata(minval[, keep])
	zerst <- stackdata(withzero[, keep])
	bpcast <- stackdata(bpca[, keep])
	lminst <- stackdata(log_min[, keep])
	lzerst <- stackdata(log_zer[, keep])
	lbpcast <- stackdata(log_bpca[, keep])
	pdf(paste(directory, "ImpMethods_CompsBetween_", meanval1, 
		"-", as.character(meanval2), "_bk", nbreaks, "_", Sys.Date(), 
		".pdf", sep = ""))
	par(mfrow = c(2, 1))
	ylim1 <- c(0, ymax1)
	xlim1 <- c(0, xmax1)
	ylim2 <- c(0, ymax2)
	xlim2 <- c(0, xmax2)
	ylim3 <- c(0, ymax3)
	xlim3 <- c(0, xmax3)
	x <- (as.numeric(as.character(zerst[, 1])))
	hist(x, breaks = nbreaks, col = "red", xlab = "With Zeros", 
		main = "Histogram of Abundances with Zeros", ylim = ylim1, 
		xlim = xlim1)
	x <- (as.numeric(as.character(lzerst[, 1])))
	hist(x, breaks = nbreaks, col = "red", xlab = "With Zeros", 
		main = "Histogram of Log2 Abundances with Zeros", ylim = ylim2, 
		xlim = xlim2)
	x <- (as.numeric(as.character(minst[, 1])))
	hist(x, breaks = nbreaks, col = "red", xlab = "With Min Val", 
		main = "Histogram of Abundances with 1/2 Minimum Value Replacement", 
		ylim = ylim1, xlim = xlim1)
	x <- (as.numeric(as.character(lminst[, 1])))
	hist(x, breaks = nbreaks, col = "red", xlab = "With Min Val", 
		main = "Histogram of Log2 Abundances with 1/2 Minimum Value Replacement", 
		ylim = ylim3, xlim = xlim2)
	x <- (as.numeric(as.character(bpcast[, 1])))
	hist(x, breaks = nbreaks, col = "red", xlab = "BPCA", main = "Histogram of Abundances with BPCA Imputation", 
		ylim = ylim1, xlim = xlim1)
	x <- (as.numeric(as.character(lbpcast[, 1])))
	hist(x, breaks = nbreaks, col = "red", xlab = "BPCA", main = "Histogram of Log2 Abundances with BPCA Imputation", 
		ylim = ylim3, xlim = xlim2)
	dev.off()
}


# diagnosticgen

diagnosticgen <- function (testobj, clindat, link1, batch, pheno, directory, ylim1, 
	ylim2, ylim3) 
{
	pca.point.pch.batch = function(model.id) {
		if (model.id == "null") {
			return("black")
		}
		else if (model.id == 1) {
			return("1")
		}
		else if (model.id == 2) {
			return("2")
		}
		else if (model.id == 3) {
			return("3")
		}
		else if (model.id == 4) {
			return("4")
		}
		else if (model.id == 5) {
			return("5")
		}
		else if (model.id == 6) {
			return("6")
		}
		else if (model.id == 7) {
			return("7")
		}
		else if (model.id == 8) {
			return("8")
		}
		else if (model.id == 9) {
			return("9")
		}
		else if (model.id == 10) {
			return("A")
		}
		else if (model.id == 11) {
			return("B")
		}
		else if (model.id == 13) {
			return("C")
		}
		else if (model.id == 14) {
			return("D")
		}
		else if (model.id == 15) {
			return("E")
		}
		else {
			return("99")
		}
	}
	pca.point.color.change = function(model.id) {
		if (is.na(model.id)) {
			return("black")
		}
		else if (model.id == 0) {
			return("green")
		}
		else if (model.id == 1) {
			return("blue")
		}
		else if (model.id == 2) {
			return("yellow")
		}
		else if (model.id == 3) {
			return("red")
		}
		else if (model.id == 4) {
			return("purple")
		}
		else {
			return("black")
		}
	}
	pcaplots <- function(dset, tlab, link1, batch, pheno) {
		dset$mid <- rownames(dset)
		test <- merge(dset, clindat, by.x = "mid", by.y = as.character(link1))
		pca1 <- prcomp(test[, 2:ncol(dset)])
		title = c(paste("PCA for ", tlab, sep = ""))
		title2 = c("Color Indicates Phenotype, Number Indicates Batch")
		xlab1 = paste("PCA 1, % Var = ", summary(pca1)$importance[2, 
			1], sep = "")
		ylab1 = paste("PCA 2, % Var = ", summary(pca1)$importance[2, 
			2], sep = "")
		temp2 <- subset(test, select = c(batch, pheno))
		colnames(temp2) <- c("batch", "pheno")
		plot(pca1$x[, 1], pca1$x[, 2], main = title, sub = title2, 
			xlab = xlab1, ylab = ylab1, pch = sapply(temp2$batch, 
				pca.point.pch.batch), col = sapply(temp2$pheno, 
				pca.point.color.change))
		return()
	}
	stackdata <- function(dset, link1, pheno, batch) {
		temp <- matrix(NA, ncol = 5, nrow = ncol(dset) * nrow(dset))
		comps <- ncol(dset)
		subjects <- nrow(dset)
		dset$mid <- rownames(dset)
		dset <- merge(dset, clindat, by.x = "mid", by.y = link1)
		for (j in 0:eval(comps - 1)) {
			for (i in 1:subjects) {
				temp[(subjects * j + i), 1] <- as.character(dset[i, 
				  j + 2])
				temp[(subjects * j + i), 2] <- as.character(dset$mid[i])
				temp[(subjects * j + i), 3] <- as.character(subset(dset, 
				  select = pheno)[i, 1])
				temp[(subjects * j + i), 4] <- as.character(subset(dset, 
				  select = batch)[i, 1])
			}
		}
		return(temp)
	}
	data_r <- stackdata(testobj$log_data, link1, pheno, batch)
	data_n <- stackdata(testobj$log_quant, link1, pheno, batch)
	data_rc <- stackdata(as.data.frame(testobj$log_data_combat), 
		link1, pheno, batch)
	data_nc <- stackdata(as.data.frame(testobj$log_quant_combat), 
		link1, pheno, batch)
	data_med <- stackdata(testobj$med_adj, link1, pheno, batch)
	data_medc <- stackdata(as.data.frame(testobj$med_combat), 
		link1, pheno, batch)
	data_ruv <- stackdata(as.data.frame(testobj$ruv_adj), link1, 
		pheno, batch)
	data_sva <- stackdata(as.data.frame(testobj$sva_adj), link1, 
		pheno, batch)
	data_crmn <- stackdata(testobj$crmn_adj, link1, pheno, batch)
	pdf(paste(directory, "Diagnostics", Sys.Date(), ".pdf", sep = ""))
	par(mfrow = c(2, 2))
	xlab1 <- link1
	xlab2 <- pheno
	xlab3 <- batch
	pcaplots(testobj$log_data, "Non-Normalized", link1, batch, 
		pheno)
	temp2 <- data_r
	boxplot(as.numeric(temp2[, 1]) ~ temp2[, 2], main = "", xlab = xlab1, 
		ylab = "Abundance", ylim = ylim2)
	boxplot(as.numeric(temp2[, 1]) ~ temp2[, 3], main = "", xlab = xlab2, 
		ylab = "Abundance", ylim = ylim2)
	boxplot(as.numeric(temp2[, 1]) ~ temp2[, 4], main = "", xlab = xlab3, 
		ylab = "Abundance", ylim = ylim2)
	pcaplots(as.data.frame(testobj$log_data_combat), "Raw with Combat", 
		link1, batch, pheno)
	temp2 <- data_rc
	boxplot(as.numeric(temp2[, 1]) ~ temp2[, 2], main = "", xlab = xlab1, 
		ylab = "Abundance", ylim = ylim2)
	boxplot(as.numeric(temp2[, 1]) ~ temp2[, 3], main = "", xlab = xlab2, 
		ylab = "Abundance", ylim = ylim2)
	boxplot(as.numeric(temp2[, 1]) ~ temp2[, 4], main = "", xlab = xlab3, 
		ylab = "Abundance", ylim = ylim2)
	pcaplots(testobj$log_quant, "Quantile Normalization", link1, 
		batch, pheno)
	temp2 <- data_n
	boxplot(as.numeric(temp2[, 1]) ~ temp2[, 2], main = "", xlab = xlab1, 
		ylab = "Abundance", ylim = ylim2)
	boxplot(as.numeric(temp2[, 1]) ~ temp2[, 3], main = "", xlab = xlab2, 
		ylab = "Abundance", ylim = ylim2)
	boxplot(as.numeric(temp2[, 1]) ~ temp2[, 4], main = "", xlab = xlab3, 
		ylab = "Abundance", ylim = ylim2)
	pcaplots(as.data.frame(testobj$log_quant_combat), "Quantile with Combat", 
		link1, batch, pheno)
	temp2 <- data_nc
	boxplot(as.numeric(temp2[, 1]) ~ temp2[, 2], main = "", xlab = xlab1, 
		ylab = "Abundance", ylim = ylim2)
	boxplot(as.numeric(temp2[, 1]) ~ temp2[, 3], main = "", xlab = xlab2, 
		ylab = "Abundance", ylim = ylim2)
	boxplot(as.numeric(temp2[, 1]) ~ temp2[, 4], main = "", xlab = xlab3, 
		ylab = "Abundance", ylim = ylim2)
	pcaplots(testobj$med_adj, "Median", link1, batch, pheno)
	temp2 <- data_med
	boxplot(as.numeric(temp2[, 1]) ~ temp2[, 2], main = "", xlab = xlab1, 
		ylab = "Abundance", ylim = ylim1)
	boxplot(as.numeric(temp2[, 1]) ~ temp2[, 3], main = "", xlab = xlab2, 
		ylab = "Abundance", ylim = ylim1)
	boxplot(as.numeric(temp2[, 1]) ~ temp2[, 4], main = "", xlab = xlab3, 
		ylab = "Abundance", ylim = ylim1)
	pcaplots(as.data.frame(testobj$med_combat), "Median with Combat", 
		link1, batch, pheno)
	temp2 <- data_medc
	boxplot(as.numeric(temp2[, 1]) ~ temp2[, 2], main = "", xlab = xlab1, 
		ylab = "Abundance", ylim = ylim1)
	boxplot(as.numeric(temp2[, 1]) ~ temp2[, 3], main = "", xlab = xlab2, 
		ylab = "Abundance", ylim = ylim1)
	boxplot(as.numeric(temp2[, 1]) ~ temp2[, 4], main = "", xlab = xlab3, 
		ylab = "Abundance", ylim = ylim1)
	pcaplots(as.data.frame(testobj$sva_adj), "SVA Adjusted", 
		link1, batch, pheno)
	temp2 <- data_sva
	boxplot(as.numeric(temp2[, 1]) ~ temp2[, 2], main = "", xlab = xlab1, 
		ylab = "Abundance", ylim = ylim2)
	boxplot(as.numeric(temp2[, 1]) ~ temp2[, 3], main = "", xlab = xlab2, 
		ylab = "Abundance", ylim = ylim2)
	boxplot(as.numeric(temp2[, 1]) ~ temp2[, 4], main = "", xlab = xlab3, 
		ylab = "Abundance", ylim = ylim2)
	pcaplots(as.data.frame(testobj$ruv_adj), "RUV Adjusted", 
		link1, batch, pheno)
	temp2 <- data_ruv
	boxplot(as.numeric(temp2[, 1]) ~ temp2[, 2], main = "", xlab = xlab1, 
		ylab = "Abundance", ylim = ylim2)
	boxplot(as.numeric(temp2[, 1]) ~ temp2[, 3], main = "", xlab = xlab2, 
		ylab = "Abundance", ylim = ylim2)
	boxplot(as.numeric(temp2[, 1]) ~ temp2[, 4], main = "", xlab = xlab3, 
		ylab = "Abundance", ylim = ylim2)
	pcaplots(testobj$crmn_adj, "CRMN", link1, batch, pheno)
	temp2 <- data_crmn
	boxplot(as.numeric(temp2[, 1]) ~ temp2[, 2], main = "", xlab = xlab1, 
		ylab = "Abundance", ylim = ylim3)
	boxplot(as.numeric(temp2[, 1]) ~ temp2[, 3], main = "", xlab = xlab2, 
		ylab = "Abundance", ylim = ylim3)
	boxplot(as.numeric(temp2[, 1]) ~ temp2[, 4], main = "", xlab = xlab3, 
		ylab = "Abundance", ylim = ylim3)
	dev.off()
}

# Katerina's function, unedited

normdata.new_kk <- function (metafin, clindat, link1, pheno, batch, ncont = 10, 
	controls = c(), ncomp = 2) 
{
	compounds <- ncol(metafin)
	subjects <- nrow(metafin)
	metan <- t(normalize.quantiles(t(metafin)))
	colnames(metan) <- colnames(metafin)
	rownames(metan) <- rownames(metafin)
	final <- as.data.frame(metafin)
	finaln <- metan
	convert <- function(dset) {
		s_dset <- as.data.frame(matrix(NA, ncol = ncol(dset), 
			nrow = nrow(dset)))
		for (j in 1:ncol(dset)) {
			s_dset[, j] <- log2(as.numeric(dset[, j]))
		}
		colnames(s_dset) <- colnames(dset)
		rownames(s_dset) <- rownames(dset)
		return(s_dset)
	}
	log_final <- convert(final)
	log_finaln <- convert(finaln)
	combat <- function(dset, pheno, batch) {
		compounds <- ncol(dset)
		temp <- dset
		temp$mid <- rownames(temp)
		test <- merge(temp, clindat, by.x = "mid", by.y = as.character(link1))
		d1 <- t(dset)
		colnames(d1) <- dset[, 1]
		d2 <- subset(test, select = pheno)
		d2a <- model.matrix(~as.factor(d2[, 1]))
		d3 <- subset(test, select = batch)
		count <- matrix(0, nrow = compounds, ncol = 1)
		d1mod <- matrix(0, nrow = compounds, ncol = subjects)
		for (j in 1:compounds) {
			d1mod[j, ] <- log2(as.numeric(d1[j, ]))
		}
		comadj <- t(ComBat(d1mod, mod = d2a, batch = d3[,1])) # Here the element "d2a" is Gold Stage phenotype and "d3" is the batch.
		colnames(comadj) <- colnames(dset)
		rownames(comadj) <- rownames(dset)
		return(comadj)
	}
	final_rc <- combat(as.data.frame(final), pheno, batch)
	final_com <- combat(as.data.frame(finaln), pheno, batch)
	svaestimate <- function(dset, pheno, batch) {
		d1 <- t(dset[, 1:compounds])
		colnames(d1) <- rownames(dset)
		rownames(d1) <- colnames(dset)
		dset$mid <- rownames(dset)
		test <- merge(dset, clindat, by.x = "mid", by.y = as.character(link1))
		d2 <- subset(test, select = c(batch))
		d3 <- subset(test, select = c(pheno))
		colnames(d3)[1] <- "pheno"
		d3m <- model.matrix(~as.factor(pheno), data = d3)
		count <- matrix(0, nrow = compounds, ncol = 1)
		d1mod <- d1
		svaadj <- sva(d1mod, d3m, method = "irw")
		return(as.matrix(svaadj$sv, ncol = svaadj$n.sv))
	}
	sva <- as.data.frame(svaestimate(log_final, pheno, batch))
	for (i in 1:ncol(sva)) {
		colnames(sva)[i] <- paste("f", i, sep = "")
	}
	ruv <- function(dset, k, ctl, pheno) {
		d1 <- t(dset[, 1:compounds])
		count <- matrix(0, nrow = compounds, ncol = 1)
		d1mod <- d1
		Y <- d1mod
		dset$mid <- rownames(dset)
		test <- merge(dset, clindat, by.x = "mid", by.y = link1)
		X <- subset(test, select = pheno)
		Z = matrix(rep(1, ncol(Y)))
		RZY = Y - Y %*% Z %*% solve(t(Z) %*% Z) %*% t(Z)
		W = svd(RZY[ctl, ])$v
		W = W[, 1:k]
		return(as.matrix(W, ncol = k))
	}
	dset <- final
	control <- matrix(NA, nrow = compounds, ncol = 6)
	for (j in 1:compounds) {
		control[j, 1] <- colnames(dset)[j]
		control[j, 2] <- j
		control[j, 3] <- mean(dset[, j], na.rm = TRUE)
		control[j, 4] <- sd(dset[, j], na.rm = TRUE)
		control[j, 5] <- sd(dset[, j], na.rm = TRUE)/mean(dset[, 
			j], na.rm = TRUE)
		control[j, 6] <- sum(dset[, j] == 0, na.rm = TRUE)
	}
	control2 <- control[order(control[, 6], control[, 5]), ]
	ctl <- as.numeric(control2[, 2][1:max(ncont, ncol(sva) + 
		5)])
	if (length(controls) > 0) {
		ctl <- controls
	}
	contout <- cbind(ctl, colnames(log_final)[ctl])
	ruv_raw <- as.data.frame(ruv(log_final, ncol(sva), ctl, pheno))
	for (i in 1:ncol(ruv_raw)) {
		colnames(ruv_raw)[i] <- paste("f", i, sep = "")
	}
	genadj <- function(dset, factors) {
		out <- sapply(1:compounds, function(j) lm(dset[, j] ~ 
			as.matrix(factors, ncol = 1))$fitted.values)
		colnames(out) <- colnames(dset)
		rownames(out) <- rownames(dset)
		return(out)
	}
	ruv_adj <- genadj(log_final, ruv_raw)
	sva_adj <- genadj(log_final, sva)
	final_shift <- final
	Y <- t(final_shift)
	Ya <- t(final_shift)[ctl, ]
	Ys <- t(final_shift)[-ctl, ]
	j <- 1
	isIS <- rep(FALSE, compounds)
	ctlo <- ctl[order(ctl)]
	for (i in 1:compounds) {
		if (j <= 10) {
			if (ctlo[j] == i) {
				isIS[i] <- TRUE
				j <- j + 1
			}
			else {
				isIS[i] <- FALSE
			}
		}
	}
	dset <- final_shift
	dset$mid <- rownames(dset)
	test <- merge(dset, clindat, by.x = "mid", by.y = link1)
	X <- subset(test, select = c(pheno))
	colnames(X)[1] <- "pheno"
	G <- model.matrix(~-1 + as.factor(X$pheno)) # Here, "G" comes from the Gold Stage phenotype.
	normed.crmn <- normalize(Y, "crmn", factors = G, standards = isIS, 
		ncomp = ncomp)
	lnormed.crmn <- convert(normed.crmn)
	normed.med <- normalize(Y, "median", factors = G, standards = isIS)
	lnormed.med <- convert(normed.med)
	final_crmn <- as.data.frame(t(lnormed.crmn))
	final_med <- as.data.frame(t(lnormed.med))
	med_com <- combat(as.data.frame(t(normed.med)), pheno, batch)
	list(log_data = log_final, log_data_combat = final_rc, log_quant = log_finaln, 
		log_quant_combat = final_com, med_adj = final_med, med_combat = med_com, 
		sva_factors = sva, sva_adj = sva_adj, ruv_factors = ruv_raw, 
		ruv_adj = ruv_adj, crmn_adj = final_crmn, controls = contout)
}


## Katerina's function, edited by me for multiple phenotypes and only ComBat

normdata.new_multpheno <- function (metafin, clindat, link1, pheno, batch, ncont = 10, 
	controls = c(), ncomp = 2, combat_log = F, do_crmn = F) 
{
	# if(length(pheno) == 1) stop("if pheno is only 1 element, use normdata.new_kk")
	if(length(pheno) > 1){
		remove_subjects <- which(apply(clindat[,pheno], 1, function(x) any(is.na(x))))
	}else if(length(pheno) == 1){
		remove_subjects <- which(is.na(clindat[,pheno]))
	}else remove_subjects <- numeric(0)
	if(length(remove_subjects) > 1){
		metafin <- metafin[-remove_subjects, ]
		clindat <- clindat[-remove_subjects, ]
	}
	compounds <- ncol(metafin)
	subjects <- nrow(metafin)
	metan <- t(normalize.quantiles(t(metafin)))
	colnames(metan) <- colnames(metafin)
	rownames(metan) <- rownames(metafin)
	final <- as.data.frame(metafin)
	finaln <- metan
	convert <- function(dset) {
		s_dset <- as.data.frame(matrix(NA, ncol = ncol(dset), 
			nrow = nrow(dset)))
		for (j in 1:ncol(dset)) {
			s_dset[, j] <- log2(as.numeric(dset[, j]))
		}
		colnames(s_dset) <- colnames(dset)
		rownames(s_dset) <- rownames(dset)
		return(s_dset)
	}
	# log_final <- convert(final)
	# log_finaln <- convert(finaln)
	combat <- function(dset, pheno, batch) {
		compounds <- ncol(dset)
		temp <- dset
		temp$mid <- rownames(temp)
		test <- merge(temp, clindat, by.x = "mid", by.y = as.character(link1))
		d1 <- t(dset)
		colnames(d1) <- dset[, 1]
		d2 <- subset(test, select = pheno)
		# d2_p <- apply(d2, 1, paste0, collapse = "")
		pheno_oc <- sapply(pheno, function(x) length(unique(d2[,x])) > 7)
		if(length(pheno) > 0){
			d2a <- eval(parse(text = sprintf("model.matrix(~%s,data=d2)", paste0(sapply(pheno, function(x, oc){
				if(x %in% pheno[oc]){
					res <- x
				}else res <- paste0("as.factor(", x, ")")
				return(res)
			}, pheno_oc), collapse = "+") )))
		}else d2a <- NULL
		d3 <- subset(test, select = batch)
		count <- matrix(0, nrow = compounds, ncol = 1)
		d1mod <- matrix(0, nrow = compounds, ncol = nrow(d3))#subjects)

		for (j in 1:compounds) {
			ifelse(combat_log, d1mod[j, ] <- log2(as.numeric(d1[j, ])), d1mod[j, ] <- as.numeric(d1[j, ]))
		}
		comadj <- t(ComBat(d1mod, mod = d2a, batch = d3[,1]))
		colnames(comadj) <- colnames(dset)
		rownames(comadj) <- rownames(dset)
		return(comadj)
	}
	final_rc <- combat(as.data.frame(final), pheno, batch)
	final_com <- combat(as.data.frame(finaln), pheno, batch)
	# svaestimate <- function(dset, pheno, batch) {
	# 	d1 <- t(dset[, 1:compounds])
	# 	colnames(d1) <- rownames(dset)
	# 	rownames(d1) <- colnames(dset)
	# 	dset$mid <- rownames(dset)
	# 	test <- merge(dset, clindat, by.x = "mid", by.y = as.character(link1))
	# 	d2 <- subset(test, select = c(batch))
	# 	d3 <- subset(test, select = c(pheno))
	# 	colnames(d3)[1] <- "pheno"
	# 	d3m <- model.matrix(~as.factor(pheno), data = d3)
	# 	count <- matrix(0, nrow = compounds, ncol = 1)
	# 	d1mod <- d1
	# 	svaadj <- sva(d1mod, d3m, method = "irw")
	# 	return(as.matrix(svaadj$sv, ncol = svaadj$n.sv))
	# }
	# sva <- as.data.frame(svaestimate(log_final, pheno, batch))
	# for (i in 1:ncol(sva)) {
	# 	colnames(sva)[i] <- paste("f", i, sep = "")
	# }
	# ruv <- function(dset, k, ctl, pheno) {
	# 	d1 <- t(dset[, 1:compounds])
	# 	count <- matrix(0, nrow = compounds, ncol = 1)
	# 	d1mod <- d1
	# 	Y <- d1mod
	# 	dset$mid <- rownames(dset)
	# 	test <- merge(dset, clindat, by.x = "mid", by.y = link1)
	# 	X <- subset(test, select = pheno)
	# 	Z = matrix(rep(1, ncol(Y)))
	# 	RZY = Y - Y %*% Z %*% solve(t(Z) %*% Z) %*% t(Z)
	# 	W = svd(RZY[ctl, ])$v
	# 	W = W[, 1:k]
	# 	return(as.matrix(W, ncol = k))
	# }
	if(do_crmn){
		dset <- final
		control <- matrix(NA, nrow = compounds, ncol = 6)
		for (j in 1:compounds) {
			control[j, 1] <- colnames(dset)[j]
			control[j, 2] <- j
			control[j, 3] <- mean(dset[, j], na.rm = TRUE)
			control[j, 4] <- sd(dset[, j], na.rm = TRUE)
			control[j, 5] <- sd(dset[, j], na.rm = TRUE)/mean(dset[, 
				j], na.rm = TRUE)
			control[j, 6] <- sum(dset[, j] == 0, na.rm = TRUE)
		}
		control2 <- control[order(control[, 6], control[, 5]), ]
		ctl <- as.numeric(control2[, 2][1:max(ncont, ncol(sva) + 
			5)])
		if (length(controls) > 0) {
			ctl <- controls
		}
		# contout <- cbind(ctl, colnames(log_final)[ctl])
		# ruv_raw <- as.data.frame(ruv(log_final, ncol(sva), ctl, pheno))
		# for (i in 1:ncol(ruv_raw)) {
		# 	colnames(ruv_raw)[i] <- paste("f", i, sep = "")
		# }
		# genadj <- function(dset, factors) {
		# 	out <- sapply(1:compounds, function(j) lm(dset[, j] ~ 
		# 		as.matrix(factors, ncol = 1))$fitted.values)
		# 	colnames(out) <- colnames(dset)
		# 	rownames(out) <- rownames(dset)
		# 	return(out)
		# }
		# ruv_adj <- genadj(log_final, ruv_raw)
		# sva_adj <- genadj(log_final, sva)
		final_shift <- final
		Y <- t(final_shift)
		# Ya <- t(final_shift)[ctl, ]
		# Ys <- t(final_shift)[-ctl, ]
		j <- 1
		isIS <- rep(FALSE, compounds)
		ctlo <- ctl[order(ctl)]
		for (i in 1:compounds) {
			if (j <= 10) {
				if (ctlo[j] == i) {
					isIS[i] <- TRUE
					j <- j + 1
				}
				else {
					isIS[i] <- FALSE
				}
			}
		}
		dset <- final_shift
		dset$mid <- rownames(dset)
		test <- merge(dset, clindat, by.x = "mid", by.y = link1)
		X <- subset(test, select = c(pheno))
		pheno_oc <- sapply(pheno, function(x) length(unique(X[,x])) > 7)
		if(length(pheno) > 0){
			G <- eval(parse(text = sprintf("model.matrix(~-1 + %s,data=X)", paste0(sapply(pheno, function(x, oc){
				if(x %in% pheno[oc]){
					res <- x
				}else res <- paste0("as.factor(", x, ")")
				return(res)
			}, pheno_oc), collapse = "+") )))
		}else G <- model.matrix(~1, data = X)
			
		# # colnames(X)[1] <- "pheno"
		# # d2_p <- apply(d2, 1, paste0, collapse = "")
		# # d2a <- model.matrix(~as.factor(d2_p))
		# G <- model.matrix(~-1 + as.factor(apply(X, 1, paste0, collapse = "")))
		normed.crmn <- normalize(Y, "crmn", factors = G, standards = isIS, 
			ncomp = ncomp)
		lnormed.crmn <- convert(normed.crmn)
		normed.med <- normalize(Y, "median", factors = G, standards = isIS)
		lnormed.med <- convert(normed.med)
		final_crmn <- as.data.frame(t(lnormed.crmn))
		final_med <- as.data.frame(t(lnormed.med))
	}else final_crmn <- final_med <- NULL
	# med_com <- combat(as.data.frame(t(normed.med)), pheno, batch)
	list(data = final, data_combat = final_rc, quant = finaln, 
		quant_combat = final_com, crmn_adj = final_crmn, med_adj = final_med)
		# data_combat is logged if using quant normalization
		#  med_combat = med_com, 
		# sva_factors = sva, sva_adj = sva_adj, ruv_factors = ruv_raw, 
		# ruv_adj = ruv_adj, controls = contout)
}	


## Katerina's function, edited to not use any phenotypes

normdata.new_nopheno <- function (metafin, clindat, link1, batch, ncont = 10, 
	controls = c(), ncomp = 2) 
{
	compounds <- ncol(metafin)
	subjects <- nrow(metafin)
	metan <- t(normalize.quantiles(t(metafin)))
	colnames(metan) <- colnames(metafin)
	rownames(metan) <- rownames(metafin)
	final <- as.data.frame(metafin)
	finaln <- metan
	convert <- function(dset) {
		s_dset <- as.data.frame(matrix(NA, ncol = ncol(dset), 
			nrow = nrow(dset)))
		for (j in 1:ncol(dset)) {
			s_dset[, j] <- log2(as.numeric(dset[, j]))
		}
		colnames(s_dset) <- colnames(dset)
		rownames(s_dset) <- rownames(dset)
		return(s_dset)
	}
	# log_final <- convert(final)
	# log_finaln <- convert(finaln)
	combat <- function(dset, pheno = NULL, mybatch) {
		compounds <- ncol(dset)
		temp <- dset
		temp$mid <- rownames(temp)
		testqq <- merge(temp, clindat, by.x = "mid", by.y = as.character(link1))
		d1 <- t(dset)
		colnames(d1) <- dset[, 1]
		d2 <- NULL# subset(testqq, select = pheno)
		d2a <- NULL# model.matrix(~as.factor(d2[, 1]))
		d3 <- subset(testqq, select = mybatch)
		# print(dim(d3))
		count <- matrix(0, nrow = compounds, ncol = 1)
		d1mod <- matrix(0, nrow = compounds, ncol = subjects)
		for (j in 1:compounds) {
			d1mod[j, ] <- as.numeric(d1[j, ]) #log2(as.numeric(d1[j, ]))
		}
		comadj <- t(ComBat(d1mod, mod = NULL, batch = d3[,1])) # Here the element "d2a" is Gold Stage phenotype and "d3" is the batch.
		colnames(comadj) <- colnames(dset)
		rownames(comadj) <- rownames(dset)
		return(comadj)
	}
	final_rc <- combat(as.data.frame(final), mybatch = batch)
	final_com <- combat(as.data.frame(finaln), pheno, mybatch = batch)
	# svaestimate <- function(dset, pheno = NULL, mybatch) {
	# 	d1 <- t(dset[, 1:compounds])
	# 	colnames(d1) <- rownames(dset)
	# 	rownames(d1) <- colnames(dset)
	# 	dset$mid <- rownames(dset)
	# 	test <- merge(dset, clindat, by.x = "mid", by.y = as.character(link1))
	# 	d2 <- subset(test, select = c(batch))
	# 	d3 <- subset(test, select = c(pheno))
	# 	colnames(d3)[1] <- "pheno"
	# 	d3m <- model.matrix(~as.factor(pheno), data = d3)
	# 	count <- matrix(0, nrow = compounds, ncol = 1)
	# 	d1mod <- d1
	# 	svaadj <- sva(d1mod, d3m, method = "irw")
	# 	return(as.matrix(svaadj$sv, ncol = svaadj$n.sv))
	# }
	# sva <- as.data.frame(svaestimate(log_final, pheno, batch))
	# for (i in 1:ncol(sva)) {
	# 	colnames(sva)[i] <- paste("f", i, sep = "")
	# }
	# ruv <- function(dset, k, ctl, pheno) {
	# 	d1 <- t(dset[, 1:compounds])
	# 	count <- matrix(0, nrow = compounds, ncol = 1)
	# 	d1mod <- d1
	# 	Y <- d1mod
	# 	dset$mid <- rownames(dset)
	# 	test <- merge(dset, clindat, by.x = "mid", by.y = link1)
	# 	X <- subset(test, select = pheno)
	# 	Z = matrix(rep(1, ncol(Y)))
	# 	RZY = Y - Y %*% Z %*% solve(t(Z) %*% Z) %*% t(Z)
	# 	W = svd(RZY[ctl, ])$v
	# 	W = W[, 1:k]
	# 	return(as.matrix(W, ncol = k))
	# }
	dset <- final
	control <- matrix(NA, nrow = compounds, ncol = 6)
	for (j in 1:compounds) {
		control[j, 1] <- colnames(dset)[j]
		control[j, 2] <- j
		control[j, 3] <- mean(dset[, j], na.rm = TRUE)
		control[j, 4] <- sd(dset[, j], na.rm = TRUE)
		control[j, 5] <- sd(dset[, j], na.rm = TRUE)/mean(dset[, 
			j], na.rm = TRUE)
		control[j, 6] <- sum(dset[, j] == 0, na.rm = TRUE)
	}
	control2 <- control[order(control[, 6], control[, 5]), ]
	ctl <- as.numeric(control2[, 2][1:max(ncont, ncol(sva) + 
		5)])
	# if (length(controls) > 0) {
	# 	ctl <- controls
	# }
	# contout <- cbind(ctl, colnames(log_final)[ctl])
	# ruv_raw <- as.data.frame(ruv(log_final, ncol(sva), ctl, pheno))
	# for (i in 1:ncol(ruv_raw)) {
	# 	colnames(ruv_raw)[i] <- paste("f", i, sep = "")
	# }
	# genadj <- function(dset, factors) {
	# 	out <- sapply(1:compounds, function(j) lm(dset[, j] ~ 
	# 		as.matrix(factors, ncol = 1))$fitted.values)
	# 	colnames(out) <- colnames(dset)
	# 	rownames(out) <- rownames(dset)
	# 	return(out)
	# }
	# ruv_adj <- genadj(log_final, ruv_raw)
	# sva_adj <- genadj(log_final, sva)
	final_shift <- final
	Y <- t(final_shift)
	# Ya <- t(final_shift)[ctl, ]
	# Ys <- t(final_shift)[-ctl, ]
	j <- 1
	isIS <- rep(FALSE, compounds)
	ctlo <- ctl[order(ctl)]
	for (i in 1:compounds) {
		if (j <= 10) {
			if (ctlo[j] == i) {
				isIS[i] <- TRUE
				j <- j + 1
			}
			else {
				isIS[i] <- FALSE
			}
		}
	}
	# dset <- final_shift
	# dset$mid <- rownames(dset)
	# test <- merge(dset, clindat, by.x = "mid", by.y = link1)
	# X <- subset(test, select = c(pheno))
	# colnames(X)[1] <- "pheno"
	G <- model.matrix(~-1 + rep(1, ncol(Y))) # Here, "G" comes from the Gold Stage phenotype.
	# G <- model.matrix(~-1 + as.factor(apply(X, 1, paste0, collapse = "")))
	normed.crmn <- normalize(Y, "crmn", factors = G, standards = isIS,
		ncomp = ncomp)
	# lnormed.crmn <- convert(normed.crmn)
	# normed.med <- normalize(Y, "median", factors = G, standards = isIS)
	# lnormed.med <- convert(normed.med)
	final_crmn <- as.data.frame(t(lnormed.crmn))
	# final_med <- as.data.frame(t(lnormed.med))
	# med_com <- combat(as.data.frame(t(normed.med)), pheno, batch)
	list(data = final, data_combat = final_rc, quant = finaln, 
		quant_combat = final_com)
}
