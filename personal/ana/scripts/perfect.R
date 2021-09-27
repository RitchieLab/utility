#!/usr/bin/env Rscript

# not perfect

###### notes
# The main thing to do in this script is clean up some of the redundancy with the 
# variable assignments. Also some parts could probably stand to be a bit more
# modular.

###### options ######
suppressPackageStartupMessages(library(optparse))
option_list <- list(make_option(c("-p", "--pheno"), 
		                type="character", 
				help="Path to phenotype file",
				dest="phenoFile", 
				default=NULL),
		    make_option(c("-g", "--geno"),
				type="character",
				help="Path to genotype file",
				dest="genoFile",
				default=NULL),
		    make_option(c("-x", "--pairs"),
				type="character",
				dest="pairsFile",
				help="Path to pairs file",
				default=NULL),
		    make_option(c("-c", "--covar"),
				type="character",
				help="Path to covariate file (optional)",
				dest="covsFile",
				default=NULL),
		    make_option(c("-f", "--full_model"),
				type="character",
				help="Full model specification: glmnet or glm\n\t\tDefault: glmnet",
				dest="full_mod",
				default="glmnet"),
		    make_option(c("-r", "--reduced_model"),
				type="character",
				help="Reduced model specification: glmnet or glm\n\t\tDefault: glmnet",
				dest="red_mod",
				default="glmnet"),
		    make_option(c("-l", "--lambda"),
				type="character", 
				help="Lambda value for ridge regression. If lambda is not supplied, minimum cv.ridge lambda will be chosen.\n\t\tDefault: NULL",
				dest="lambda",
				default=NULL),
		    make_option(c("--penalty_main"),
				type="character",
				help="Penalty factors for main effects in a comma separated list\n\t\tDefault: 0,0",
				dest="penalty_main",
				default="0,0"),
		    make_option(c("--penalty_int"),
				type="character",
				help="Penalty factor for interaction terms, will be applied to each term\n\t\tDefault: 1",
				dest="penalty_int",
				default="1"),
		    make_option(c("--penalty_covar"),
				type="character", 
				help="Penalty factors for covariates in a comma separated list\n\t\tDefault: NULL",
				dest="penalty_covs",
				default=NULL),
		    make_option(c("-o", "--output"),
			        type="character",
			        help="Output file prefix",
				dest="outPrefix",
				default="PERFECT.results"))

opt <- parse_args(OptionParser(option_list=option_list))

phenoFile = opt$phenoFile
genoFile = opt$genoFile
pairsFile = opt$pairsFile
covsFile = opt$covsFile
outPrefix = opt$outPrefix
red_mod = opt$red_mod
full_mod = opt$full_mod
penalty_main = as.integer(unlist(strsplit(opt$penalty_main, ",")))
penalty_int = as.integer(opt$penalty_int)

###### end set options ######



###### load libraries ######
cat("loading R packages ...\n")
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(caret))  
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(dummies))
cat("... done\n")

###### end load libraries ######



###### parameter checks ######
# covariate parameter agreement
if(is.null(opt$penalty_covs)){
        if(!is.null(opt$covsFile)){
                stop("Covariate file provided, but no penalty factors provided")
        } else {
                penalty_covs <- NULL
        }
} else {
        if(!is.null(opt$covsFile)){
                stop("Penalty factor set for covariates, but no covariate file provided.")
        } else {
                penalty_covs = as.integer(unlist(strsplit(opt$penalty_covs, ",")))
        }
}

# lambda type
if(is.null(opt$lambda)){ lambda <- NULL } else { lambda <- as.numeric(opt$lambda) }


# function to check equality
all.identical <- function(x) all(mapply(identical, head(x, 1), tail(x, -1)))

# model validity
if("glmnet" %in% c(full_mod, red_mod)  & isTRUE(all.identical(c(penalty_main,0))) & isTRUE(all.identical(c(penalty_int,0)))){
	if(is.null(penalty_covs)){
		stop("'glmnet' was specificied, but penalty factors for all model terms are 0.")	
	} else {
		if(isTRUE(all.identical(c(penalty_covs,0)))){
			stop("'glmnet' was specificied, but penalty factors for all model terms are 0.")
		}
	}
	
}

if(red_mod=="glmnet" & isTRUE(all.identical(c(penalty_main,0)))){
	if(is.null(penalty_covs)){
		stop("'glmnet' was specified as the reduced model, but penalty factors for the model terms are all 0.")		
	} else {
		if(isTRUE(all.identical(c(penalty_covs,0)))){
			stop("'glmnet' was specified as the reduced model, but penalty factors for the model terms are all 0.")
		}
	}
} 

if(red_mod!=full_mod){
	print("Warning: Full and reduced models are not using the same method. Results may be non-comparable if 'lambda' !=0")
}

# warning about glm penalty
if((red_mod=="glmnet" & full_mod=="glm") | (red_mod=="glm" & full_mod=="glmnet")){
	if(isFALSE(all.identical(c(penalty_main,0)))){
		print("Warning: main effect penalties are being ignored in 'glm'")
	} 

	if(!is.null(penalty_covs)){
		if(isFALSE(all.identical(c(penalty_main,0)))){
			print("Warning: covariate penalties are being ignored in 'glm'")
		}
	}
}

if((full_mod=="glm" | red_mod=="glm") & !is.null(lambda)){
	print("Warning: Ignoring user supplied lambda for glm model(s)")
	if(full_mod=="glm" & red_mod=="glm"){
		lambda <- NULL # will avoid unnecessary lambda calculations
	}
}

###### end parameter checks ######



###### read in input files ######
# read the phenotype data
cat("reading phenotype data from '",phenoFile,"' ...\n")
pheno = read.table(phenoFile, header=TRUE, sep=" ", quote="", comment.char="", stringsAsFactors=FALSE)
# check for ID column
if(length(grep("ID", names(pheno)))!=1){
	stop("Phenotype file must contain a single 'ID' column")
}
pheno$ID <- as.character(pheno$ID)
cat("... done:",nrow(pheno),"samples,",ncol(pheno)-1,"phenotypes\n")
# check that no phenotypes have -9 as missing value
if(any(pheno == -9, na.rm=TRUE)){
	stop("Warning: Value '-9' found in phenotype file. Please ensure missing values are encoded as NA") 
}
if(ncol(pheno)>2){
	print("Warning: Only samples that contain no missing values across all phenotypes are considered. Consider providing a single phenotype.")
}

# read the covariate data
if(!is.null(covsFile)){
	cat("reading covariate data from '",covsFile,"' ...\n")
	covar = read.table(covsFile, header=TRUE, sep=" ", quote="", comment.char="", stringsAsFactors=FALSE)
	# check for ID column
	if(length(grep("ID", names(covar)))!=1){
        	stop("Covariate file must contain a single 'ID' column")
	}
	covar$ID <- as.character(covar$ID)
	# check that no phenotypes have -9 as missing value
	if(any(covar == -9, na.rm=TRUE)){
        	print("Warning: Value '-9' found in covariate file. Treating as a real value.")
	}
	cat("... done:",nrow(covar),"samples,",ncol(covar)-1,"covariates\n")
} else {
	print("No covariate file provided")
	covar <- NULL
}

# read the genotype data header row to identify available variants
cat("reading genotype header from '",genoFile,"' ...\n")
colFilter = read.table(genoFile, header=TRUE, sep=" ", quote="", nrows=1, comment.char="", stringsAsFactors=FALSE)
# check for ID column
if(length(grep("ID", names(colFilter)))!=1){
        stop("Genotype file must contain a single 'ID' column")
}
cat("... done:",length(colFilter)-1,"available variants\n")

# read the list of desired variant pairs and identify variants used in any pair
cat("reading variant pairs from '",pairsFile,"' ...\n")
pairs = read.table(pairsFile, header=TRUE, sep=" ", quote="", comment.char="", stringsAsFactors=FALSE)
variants = make.names(unique(union(pairs[,1], pairs[,2])), unique=TRUE)
cat("... done:",nrow(pairs),"pairs of",length(variants),"unique variants\n")

# check for needed variants not found in the header, and prepare the column filter
missing = setdiff(variants, names(colFilter))
if (length(missing) > 0) {
	cat("ERROR:",length(missing),"needed variant(s) not found in data file:\n")
	print(missing)
	stop("cannot continue.")
}
colFilter[,] = "NULL"
colFilter[, "ID"] = "character"
colFilter[,variants] = "integer"

# read the full genotype data frame, but including only the needed columns
cat("reading genotype data ...\n")
geno = read.table(genoFile, header=TRUE, sep=" ", quote="", colClasses=as.character(colFilter[1,]), comment.char="")
cat("... done:",nrow(geno),"samples,",ncol(geno)-1,"variants\n")
# check that no phenotpyes have missnig value
if(any(geno == -9, na.rm=TRUE)){
	stop("Value '-9' found in genotype file. Acceptable genotype values are: 0,1,2,NA.")
}

###### end read input files ######


###### data checks & processing ######
# check for sample count match
myl <- list(geno$ID, covar$ID, pheno$ID)
if(max(lengths(lapply(1:length(myl), function(n) setdiff(myl[[n]], unlist(myl[-n])))))>0){
	#cat("Sample count mismatch (",nrow(pheno),"pheno vs",nrow(geno),"geno vs",nrow(covar),"covar)\n")
	print("Warning: Not all samples found in all input files... taking intersect")
	#stop("cannot continue.")
}

# make sure data is in same order
### WIP this could be improved so that it retains samples with differential phenotype missingness
merge_df <- pheno[complete.cases(pheno), ] %>%
		inner_join(geno[complete.cases(geno), ], by="ID")
if(!is.null(covar)){
	merge_df <- merge_df %>%
			inner_join(covar[complete.cases(covar), ], by="ID")

}

print(paste(nrow(merge_df), "samples remaining after excluding missing data"))

# merge covar if avail
if(!is.null(covar)){
	covar_cols <- names(covar)[names(covar)!="ID"]
} else {
	covar_cols <- "nothing to remove"
}
pheno_cols <- names(pheno)[names(pheno)!="ID"]
geno_cols <- names(geno)[names(geno)!="ID"]

pheno <- merge_df[, pheno_cols, drop=FALSE]
covar <- merge_df[, names(covar)[names(covar)!="ID"], drop=FALSE] # allows for NULL
geno <- merge_df[, geno_cols, drop=FALSE] 

###### end data checks & processing ######


###### permute phenotypes ######
# pre-shuffle each phenotype rather than re-doing it for every SNP pair
cat("generating permuted phenotypes ...\n")
phenopermu = list()
for (phenoNum in 1:ncol(pheno)) {
	fullpheno <- cbind(pheno, covar)
	#cat("... pheno #",phenoNum,"has",length(phenotypes),"non-missing samples ...\n")
	phenopermu[[phenoNum]] = list() # 1000 permutations for each phenotype
	for (i in 1:1000) {
		set.seed(i) 
		complete_pheno <- fullpheno[complete.cases(fullpheno[,c(names(fullpheno)[phenoNum], names(covar))]), c(names(fullpheno)[phenoNum], names(covar)), drop = FALSE]
		rows = sample(nrow(complete_pheno))
		phenopermu[[phenoNum]][[i]] = complete_pheno[rows, ]
	}
}
cat("... done\n")
###### end permute phenotypes ######


###### MAIN FUNCTION ######
# define pairwise processing loop function
process <- function(fullpheno, phenoNum, geno1, geno2) {
	#set.seed(123)
	
	### MAIN EFFECTS & INTXN TERM - FULL MODEL ###
	
	# combine phenotype column with selected genotype columns and add interaction term (3x3 for genotype data)

	extract_cols <- c(names(fullpheno)[phenoNum], names(covar))
	samples = complete.cases(fullpheno[,extract_cols, drop=FALSE])
	complete_pheno <- fullpheno[samples ,extract_cols, drop=FALSE] # this is a bit redundant becase of the complete.cases() in the merge, but I will leave it in as it would be needed to build in support to allow for differential missingness in phenotypes
	
	######################################
	# For testing                        #
	# pairNum <- 1                       #
	# genos = make.names(pairs[pairNum,])#
	# geno1 <- genos[1]		     #
	# geno2 <- genos[2]		     #
	######################################

	df = cbind(complete_pheno,
		   data.frame(SNP1=geno[samples,geno1], 
		              SNP2=geno[samples,geno2], 
		              Intxn_Term=paste(geno[samples,geno1], "-", geno[samples,geno2])))

	### FIT THE PENALIZED REGRESSION MODEL FOR FULL MODEL  ###

	
	# Dummy code categorical predictor variables
	x <- as.matrix(dummy.data.frame(df)[,-1]) # the warning is an update that wasn't fixed by dummies package but shouldn't change anything
	
	# Outcome
	y <- df[,1]
	
	# Find the best lambda using cross-validation with technique called "Warm Start" where every lambda in the grid (cv$lambda) is run
	# The default is 10 fold crossvalidation
	if(is.null(lambda)){
		set.seed(123)
		cv.ridge <- cv.glmnet(x, y, alpha = 0, family = "binomial") #can't set penalty = 0 
		lambda <- cv.ridge$lambda.min
	}
	
	
	#Fit the final model (with the optimal lambda) on all data, and set the penalty to 0 for the main effects
		
	# Full model based on glm or glmnet
	if(full_mod=="glmnet"){
		full.model <- glmnet(x, y, alpha = 0, family = "binomial", lambda = lambda, penalty.factor = c(penalty_covs, # covariates
													       penalty_main, # main effect
													       rep(penalty_int, length(colnames(x)[grep("Intxn_Term", colnames(x))])))) # interaction terms 
	
	
		# Make dataframe to append  permutation output
		betas <- data.frame(as.matrix(full.model$beta))
		dev <- data.frame(full.model$nulldev*(1-full.model$dev.ratio)) # deviance
		#nulldev <- data.frame(ridge.model$nulldev)

	} else {
                full.model <- glm(y~x, family = "binomial")
		if(isTRUE(full.model$converged)){
			betas <- data.frame(as.matrix(full.model$coefficients))
			dev <- data.frame(full.model$null.deviance*(1-full.model$dev.ratio))
			#nulldev <- data.frame(ridge.model$null.deviance)
		} else {
			stop("Full model did not converge")
		}
        } 
	# Remove covariate columns if needed
        betas <- data.frame(betas[!grepl(paste0(c(covar_cols,"Intercept"), collapse="|"), row.names(betas)), , drop=FALSE])

	### MAIN EFFECTS ONLY - REDUCED MODEL ###
	
	# drop the interaction term from the data frame
	df_reduced <- subset(df, select = -Intxn_Term)
	
	# Dummy code categorical predictor variables
	x_reduced <- as.matrix(dummy.data.frame(df_reduced)[,-1])
	
	# Reduced model based on glm or glmnet
	if(red_mod=="glmnet"){
		reduced <- glmnet(x_reduced, y, alpha = 0, family = "binomial", lambda = lambda, penalty.factor = c(penalty_covs, penalty_main)) 
	
		betas_reduced <- data.frame(as.matrix(reduced$beta))
		dev_reduced <- data.frame(reduced$nulldev*(1-reduced$dev.ratio))


	} else {
		reduced <- glm(y ~ x_reduced, family = "binomial")
		if(isTRUE(reduced$converged)){
                	betas_reduced <- data.frame(as.matrix(reduced$coefficients))
			dev_reduced <- data.frame(reduced$deviance)
		} else {
			stop("Reduced model did not converge")
		}
	}
	# Remove covariates if necessary
	betas_reduced <- betas_reduced[!grepl(paste0(c(covar_cols,"Intercept"), collapse="|"), row.names(betas_reduced)), , drop=FALSE]	
	
	### Run permutation test: Extract Beta, dev.ratio, nulldev for each permutation ###
	
	# Scramble Y, runs glm or glmnet 
	# For glm lambda is unique for each permuted dataset but same in corresponding full and reduced models)
	
	for (i in 1:1000) { # append to original dataframe
		y_perm = phenopermu[[phenoNum]][[i]]
	
		# Choose lambda separately for each perm
		set.seed(i)  

		# Don't run if not being used
		if(full_mod=="glmnet" | red_mod=="glmnet"){
			cv.ridge_perm <- cv.glmnet(x, y_perm, alpha = 0, family = "binomial")
		}

		# Run
		if(full_mod=="glmnet"){
			full_model_perm <- glmnet(x, y_perm, alpha = 0, family = "binomial", lambda = cv.ridge_perm$lambda.min, penalty.factor=c(penalty_covs, penalty_main))
			tmp <- data.frame(as.matrix(full_model_perm$beta))
			dev[i+1] <- data.frame(as.matrix(full_model_perm$nulldev*(1-full_model_perm$dev.ratio)))
                	#nulldev[i] <- data.frame(as.matrix(full_model_perm$nulldev))
		} else {
			full_model_perm <- glm(y_perm ~ x, family = "binomial")
			if(isTRUE(full_model_perm$converged)){
                        	tmp <- data.frame(as.matrix(full_model_perm$coefficients))
                        	dev[i+1] <- data.frame(full_model_perm$deviance)
                        	#nulldev[i] <- data.frame(as.matrix(ridge.model_perm$null.deviance))
			} else {
				print(paste("Full model", i, "did not converge"))
				tmp <- NULL; dev[i+1] <- NULL #; nulldev_perm[i] <- NULL
			}
		}
	        # Save outputs
		if(!is.null(tmp)){
			betas[i+1] <- tmp[!grepl(paste0(c(covar_cols,"Intercept"), collapse="|"), row.names(tmp)), , drop=FALSE]
		} else {
			betas[i+1] <- rep(NA, ncol(x))
		}

		### Permutation test for reduced model ###
		# Run
		if(red_mod=="glmnet"){
			reduced_perm <- glmnet(x_reduced, y_perm, alpha = 0, family = "binomial", lambda = cv.ridge_perm$lambda.min, penalty.factor=c(penalty_covs, penalty_main))
			dev_reduced[i+1] <- data.frame(as.matrix(reduced_perm$nulldev*(1-reduced_perm$dev.ratio)))
                        tmp <- data.frame(as.matrix(reduced_perm$beta)) 
		} else {
			reduced_perm <- glm(y_perm ~ x_reduced, family = "binomial")
			if(isTRUE(reduced_perm$converged)){
				dev_reduced[i+1] <- data.frame(as.matrix(reduced_perm$deviance))
                		tmp <- data.frame(as.matrix(reduced_perm$coefficients))	
			} else {
				deviance_reduced_perm[i+1] <- NULL; tmp <- NULL			
			}
		}
		# Save outputs
		if(!is.null(tmp)){
			betas_reduced[i+1] <- tmp[!grepl(paste0(c(covar_cols,"Intercept"), collapse="|"), row.names(tmp)), , drop=FALSE]
		} else {
			betas_reduced[i+1] <- rep(NA, ncol(x_perm_reduced))
		}
	}
	
	### Outputs of interest for test ###
	### WIP I am not sure the rationale behind capturing the betas above as it doesn't look like they are being exported,
	### but will keep in case it is needed in the future	

	### Calculate P-value ###
	
	# Calculate p-value: calculate the number of times the permuted test statistic exceeds the un-permuted and divide by the number of permutations (1000) to get a single p-value corresponding to the interaction term
	
	# Calculate empirical LRT
	eLRT <- as.data.frame(cbind(t(dev), t(dev_reduced)))
	names(eLRT) <- c("full_deviance", "reduced_deviance")
	# deviance(full) = 2(loglikelihood_saturated-loglikelihood_full). deviance(reduced) = 2(loglikelihood_saturated-loglikelihood_reduced)
	# => eLRT = 2(loglikelihood_full-loglikelihood_reduced) = deviance(reduced)-deviance(full)
	eLRT$eLRT <- abs(eLRT$reduced_deviance - eLRT$full_deviance)	
	eLRT_real <- eLRT[1,3]	
	# Identify how many rows are =/> than the unpermuted and divide by total converged permutations
	nconverged <- nrow(eLRT[!is.na(eLRT$eLRT), ])-1
	pvalue <- nrow(eLRT[eLRT$eLRT>=eLRT_real,])/nconverged
	print(paste(nconverged, "permuted models converged and will be used for p-value calculation"))
	# collect results
	return(cbind(pvalue, eLRT_real))
} # process()

###### loop over the phenotypes
numcores = min(16, nrow(pairs))
cat("using",numcores,"cores\n")
registerDoParallel(cores=numcores)
for (phenoNum in 1:ncol(pheno)) {
	phenoName <- colnames(pheno)[phenoNum]
	cat("processing phenotype",phenoNum,"(",phenoName,") ...\n")
	# loop over the desired variant pairs
	results <- foreach (pairNum = 1:nrow(pairs), .combine=rbind, .multicombine=TRUE) %dopar% {
		cat("processing variant pair",pairNum,"(",pairs[pairNum,1],"-",pairs[pairNum,2],") ...\n")
		genos = make.names(pairs[pairNum,])
		row = process(fullpheno, phenoNum, genos[1], genos[2])
		data.frame(pairs[pairNum,1], pairs[pairNum,2], row[1], row[2])
	}
	cat("... all pairs complete.\n")
	
	# print results
	outFile <- paste(outPrefix,"-",phenoName,".txt", sep="")
	cat("writing '",outFile,"' ...\n")
	colnames(results) <- c("SNP1","SNP2","pvalue","eLRT")
	write.table(results, file=outFile, row.names=FALSE, col.names=TRUE)
	cat("... done.\n")
}
stopImplicitCluster()
cat("... all phenotypes complete.\n")

#warnings(): # In model.matrix.default(~x - 1, model.frame(~x - 1), contrasts = FALSE) : non-list contrasts argument ignored

# 4 tests in 1x serial: 4.6 minutes, 330 MB
# 100 tests in 16x parallel: 5.7 minutes, 320 MB
