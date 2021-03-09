# read input filenames from command line arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
        stop("usage: script.R <pheno.wsv> <geno.wsv> <pairs.wsv> <prefix>")
}
phenoFile = args[1]
genoFile = args[2]
pairsFile = args[3]
outPrefix = args[4]

# load libraries
cat("loading R packages ...\n")
library(doParallel)
library(tidyverse)
library(caret)  
library(glmnet)
library(dummies)
cat("... done\n")

# read the phenotype data
cat("reading phenotype data from '",phenoFile,"' ...\n")
pheno = read.table(phenoFile, header=TRUE, sep=" ", quote="", colClasses="integer", comment.char="")
cat("... done:",nrow(pheno),"samples,",ncol(pheno),"phenotypes\n")

# read the genotype data header row to identify available variants
cat("reading genotype header from '",genoFile,"' ...\n")
colFilter = read.table(genoFile, header=TRUE, sep=" ", quote="", colClasses="character", nrows=1, comment.char="")
cat("... done:",length(colFilter),"available variants\n")

# read the list of desired variant pairs and identify variants used in any pair
cat("reading variant pairs from '",pairsFile,"' ...\n")
pairs = read.table(pairsFile, header=FALSE, sep=" ", quote="", colClasses="character", comment.char="")
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
colFilter[,variants] = "integer"

# read the full genotype data frame, but including only the needed columns
cat("reading genotype data ...\n")
geno = read.table(genoFile, header=TRUE, sep=" ", quote="", colClasses=as.character(colFilter[1,]), comment.char="")
cat("... done:",nrow(geno),"samples,",ncol(geno),"variants\n")

# check for sample count match
if (nrow(pheno) != nrow(geno)) {
	cat("ERROR: sample count mismatch (",nrow(pheno),"pheno vs",nrow(geno),")\n")
	stop("cannot continue.")
}

# pre-shuffle each phenotype rather than re-doing it for every SNP pair
cat("generating permuted phenotypes ...\n")
phenopermu = list()
for (phenoNum in 1:ncol(pheno)) {
	samplefilter = pheno[,phenoNum] != -9
	phenotypes = pheno[samplefilter,phenoNum]
	cat("... pheno #",phenoNum,"has",length(phenotypes),"non-missing samples ...\n")
	phenopermu[[phenoNum]] = list()
	for (i in 1:1000) {
		set.seed(i) #ask alex about this seed
		rows = sample(length(phenotypes))
		phenopermu[[phenoNum]][[i]] = phenotypes[rows]
	}
}
cat("... done\n")

# define pairwise processing loop function
process <- function(phenoNum, geno1, geno2) {
	set.seed(123)
	
	### MAIN EFFECTS & INTXN TERM - FULL MODEL ###
	
	# combine phenotype column with selected genotype columns and add interaction term (3x3 for genotype data)
	# while dropping sample rows that are missing for the phenotype
	samples = pheno[,phenoNum] != -9
	df = data.frame(Phenotype=pheno[samples,phenoNum], SNP1=geno[samples,geno1], SNP2=geno[samples,geno2], Intxn_Term=paste(geno[samples,geno1], "-", geno[samples,geno2]))
	
	### FIT THE PENALIZED REGRESSION MODEL FOR INTERACTION TERM ONLY ###

	#keep only phenotype and interaction term columns
	#df[,2:3]<- list(NULL)
	
	# Dummy code categorical predictor variables
	x <- dummy.data.frame(df)[,-1]
	x <- as.matrix(x)
	
	# Dummy encode using model.matrix, this reduces 0-0 intxn term
	#x <- model.matrix(Phenotypes~., df)[,-1]
	#[,-1] removes the intercept column
	
	# Convert the outcome (class) to a numerical variable
	y <- df[,1]
	
	# Find the best lambda using cross-validation with technique called "Warm Start" where every lambda in the grid (cv$lambda) is run
	# The default is 10 fold crossvalidation
	cv.ridge <- cv.glmnet(x, y, alpha = 0, family = "binomial") #can't set penalty = 0 
	#plot(cv.ridge)
	cv.ridge$lambda.min 
	#coef(cv.ridge, cv.ridge$lambda.min)
	
	#Try 5 fold cross validation
	#cv.ridge <- cv.glmnet(x, y, alpha = 0, family = "binomial", nfolds = 5)
	#plot(cv.ridge)
	#cv.ridge$lambda.min 
	#coef(cv.ridge, cv.ridge$lambda.min)
	
	#Fit the final model (with the optimal lambda) on all data, and set the penalty to 0 for the main effects
	
	# Dummy version
	ridge.model <- glmnet(x, y, alpha = 0, family = "binomial", lambda = cv.ridge$lambda.min, penalty.factor = 0)
	
	# Model.matrix version
	#ridge.model.x <- glmnet(x, y, alpha = 0, family = "binomial", lambda = cv.ridge$lambda.min, penalty.factor = 0)
	
	# Make matrix that I can append to add the permutation betas
	betas <- as.matrix(ridge.model$beta)
	betas <- data.frame(betas)
	
	# Make matrix that I can append to add the permutation dev.ratios
	dev.ratio <- as.matrix(ridge.model$dev.ratio)
	dev.ratio <- data.frame(dev.ratio)
	
	# Make matrix that I can append to add the permutation nulldev
	nulldev <- as.matrix(ridge.model$nulldev)
	nulldev <- data.frame(nulldev)
	
	### MAIN EFFECTS ONLY - REDUCED MODEL ###
	
	# drop the interaction term from the data frame
	df_reduced <- subset(df, select = -Intxn_Term)
	
	# Dummy code categorical predictor variables
	x_reduced <- dummy.data.frame(df_reduced)[,-1]
	x_reduced <- as.matrix(x_reduced)
	
	# Dummy encode using model.matrix, this reduces 0-0 intxn term
	#x_reduced <- model.matrix(Phenotypes~., df_reduced)[,-1]
	#[,-1] removes the intercept column
	# Convert the outcome (class) to a numerical variable
	y_reduced <- df_reduced[,1]
	
	#run non-penalized regression using glm
	reduced <- glm(y_reduced ~ x_reduced, data = df_reduced, family = binomial) 
	
	#run non-penalized regression using glmnet
	#reduced_glment <- glmnet(x_reduced_dummy, y_reduced, data = df_reduced, family = "binomial") 
	
	#Make matrix that I can append to add the permutation betas
	deviance_reduced <- data.frame(reduced$deviance)
	#deviance in glm = residual deviance which is = 2(LL(Saturated Model) - LL(Proposed Model)) df = df_Sat - df_Proposed
	betas_reduced <- data.frame(reduced$coefficients)
	
	### Run permutation test: Extract Beta, dev.ratio, nulldev for each permutation ###
	
	#Scramble Y, runs glm (lambda is unique for each permuted dataset but same in corresponding full and reduced models)
	
	# Initializing full model vectors
	betas_perm <- data.frame(betas)
	dev.ratio_perm <- data.frame(dev.ratio)
	nulldev_perm <- data.frame(nulldev)
	
	#Initializing reduced model vectors
	betas_reduced_perm <- data.frame(betas_reduced)   
	deviance_reduced_perm <- data.frame(deviance_reduced)
	
	for (i in 1:1000) {
		n <- nrow(df)
		df[,1] = phenopermu[[phenoNum]][[i]]
		# Dummy code categorical predictor variables
		x_perm <- dummy.data.frame(df)[,-1]
		x_perm <- as.matrix(x_perm)
		#x_perm <- model.matrix(Phenotypes~., df)[,-1]
		# Convert the outcome (class) to a numerical variable
		y_perm <- df[,1]
		# Choose lambda separately for each perm
		#set.seed(i) #ask alex about this seed
		cv.ridge_perm <- cv.glmnet(x_perm, y_perm, alpha = 0, family = "binomial")
		cv.ridge_perm$lambda.min 
		coef(cv.ridge_perm, cv.ridge_perm$lambda.min)
		# Run
		#ridge.model_perm <- glmnet(x_perm, y_perm, alpha = 0, family = "binomial", lambda = cv.ridge_perm$lambda.min, penalty.factor = 0)
		ridge.model_perm <- glmnet(x_perm, y_perm, alpha = 0, family = "binomial", lambda = cv.ridge_perm$lambda.min)
		# Save outputs
		betas_perm[i] <- as.matrix(ridge.model_perm$beta)
		betas_perm[i] <- data.frame(betas_perm[i])
		dev.ratio_perm[i] <- as.matrix(ridge.model_perm$dev.ratio)
		dev.ratio_perm[i] <- data.frame(dev.ratio_perm[i])
		nulldev_perm[i] <- as.matrix(ridge.model_perm$nulldev)
		nulldev_perm[i] <- data.frame(nulldev_perm[i])
		
		### Permutation test for reduced model ###
		# Keep same x matrix as above but only keep the main effects
		x_perm_reduced <- x_perm[,1:2]
		# Run
		reduced_perm <- glm(y_perm ~ x_perm_reduced, data = df_reduced, family = "binomial") 
		# Save outputs
		deviance_reduced_perm[i] <- as.matrix(reduced_perm$deviance)
		deviance_reduced_perm[i] <- data.frame(deviance_reduced_perm[i])
		betas_reduced_perm[i] <- as.matrix(reduced_perm$coefficients)
		betas_reduced_perm[i] <- data.frame(betas_reduced_perm[i])
	}
	
	### note ###
	# due to Error: from glmnet Fortran code (error code 10000); All penalty factors are <= 0 I decided to remove the penalty factor
	
	### Make matrix with permuted values ###
	
	# Full model
	betas_perm <- as.matrix(betas_perm)
	betas_perm <- data.frame(betas_perm)
	dev.ratio_perm <- as.matrix(dev.ratio_perm)
	dev.ratio_perm <- data.frame(dev.ratio_perm)
	nulldev_perm <- as.matrix(nulldev_perm)
	nulldev_perm <- data.frame(nulldev_perm)
	
	# Reduced model 
	betas_reduced_perm <- as.matrix(betas_reduced_perm)
	betas_reduced_perm <- data.frame(betas_reduced_perm)
	deviance_reduced_perm <- as.matrix(deviance_reduced_perm)
	deviance_reduced_perm <- data.frame(deviance_reduced_perm)
	
	# Bind the non permuted and permuted outputs and format column names
	
	# Full model
	# note: change headers to reflect permuted and non permuted
	all_betas <- cbind(betas,betas_perm)
	colnames(all_betas)[1] <- 'Non-Permuted'
	all_dev.ratio <- cbind(dev.ratio,dev.ratio_perm)
	colnames(all_dev.ratio)[1] <- 'Non-Permuted'
	colnames(all_dev.ratio)[2] <- 'V1'
	all_nulldev <- cbind(nulldev,nulldev_perm)
	colnames(all_nulldev)[1] <- 'Non-Permuted'
	colnames(all_nulldev)[2] <- 'V1'
	
	# Reduced model
	all_betas_reduced <- cbind(betas_reduced,betas_reduced_perm)
	colnames(all_betas_reduced)[1] <- 'Non-Permuted'
	colnames(all_betas_reduced)[2] <- 'V1'
	all_deviance_reduced <- cbind(deviance_reduced,deviance_reduced_perm)
	colnames(all_deviance_reduced)[1] <- 'Non-Permuted'
	colnames(all_deviance_reduced)[2] <- 'V1'
	
	### Outputs of interest for test ###
	
	# Reduced model calculate deviance
	all_deviance_reduced.t <- as.data.frame(t(as.matrix(all_deviance_reduced)))
	colnames(all_deviance_reduced.t)[1] <- 'deviance_reduced'
	
	# Full model calculate deviance
	# Transpose the deviance ratio and null deviance to make one matrix that we can use to find deviance
	all_dev.ratio.t <- as.data.frame(t(as.matrix(all_dev.ratio)))
	colnames(all_dev.ratio.t)[1] <- 'dev.ratio'
	all_nulldev.t <- as.data.frame(t(as.matrix(all_nulldev)))
	colnames(all_nulldev.t)[1] <- 'nulldev'
	dev <- cbind(all_dev.ratio.t,all_nulldev.t)
	# Calculate our parameter of interest (measure of deviance: dev.ratio = 1-dev/nulldev)
	dev[,3] <- (1-dev[,1])*dev[,2]
	colnames(dev)[3] <- 'deviance_full'
	# dev[,3] contains the deviances for each permuted dataset
	
	### Calculate P-value ###
	
	# Calculate p-value: calculate the number of times the permuted test statistic exceeds the un-permuted and divide by the number of permutations (1000) to get a single p-value corresponding to the interaction term
	
	# Calculate empirical LRT
	dev2 <- cbind(dev$deviance_full, all_deviance_reduced.t)
	eLRT <- as.data.frame(dev2[,2]-dev2[,1])
	
	# Take absolute value of eLRT
	eLRT.abs <- abs(eLRT)
	
	# Identify how many rows are =/> than the first row value (unpermuted)
	count <- as.data.frame(rowSums(eLRT.abs >= eLRT.abs[1,]))
	# Find sum of column and subtract 1 for row 1, divide by total tests
	count_colSum <- (colSums(count) - 1)
	pvalue <- (colSums(count) - 1)/1000
	
	# collect results
	return(cbind(pvalue, eLRT.abs[1,]))
} # process()

# loop over the phenotypes
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
		row = process(phenoNum, genos[1], genos[2])
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
