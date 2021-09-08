#!/usr/bin/env Rscript

# not perfect

# options
library(optparse)
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
pheno = read.table(phenoFile, header=TRUE, sep=" ", quote="", comment.char="", stringsAsFactors=FALSE)
pheno$ID <- as.character(pheno$ID)
cat("... done:",nrow(pheno),"samples,",ncol(pheno),"phenotypes\n")
# check that no phenotypes have -9 as missing value
if(any(pheno == -9, na.rm=TRUE)){
        stop("Value '-9' found in phenotype file. Please ensure missing values are encoded as NA")
}

# read the covariate data
cat("reading covariate data from '",covsFile,"' ...\n")
covar = read.table(covsFile, header=TRUE, sep=" ", quote="", comment.char="", stringsAsFactors=FALSE)
covar$ID <- as.character(covar$ID)
if(any(covar == -9, na.rm=TRUE)){
        warnings("Value '-9' found in covariate file. Treating as a real value.")
}
cat("... done:",nrow(covar),"samples,",ncol(covar),"covariates\n")

# read the genotype data header row to identify available variants
cat("reading genotype header from '",genoFile,"' ...\n")
colFilter = read.table(genoFile, header=TRUE, sep=" ", quote="", nrows=1, comment.char="", stringsAsFactors=FALSE)
cat("... done:",length(colFilter),"available variants\n")

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
cat("... done:",nrow(geno),"samples,",ncol(geno),"variants\n")

# check for sample count match
if (!identical(nrow(covar), nrow(geno), nrow(pheno))) {
        cat("Sample count mismatch (",nrow(pheno),"pheno vs",nrow(geno),"geno vs",nrow(covar),"covar)\n")
        stop("cannot continue.")
}


# make sure data is in same order
merge_df <- pheno[complete.cases(pheno), ] %>%
                inner_join(covar[complete.cases(covar), ], by="ID") %>%
                inner_join(geno[complete.cases(geno), ], by="ID")

covar_cols <- names(covar)[names(covar)!="ID"]
pheno_cols <- names(pheno)[names(pheno)!="ID"]
geno_cols <- names(geno)[names(geno)!="ID"]

pheno <- merge_df[, pheno_cols, drop=FALSE]
covar <- merge_df[, covar_cols, drop=FALSE]
geno <- merge_df[, geno_cols, drop=FALSE]

# pre-shuffle each phenotype rather than re-doing it for every SNP pair
cat("generating permuted phenotypes ...\n")
phenopermu = list()
for (phenoNum in 1:ncol(pheno)) {
        #samplefilter = !is.na(pheno[,phenoNum])
        #covfilter = complete.cases(covar)
        # I think just permute the phenotypes with missingness then apply terminal complete cases as they will always be the same ones missisng?
        fullpheno <- cbind(pheno, covar)
        #phenotypes = pheno[samplefilter,phenoNum]
        #cat("... pheno #",phenoNum,"has",length(phenotypes),"non-missing samples ...\n")
        phenopermu[[phenoNum]] = list() # 1000 permutations for each phenotype
        for (i in 1:1000) {
                set.seed(i) #ask alex about this seed
                complete_pheno <- fullpheno[complete.cases(fullpheno[,c(names(fullpheno)[phenoNum], names(covar))]), c(names(fullpheno)[phenoNum], names(covar))]
                rows = sample(nrow(complete_pheno))
                phenopermu[[phenoNum]][[i]] = complete_pheno[rows, ]
        }
}
cat("... done\n")

# define pairwise processing loop function
process <- function(fullpheno, phenoNum, geno1, geno2) {
        set.seed(123)

        ### MAIN EFFECTS & INTXN TERM - FULL MODEL ###

        # combine phenotype column with selected genotype columns and add interaction term (3x3 for genotype data)
        # while dropping sample rows that are missing for the phenotype
        extract_cols <- c(names(fullpheno)[phenoNum], names(covar))
        samples = complete.cases(fullpheno[,extract_cols])
        complete_pheno <- fullpheno[samples ,extract_cols]
        df = cbind(complete_pheno,
                   data.frame(SNP1=geno[samples,geno1],
                              SNP2=geno[samples,geno2],
                              Intxn_Term=paste(geno[samples,geno1], "-", geno[samples,geno2])))
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
        ridge.model <- glmnet(x, y, alpha = 0, family = "binomial", lambda = cv.ridge$lambda.min, penalty.factor = c(rep(0, ncol(covar)), # covariates
                                                                                                                     rep(0, 2), # main effect
                                                                                                                     rep(1, length(colnames(x)[grep("Intxn_Term", colnames(x))])))) # interaction terms 

        # Model.matrix version
        #ridge.model.x <- glmnet(x, y, alpha = 0, family = "binomial", lambda = cv.ridge$lambda.min, penalty.factor = 0)

        # Make matrix that I can append to add the permutation betas
        betas <- data.frame(as.matrix(ridge.model$beta))
        betas <- data.frame(betas[!grepl(paste0(covar_cols, collapse="|"), row.names(betas)), , drop=FALSE])

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
        y_reduced <- df_reduced[,1] # identical to y — redundant

        #run non-penalized regression using glm
        reduced <- glm(y_reduced ~ x_reduced, data = df_reduced, family = binomial)

        #run non-penalized regression using glmnet
        #reduced_glment <- glmnet(x_reduced_dummy, y_reduced, data = df_reduced, family = "binomial") 

        #Make matrix that I can append to add the permutation betas
        deviance_reduced <- data.frame(reduced$deviance)
        #deviance in glm = residual deviance which is = 2(LL(Saturated Model) - LL(Proposed Model)) df = df_Sat - df_Proposed
        betas_reduced <- data.frame(reduced$coefficients)
        betas_reduced <- betas_reduced[!grepl(paste0(covar_cols, collapse="|"), row.names(betas_reduced)), , drop=FALSE]
        ### Run permutation test: Extract Beta, dev.ratio, nulldev for each permutation ###

        #Scramble Y, runs glm (lambda is unique for each permuted dataset but same in corresponding full and reduced models)

        # Initializing full model vectors
        betas_perm <- data.frame(betas[!row.names(betas) %in% covar_cols, ,drop=FALSE])
        dev.ratio_perm <- data.frame(dev.ratio)
        nulldev_perm <- data.frame(nulldev)

        #Initializing reduced model vectors
        betas_reduced_perm <- data.frame(betas_reduced)
        betas_reduced_perm <- betas_reduced_perm[!grepl(paste0(covar_cols, collapse="|"), row.names(betas_reduced_perm)), , drop=FALSE]
        deviance_reduced_perm <- data.frame(deviance_reduced)

        for (i in 1:1000) {
                # what is this n for?
                #n <- nrow(df)
                df[,extract_cols] = phenopermu[[phenoNum]][[i]]
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
                tmp <- data.frame(as.matrix(ridge.model_perm$beta))
                betas_perm[i] <- tmp[!row.names(tmp) %in% covar_cols, , drop=FALSE]
                dev.ratio_perm[i] <- data.frame(as.matrix(ridge.model_perm$dev.ratio))
                nulldev_perm[i] <- data.frame(as.matrix(ridge.model_perm$nulldev))

                ### Permutation test for reduced model ###
                # Keep same x matrix as above but only keep the main effects
                x_perm_reduced <- x_perm[, c(covar_cols, "SNP1", "SNP2")]
                # Run
                reduced_perm <- glm(y_perm ~ x_perm_reduced, data = df_reduced, family = "binomial")
                # Save outputs
                deviance_reduced_perm[i] <- data.frame(as.matrix(reduced_perm$deviance))
                tmp <- data.frame(as.matrix(reduced_perm$coefficients))
                betas_reduced_perm[i] <- tmp[!grepl(paste0(covar_cols, collapse="|"), row.names(tmp)), , drop=FALSE]
        }
        
          ### note ###
        # due to Error: from glmnet Fortran code (error code 10000); All penalty factors are <= 0 I decided to remove the penalty factor

        ### Make matrix with permuted values ###

        # Full model
        #betas_perm <- as.matrix(betas_perm)
        #betas_perm <- data.frame(betas_perm)
        #dev.ratio_perm <- as.matrix(dev.ratio_perm)
        #dev.ratio_perm <- data.frame(dev.ratio_perm)
        #nulldev_perm <- as.matrix(nulldev_perm)
        #nulldev_perm <- data.frame(nulldev_perm)

        # Reduced model 
        #betas_reduced_perm <- as.matrix(betas_reduced_perm)
        #betas_reduced_perm <- data.frame(betas_reduced_perm)
        #deviance_reduced_perm <- as.matrix(deviance_reduced_perm)
        #deviance_reduced_perm <- data.frame(deviance_reduced_perm)

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
