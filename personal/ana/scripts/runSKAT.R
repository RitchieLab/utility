#!/usr/bin/env Rscript
#All settings for covariates are ignored so make sure to change them if running with covariates

################ MAIN RUNNING SCRIPT ########################

library(optparse,quietly=TRUE)

opt_list <- list(
		make_option(c("-p", "--prefix"), action="store", default=NULL,
				help="Prefix of input files"),
		make_option(c("-s", "--setid"), action="store", default=NULL,
				help="SKAT setid File"),
		make_option(c("-o", "--out"), action="store", default=NULL,
				help="output file"),
		make_option(c("-b", "--bins"), action="store", default=NULL,
                                help="name of bins file"))

op <- OptionParser(option_list=opt_list)

opts <- parse_args(op)

#load SKAT package
library(SKAT)

#Load necessary files
bed.fn<-paste(opts$prefix,".bed", sep="")
bim.fn<-paste(opts$prefix,".bim", sep="")
fam.fn<-paste(opts$prefix,".fam", sep="")
cov.fn<-paste(opts$prefix,".cov", sep="")
SSD.fn<-paste(opts$out,".SSD",sep="")
INFO.fn<-paste(opts$out,".INFO",sep="")
#SSD.fn<-tempfile()
#INFO.fn<-tempfile()
SKAT.setid<-paste(opts$setid, sep="")

#Generate SSD Files
Generate_SSD_SetID(bed.fn, bim.fn, fam.fn, SKAT.setid, SSD.fn, INFO.fn)

# Read Plink FAM
# flag1=1 = unaffected/accected coded as 0/1
#flag1=0: unaffected/affected coded as 1/2
# Is.binary= TRUE if pheno is binary
FAM.data<- Read_Plink_FAM(fam.fn, Is.binary=TRUE, flag1=0)

#Read Covariate file
# Is.binary= TRUE if pheno is binary
# flag1=1 = unaffected/accected coded as 0/1
FAM.cov<-Read_Plink_FAM_Cov(fam.fn, cov.fn, Is.binary=TRUE, flag1=1,cov_header=TRUE)

#Specify covariates from file
#X1 = FAM.cov[,5]
#X2 = as.factor(FAM.cov[,7])
#X3 = FAM.cov[,8]
X1<-FAM.cov$SequencedGender
X2<-FAM.cov$BIRTH_YEAR
X3<-FAM.cov$PC1
X4<-FAM.cov$PC2 
X5<-FAM.cov$PC2
X6<-FAM.cov$PC4
# Run Null Model w/o covariates
# Change out="C" for continuous
#obj.fn<-SKAT_Null_Model(FAM.data$Phenotype ~1, out="D", Adjustment=FALSE)

# Run Null Model w/ covariates
# Change out_type="C" for continuous, "D" for dichotomous
#obj.fn<-SKAT_Null_Model(FAM.data$Phenotype ~ X1 + X2 + X3 +X4 +X5 + X6, out_type="D", Adjustment=FALSE)
Phenotype <- FAM.data$Phenotype
obj.fn<-SKAT_Null_Model_ChrX(FAM.data$Phenotype ~ X1 + X2 + X3 +X4 +X5 + X6, SexVar="X1", out_type="D", Adjustment=FALSE)
# Open SSD files
SSD.INFO<-Open_SSD(SSD.fn, INFO.fn)

#Run SKAT
#to use Mad Brown Weight use weights.beta=c(0.5,0.5)
#for no weighting, specify weightx.beta=c(1,1) or kernel="linear"
output<-SKAT.SSD.All(SSD.INFO, obj.fn, weights.beta=c(1,1))
#output<-SKAT.SSD.All(SSD.INFO, obj.fn, weights.beta=c(0.5,0.5))

#Create df of results
output.df=output$results

#Write table
write.table(output.df, file= opts$out, col.names=TRUE, row.names=FALSE, quote=FALSE)

#Close SSD files
Close_SSD()
warnings()
# Add Bonferroni Correction to results
#SKAT_results<-read.table(file=paste(opts$out),header=T)
#bonf_adjust<-p.adjust(SKAT_results$P.value, "bonferroni")
#bonf_column<-cbind(c(bonf_adjust))
#colnames(bonf_column) <- "P.value_bonf_adj"
#table<-data.frame(SKAT_results, bonf_column)
#write.table(table, file=paste(opts$out, ".corrected_pval", sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE)

# Output Logistic Regression Summary
# if phe coded 1/2: add -1 to the phenotype
#bins_data<-read.csv(file=opts$bins, head = TRUE)
#df_bins <-data.frame(counts =bins_data[-(1:6),3], phenotype= (bins_data[-(1:6),2]))
#f<- glm(formula= (bins_data[-(1:6),2]) ~ bins_data[-(1:6),3], family= "binomial")
#f<- glm(formula= (bins_data[-(1:6),2]) ~ bins_data[-(1:6),3] + X1 + X2 + X3, family= "binomial")
#capture.output(print(summary(f), prmsd=TRUE, digits=1), file=paste(opts$out, ".regression_summary", sep=""))


