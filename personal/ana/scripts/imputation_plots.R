#!/usr/bin/env Rscript

###### description ######
# this script creates MAF and R2 summary plots 
# for imputed genetic data

suppressPackageStartupMessages(library(optparse))

###### options ######
option_list <- list(make_option(c("-f", "--file"),
		    		type="character",
		    		help="Path to input file with the following columns: SNP, MAF, R2, CHR\n\t\tOrder does not matter but the header should have these column names.",
		    		dest="inputfile",
		   		default=NULL),
		    make_option(c("-o", "--output"),
				type="character",
				help="Output file name",
				dest="outputfile",
				default="plot"),
		    make_option(c("-i", "--image"),
				type="character",
				help="Image type (png, pdf, svg, etc.)\n\t\tdefault: png",
				dest="image",
				default="png"))

opt <- parse_args(OptionParser(option_list=option_list))

inputfile <- opt$inputfile
outputfile <- opt$outputfile
image <- opt$image

###### main ######

suppressPackageStartupMessages(library(ggplot2))

data <- data.table::fread(inputfile, header=TRUE)

# check column names
set <- setdiff(c("CHR", "SNP", "MAF", "R2"), names(data))
if(length(set)!=0){ stop(sprintf("The following columns were not found: %s", paste(set, collapse=", "))) }

data <- as.data.frame(data)

# check data classes
if(class(data$MAF)!="numeric"){
	print("Warning: coercing MAF to numeric")
	data$MAF <- as.numeric(data$MAF)
}

if(class(data$R2)!="numeric"){
	print("Warning: coercing R2 to numeric")
	data$R2 <- data$R2
}

###Rsq by Chromosome
p1 <-  ggplot(data=data, aes(x=factor(CHR, levels=c(1:22, "X")), y=R2)) + 
		geom_boxplot() + 
		xlab("CHR") + 
		ggtitle("Average Rsq per Chromosome") + 
		theme_minimal()
ggsave(p1, filename=paste0(outputfile, "_Rsq_by_CHR", image), dpi=300, height=7, width=14, units="in")

###Distribution of Rsq
p2 <- ggplot(data=data, aes(x=R2)) + 
	geom_density() + 
	theme_minimal() + 
	ggtitle("Distribution of Rsq")
ggsave(p2, filename=paste0(outputfile, "_Rsq_dist", image), dpi=300, height=7, width=7, units="in")

###Rsq by MAF
#p3 <- ggplot(data=data, aes(x=MAF, y=R2)) + geom_hex(bins=100) + theme_minimal() + ggtitle("Rsq by Minor Allele Frequency")
#ggsave(p3, filename=paste0(file, "_Rsq_by_MAF.png"), dpi=300, height=7, width=7, units="in", type="cairo-png")
