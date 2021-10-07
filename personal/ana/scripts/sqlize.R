#!/usr/bin/env Rscript

###### description ######
# this script calls the sqlizeR package 
# to create input statements from text files
# usage: ./sqlize.R --help

suppressPackageStartupMessages(library(optparse))

###### options ######
option_list <- list(make_option(c("-f", "--file"),
                                type="character",
                                help="Path to tab delimited input file where missing values are encoded as NA",
                                dest="inputfile",
                                default=NULL),
                    make_option(c("-t", "--table"),
                                type="character",
                                help="SQL table name",
                                dest="tablename",
                                default=NULL),
                    make_option(c("--na"),
                                type="character",
                                help="String to replace NA values\n\t\tDefault: NULL",
                                dest="na_string",
                                default="NULL"),
                    make_option(c("-o", "--output"),
                                type="character",
                                help="Output file prefix",
                                dest="outputfile",
                                default="sqlizeR"))

opt <- parse_args(OptionParser(option_list=option_list))

inputfile <- opt$inputfile
tablename <- opt$tablename
na_string <- opt$na_string
outputfile <- opt$outputfile

if(is.null(tablename)) { stop("Please provide a table name using the --table flag") }

###### main ######

if (!require("sqlizeR")) devtools::install_github("PMBB-Informatics-and-Genomics/sqlizeR")

d <- data.table::fread(inputfile, sep="\t", quote="", fill=TRUE)

if(is.null(na_string)){
        sqlizeR::sqlize(df=d, table=tablename, file=outputfile, save=TRUE)
} else {
        sqlizeR::sqlize(df=d, table=tablename, na_string=na_string, file=outputfile, save=TRUE)
}
