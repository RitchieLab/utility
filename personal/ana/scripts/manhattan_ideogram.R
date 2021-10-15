library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(stringr)

#Change filename
#Change icd9_cat -> icd9_top_cat
#Change Clinical_Lab -> ICD9_Category

df <- read.csv("megaphewas_1e-04_icd9_cat.tsv", sep="\t", header=TRUE)

chrpos <- read.csv("hg19.chrom.sizes", sep="\t", header=TRUE)
ideo <- read.csv("cytoBandIdeo.txt", sep="\t", header=FALSE)
names(ideo) <- c("chr", "xmin", "xmax", "pos", "color")

plot_list_new <- list()

df$gwas_trait2 <- ifelse(df$gwas_trait=="Y", "Yes", "No")
#Change to overall_pval_lrt for ICD9
df$overall_pval_log <- -log10(df$overall_pval_db)
#For clinical lab file
#Comment out for ICD9
df$pos <- as.numeric(as.character(df$pos))

#Change to 15 for ICD9
colourCount <- 25
getPalette = colorRampPalette(brewer.pal(11, "Spectral"))

for(i in 1:22){

	findchr <- paste("chr", i, sep="")

	ideo1 <- ideo[ideo$chr==findchr, ]
	ideo1.2 <- ideo1[ideo1$color!="gneg", ]
	acen1 <- ideo1[ideo1$color=="acen",]
	dat1 <- df[df$chr==i,]

	xmaxchrpos <- chrpos$length[chrpos$chr==i]

	#Calculate offset so chromosomes align at the bottom
	offset <- chrpos[1,2] - xmaxchrpos

	#Position width of chromosome
	ymaxchr=floor(min(dat1$overall_pval_log))
	#xlimit <- max(max(dat1$overall_pval_log), 11)
	chrhgt <- (max(max(dat1$overall_pval_log), 11)-ymaxchr)*0.0625
	#chrght <- 0.7/xlimit
	yminchr=ymaxchr-chrhgt

	#Offset data
	acen1$xmin <- acen1$xmin + offset
	acen1$xmax <- acen1$xmax + offset
	ideo1.2$xmin <- ideo1.2$xmin + offset
	ideo1.2$xmax <- ideo1.2$xmax + offset
	dat1$pos <- dat1$pos + offset
	
	label_pos_x <- max(dat1$overall_pval_log)
	dat1$result_annotation <- paste( unlist(str_split(dat1$result_annotation," ",n=2))[1],  unlist(str_split(dat1$result_annotation," ",n=2))[2], sep="\n")
	dat1$result_annotation[dat1$overall_pval_log!=label_pos_x] <- NA
	dat1$result_annotation <- as.character(dat1$result_annotation)
	#label_text <- as.character(dat1$result_annotation[which(dat1$ovoerall_pval_log==label_pos_x)]	

	p <- ggplot() + geom_point(data=dat1, aes(y=pos, x=overall_pval_log, colour=icd9_cat, shape=gwas_trait2, alpha=1/5), size=3) + scale_colour_manual(values = getPalette(colourCount)) + theme_bw() + theme(legend.position="none", panel.background = element_rect(fill = 'white', colour = 'white'), panel.grid.major = element_line(colour = "gray34"), panel.grid.minor=element_line(colour="gray34", size=.1), panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(), panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank(), plot.margin=unit(c(1,1,1,-1), "lines"), axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_text(size=15))
	p <- p + coord_cartesian(xlim=c(yminchr, max(11, max(dat1$overall_pval_log))))
	p <- p + geom_rect(data=ideo1.2, aes(ymin=xmin, ymax=xmax, xmin=yminchr, xmax=ymaxchr), fill="gray80", alpha=0.75)
	p <- p + geom_rect(data=acen1, aes(ymin=xmin, ymax=xmax, xmin=yminchr, xmax=ymaxchr), fill="red", alpha=0.3)

	p <- p + geom_rect(mapping=aes(ymax=xmaxchrpos+offset, ymin=0+offset, xmax=ymaxchr, xmin=yminchr), colour="black", alpha=0.01) + ylim(0, chrpos[1,2]) + theme(axis.title.x=element_blank())
	p <- p + theme(panel.border=element_blank())

	p <- p + geom_segment(aes(x=11, xend=11, y=0+offset, yend=xmaxchrpos+offset), colour="red", linetype=2)

	#p <- p + geom_text(data=dat1, aes(y=pos, x=overall_pval_log, label=result_annotation), vjust=-0.5, hjust="inward", size=2.7)


	p <- p + scale_y_reverse(lim=c(chrpos[1,2], 0))

	#If you want to save the plot
	#assign(paste("n", i, sep=""), p)

	plot_list_new[[i]] <- grid.arrange(p, ncol=1)

	#Save each chromosome individually
	#ggsave(p, filename=paste("n", i, ".png", sep=""), dpi=300, units="in")

}

#Save plot with all chromosomes
#Each chromosome is 5x3
#png("MegaPheWAS_manhattan_ideogram_no_labels_p1E-04_icd9_cat_annotated.png", width=33, height=10, units="in", res=300, pointsize=3)

#do.call("grid.arrange", c(plot_list_new, ncol=11))

#dev.off()

png("MegaPheWAS_manhattan_ideogram_no_labels_p1E-04_icd9_cat_annotated.png", width=15, height=25, units="in", res=300, pointsize=3)

do.call("grid.arrange", c(plot_list_new, ncol=5))
dev.off()

#Make legend that you have to crop yourself
#dat1$GWAS_Catalog <- ifelse(dat1$gwas_trait=="Y", "Yes", "No")
df$GWAS_Catalog <- df$gwas_trait2
#Change to column 7 for ICD9
colnames(df)[4] <- "Clinical_Lab"
p <- ggplot() + geom_point(data=df, aes(y=pos, x=overall_pval_log, colour=Clinical_Lab, shape=GWAS_Catalog)) + scale_colour_manual(values=getPalette(colourCount)) +  theme_bw()
ggsave(p, filename="icd9_cat_legend.png", height=6, width=18, units="in", dpi=300)

