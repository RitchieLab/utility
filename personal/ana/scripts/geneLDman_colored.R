#' gldman
#'
#' Create Manhattan plots for GWAS
#' @param d data frame, must contain SNP, CHR (Gene), POS, pvalue columns, optional Shape and Color
#' @param line optional pvalue threshold to draw red line at
#' @param log10 plot -log10() of pvalue column, boolean
#' @param yaxis label for y-axis, automatically set if log10=TRUE
#' @param opacity opacity of points, from 0 to 1, useful for dense plots
#' @param annotate_snp vector of SNPs to annotate
#' @param annotate_p pvalue threshold to annotate
#' @param title optional string for plot title
#' @param chrcolor1 first alternating color for chromosome
#' @param chrcolor2 second alternating color for chromosome
#' @param highlight_snp vector of snps to highlight
#' @param highlight_p pvalue threshold to highlight
#' @param highlighter color to highlight
#' @param groupcolors named vector of colors for data in 'Color' column
#' @param background variegated or white
#' @param chrblocks boolean, turns on x-axis chromosome marker blocks
#' @param file file name of saved image
#' @param hgt height of plot in inches
#' @param wi width of plot in inches
#' @param res resolution of plot in pixels per inch
#' @import ggplot2
#' @return png image
#' @export
#' @family GWAS functions
#' @family static plotting functions
#' @seealso \code{\link{igman}}, \code{\link{agman}}, \code{\link{pheman}}, \code{\link{eman}}
#' @examples
#' data(gwas)
#' gldman(d=gwas[gwas$Frame=="Additive", 1:4], line=0.0005, title="GWAS Example: Additive")

  library(ggplot2)
  # Assuming df in format        
  # SNP CHR BP P gwassnp Gene R2
  df <- as.data.frame(data.table::fread("ld_snps_pos_p_densearea_r2_wgene.txt"))
  names(df) <- c("SNP", "chrom", "POS", "pvalue", "gwassnp", "CHR", "R2")
  levs <- unique(df[order(df$chrom, df$POS), 6])
  d <- df[, c(1,6,3,4,7)]

  #Sort data
  #d$Gene <- factor(d$Gene, levels = levs)
  #d2 <- d[order(d$Gene),]
  #d2$CHR <- factor(d2$CHR, levels=unique(d2$CHR))
  d$CHR <- factor(d$CHR, levels = levs)
  d <- d[, colnames(d)!="Gene"]
  d_order <- d[order(d$CHR, d$POS), ]
  #d_order$pos_index <- seq.int(nrow(d_order))
  d_order$pos_index <- ave(d_order$POS, d_order$CHR, FUN=seq_along)
  d_order_sub <- d_order[colnames(d_order) %in% c("SNP", "CHR", "POS", "pvalue", "pos_index")]

  #Set up dataframe with color and position info
  maxRows <- by(d_order_sub, d_order_sub$CHR, function(x) x[which.max(x$pos_index),])
  minRows <- by(d_order_sub, d_order_sub$CHR, function(x) x[which.min(x$pos_index),])
  milimits <- do.call(rbind, minRows)
  malimits <- do.call(rbind, maxRows)
  lims <- merge(milimits, malimits, by="CHR")
  names(lims) <- c("Color", "snpx", "px", "posx", "posmin", "snpy", "py", "posy", "posmax")
  lims$av <- (lims$posmin + lims$posmax)/2
  lims <- lims[order(lims$Color),]
  lims$shademap <- rep(c("shade_ffffff", "shade_ebebeb"), length.out=nrow(lims), each=1)
  
  lims$nposmin <- seq_along(lims$Color) - 1
  lims$nposmax <- seq_along(lims$Color)
  lims$navg <- (lims$nposmin + lims$nposmax)/2
  
  #d_order2 <- merge(d_order, lims[, c(1,5,9,12)], by.x="CHR", by.y="Color")
  #d_order2$xpos <- d_order2$nposmin + (d_order2$pos_index/d_order2$posmax)
  d_order2 <- merge(d_order, lims[, c(1,5,9,12)], by.x="CHR", by.y="Color")
  d_order2$xpos <- ifelse( d_order2$posmax-d_order2$posmin <= 0, d_order2$nposmin, 
                           d_order2$nposmin + ((d_order2$pos_index-d_order2$posmin)/(d_order2$posmax-d_order2$posmin)))
  
  d_order <-  d_order2[, c(2,1,3:5,10)]
  colnames(d_order)[6] <- "pos_index"
  
  names(lims) <- c("Gene", "snpx", "px", "posx", "oposmin", "snpy", "py", "posy", 
                   "oposmax", "oav", "shademap", "posmin", "posmax", 'av')
  
  
  yaxismin <- min(d_order$pval)
  d_order$pval <- -log10(d_order$pvalue)
  yaxislab <- expression(paste("-log"[10], "(p-value)", sep=""))
  redline <- -log10(5e-8)
  cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  p <- ggplot() + geom_rect(data = lims, aes(xmin = posmin, xmax = posmax, ymin = yaxismin, ymax = Inf, fill=factor(shademap)), alpha = 0.5)
  p <- p + geom_point(data=d_order, aes(x=pos_index, y=pval, color=R2), alpha=0.5,  stroke=1, shape=16)
  p <- p + scale_fill_manual(name = "Color", values = c("#FFFFFF", "#EBEBEB"), guides(alpha=FALSE))
  p <- p + scale_x_continuous(breaks=lims$av, labels=lims$Gene, expand=c(0,0))
  p <- p + theme(axis.text.x=element_text(angle=90), panel.grid.minor.x = element_blank(), panel.grid.major.x=element_blank(), axis.title.x=element_blank(), legend.position="bottom", legend.title=element_blank())
  p <- p + ggtitle(title) + ylab(yaxislab)
  p <- p + geom_hline(yintercept = redline, colour="red")
  p <- p + scale_y_continuous(limits=c(yaxismin, max(max(d_order$pval), -log10(5e-8))), expand=expand_scale(mult=c(0,0.1)))
  p <- p + scale_color_gradientn(name='R2', colours=cbbPalette[c(5,2,4,6)])
  p <- p + ggtitle("Dense Area")
  p <- p + guides(fill="none")
  ggsave(p, filename = "ld_snps_pos_p_densearea_r2_wgene.png", height=7, width=12, dpi=600)
  
  
