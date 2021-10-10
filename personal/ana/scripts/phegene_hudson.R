#' phegene
#'
#' Create Manhattan plots for PheWAS
#' @param d data frame, must contain PHE, SNP, CHR, POS, pvalue columns, optional Shape
#' @param phegroup optional grouping file for phenotypes, must contain PHE and Group columns
#' @param line optional vector of y values to draw red line at
#' @param log10 plot -log10() of pvalue column, boolean
#' @param yaxis label for y-axis, automatically set if log10=TRUE
#' @param ymax set the upper limit for the y-axis if not automatically scaled
#' @param opacity opacity of points, from 0-1, useful for dense plots
#' @param annotate_snp vector of SNPs to annotate
#' @param annotate_p pvalue threshold to annotate
#' @param title optional string for plot title
#' @param chrcolor1 first alternating color for chromosome
#' @param chrcolor2 second alternating color for chromosome
#' @param highlight_snp vector of snps to highlight
#' @param highlight_p pvalue threshold to highlight
#' @param highlighter color to highlight
#' @param groupcolors named vector of colors for data in 'PHE' or 'Group' column
#' @param chrblocks boolean, turns on x-axis chromosome marker blocks
#' @param background variegated or white
#' @param file file name of saved image
#' @param hgt height of plot in inches
#' @param wi width of plot in inches
#' @param res resolution of plot in pixels per inch
#' @import ggplot2
#' @return png image
#' @export
#' @family PheWAS functions
#' @family static plotting functions
#' @seealso \code{\link{ipheman}}, \code{\link{apheman}}, \code{\link{gman}}, \code{\link{eman}}
#' @examples
#' data(phewas)
#' pheman(d=phewas, title="PheWAS Example")

phegene <- function(top, bottom, phegroup, line, levs, log10=TRUE, yaxis, ymax, opacity=1, annotate_snp, annotate_p, highlight_snp, highlight_p, highlighter="red", title=NULL, chrcolor1="#AAAAAA", chrcolor2="#4D4D4D", groupcolors, background="variegated", chrblocks=TRUE, file="phegene", hgt=7, wi=12, res=300 ){
  
  #Sort data
  top$Data <- "Top"
  bottom$Data <- "Bottom"
  d <- rbind(top, bottom)
  d$POS <- as.numeric(as.character(d$POS))
  d$CHR <- factor(d$CHR, levels=levs)
  if(!missing(phegroup)){
    print("Only phenotypes with grouping information will be plotted")
    d_phe <- merge(phegroup, d, by="PHE")
    names(d_phe)[names(d_phe)=="Group"] <- "Color"
  } else {
    d_phe <- d
    names(d_phe)[names(d_phe)=="PHE"] <- "Color"
  }
  #d_order <- d_phe[order(d_phe$CHR, d_phe$POS), ]
  d_order <- d_phe[order(d_phe$CHR), ]
  #d_order$pos_index <- seq.int(nrow(d_order))
  d_order$pos_index <- as.numeric(d_order$CHR)
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
  
  #lims$nposmin <- (seq_along(lims$Color) - 1)
  #lims$nposmax <- seq_along(lims$Color)
  #lims$navg <- (lims$nposmin + lims$nposmax)/2
  
  #d_order2 <- merge(d_order, lims[, c(1,5,9,12)], by.x="CHR", by.y="Color")
  #d_order2$xpos <- ifelse( d_order2$posmax-d_order2$posmin <= 0, d_order2$nposmin, d_order2$nposmin + ((d_order2$pos_index-d_order2$posmin)/(d_order2$posmax-d_order2$posmin)))
  
  #d_order <-  d_order2[, c(2,1,3:5,10)]
  #colnames(d_order)[6] <- "pos_index"
  
  #names(lims) <- c("Color", "snpx", "px", "posx", "oposmin", "snpy", "py", "posy",
  #                 "oposmax", "oav", "shademap", "posmin", "posmax", 'av')
  
  #Set up colors
  nchrcolors <- nlevels(factor(lims$Color))
  if(!missing(groupcolors)){
    dcols <- c(rep(x=c(chrcolor1, chrcolor2), length.out=nchrcolors, each=1), "#FFFFFF", "#EBEBEB")
    names(dcols) <-c(levels(factor(lims$Color)), "shade_ffffff", "shade_ebebeb")
    newcols <- c(dcols, groupcolors)
  } else {
    ngroupcolors <- nlevels(factor(d_order$Color))
    if(ngroupcolors > 15){
      if (!requireNamespace(c("RColorBrewer"), quietly = TRUE)==TRUE) {
        stop("Please install RColorBrewer to add color attribute for more than 15 colors.", call. = FALSE)
      } else {
        getPalette = grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))
        newcols <- c(rep(x=c(chrcolor1, chrcolor2), length.out=nchrcolors, each=1), getPalette(ngroupcolors), "#FFFFFF", "#EBEBEB")
      }
    } else {
      pal <- pal <- c("#009292", "#920000", "#490092", "#db6d00", "#24ff24",
                      "#ffff6d", "#000000", "#006ddb", "#004949","#924900",
                      "#ff6db6", "#6db6ff","#b66dff", "#ffb6db","#b6dbff")
      newcols <- c(rep(x=c(chrcolor1, chrcolor2), length.out=nchrcolors, each=1), pal[1:ngroupcolors], "#FFFFFF", "#EBEBEB")
    }
    names(newcols) <-c(levels(factor(lims$Color)), levels(factor(d_order$Color)), "shade_ffffff", "shade_ebebeb")
  }
  
  #Info for y-axis
  if(log10==TRUE){
    d_order$pval <- ifelse(d_order$Data=="Top", -log10(d_order$pvalue), log10(d_order$pvalue))
    yaxislab <- expression(paste("+/-log"[10], "(p-value)", sep=""))
    #if(!missing(line)) {redline <- -log10(line)}
  } else {
    d_order$pval <- d_order$pvalue
    yaxislab <- yaxis
    #if(!missing(line)) {redline <- line}
  }
  
  #Allow more than 6 shapes
  if("Shape" %in% names(d)){
    allshapes <- c("\u25BC", "\u25B2",16,15,17,3,7,8,0:2,4:6,9:14,18:25,33:127)
    shapevector <- allshapes[1:nlevels(as.factor(d$Shape))]
  }
  
  #Theme options
  yaxismin <- min(d_order$pval)
  backpanel <- ifelse(background=="white", "NULL", "geom_rect(data = lims, aes(xmin = posmin-0.5, xmax = posmax+0.5, ymin = -Inf, ymax = Inf, fill=factor(shademap)), alpha = 0.5)" )
  
  #Start plotting
  p <- ggplot() + eval(parse(text=backpanel))
  #Add shape info if available
  if("Shape" %in% names(d)){
    p <- p + geom_point(data=d_order, aes(x=pos_index, y=pval, color=factor(Color), shape=factor(Shape)), alpha=opacity) + scale_shape_manual(values=shapevector)
  } else {
    p <- p + geom_point(data=d_order, aes(x=pos_index, y=pval, color=factor(Color)), alpha=opacity)
  }
  p <- p + scale_x_continuous(breaks=lims$av, labels=lims$Color, expand=c(0,0))
  if(chrblocks==TRUE){p <- p + geom_rect(data = lims, aes(xmin = posmin, xmax = posmax, ymin = -Inf, ymax = min(d_order$pval), fill=as.factor(Color)), alpha = 1)}
  #Add legend
  p <- p + scale_colour_manual(name = "Color", values = newcols) + scale_fill_manual(name = "Color", values = newcols, guides(alpha=FALSE))
  p <- p + theme(axis.text.x=element_text(angle=90), panel.grid.minor.x = element_blank(), panel.grid.major.x=element_blank(), axis.title.x=element_blank(), legend.position="bottom", legend.title=element_blank())
  #Highlight if given
  if(!missing(highlight_snp)){
    if("Shape" %in% names(d)){
      p <- p + geom_point(data=d_order[d_order$SNP %in% highlight_snp, ], aes(x=pos_index, y=pval, shape=Shape), colour=highlighter) + scale_shape_manual(values=shapevector)
      p <- p + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p <- p + geom_point(data=d_order[d_order$SNP %in% highlight_snp, ], aes(x=pos_index, y=pval), colour=highlighter)
    }
  }
  if(!missing(highlight_p)){
    if("Shape" %in% names(d)){
      p <- p + geom_point(data=d_order[d_order$pvalue < highlight_p, ], aes(x=pos_index, y=pval, shape=Shape), colour=highlighter) + scale_shape_manual(values=shapevector)
      p <- p + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p <- p + geom_point(data=d_order[d_order$pvalue < highlight_p, ], aes(x=pos_index, y=pval), colour=highlighter)
    }
  }
  if(!missing(annotate_p)){
    if (!requireNamespace(c("ggrepel"), quietly = TRUE)==TRUE) {
      print("Consider installing 'ggrepel' for improved text annotation")
      p <- p + geom_text(data=d_order[d_order$pvalue < annotate_p,], aes(pos_index,pval,label=SNP))
    } else {
      p <- p + ggrepel::geom_text_repel(data=d_order[d_order$pvalue < annotate_p,], aes(pos_index,pval,label=SNP))
    }
  }
  if(!missing(annotate_snp)){
    if (!requireNamespace(c("ggrepel"), quietly = TRUE)==TRUE){
      print("Consider installing 'ggrepel' for improved text annotation")
      p <- p + geom_text(data=d_order[d_order$SNP %in% annotate_snp,], aes(pos_index,pval,label=SNP))
    } else {
      p <- p + ggrepel::geom_text_repel(data=d_order[d_order$SNP %in% annotate_snp,], aes(pos_index,pval,label=SNP))
    }
  }
  #Add title and y axis title
  p <- p + ggtitle(title) + ylab(yaxislab)
  #Add pvalue threshold line
  if(!missing(line)){
    for(l in line){
      p <- p + geom_hline(yintercept = l, colour="red")
    }
  }
  
  #Theme
  if(!missing(ymax)){
    yaxismax <- ymax
  } else {
    yaxismax <- max(d_order$pval)
  }
  
  if(chrblocks==TRUE){
    p <- p+ylim(c(yaxismin,yaxismax))
  } else {
    p <- p+scale_y_continuous(limits=c(yaxismin, yaxismax), expand=expand_scale(mult=c(0,0.1)))
  }
  if(background=="white"){p <- p + theme(panel.background = element_rect(fill="white"))}
  p <- p + guides(fill="none") 
  #Save
  print(paste("Saving plot to ", file, ".png", sep=""))
  ggsave(p, filename=paste(file, ".png", sep=""), dpi=res, units="in", height=hgt, width=wi)
  print(p)
  
  return(p)
  
}
