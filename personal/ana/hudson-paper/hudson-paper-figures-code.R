###### Figure 1A ######
sbp <- data.table::fread("27618452-GCST006259-EFO_0006335.h.tsv")
dbp <- data.table::fread("27618452-GCST006258-EFO_0006336.h.tsv")
cols <- c("hm_rsid", "hm_chrom", "hm_pos", "p_value")

top <- as.data.frame(sbp[!is.na(sbp$hm_chrom), ..cols])
bottom <- as.data.frame(dbp[!is.na(dbp$hm_chrom), ..cols])
names(top) <- c("SNP", "CHR", "POS", "pvalue")
names(bottom) <- names(top)

hudson::gmirror(top = top,
                bottom = bottom,
                tline = 0.05/nrow(top),
                bline = 0.05/nrow(bottom),
                highlight_p = c(0.05/nrow(top),
                                0.05/nrow(bottom)),
                highlighter = "green",
                annotate_snp = "rs11105354",
                toptitle = "Ehret et al Blood Pressure: SBP",
                bottomtitle = "Ehret et al Blood Pressure: DBP",
                background = "white",
                file="Figure_1A")

###### Figure 1B ######
hdl <- data.table::fread("AGEN_lipids_hapmap_hdl_m2.txt")
ldl <- data.table::fread("AGEN_lipids_hapmap_ldl_m2.txt")
tc <- data.table::fread("AGEN_lipids_hapmap_tc_m2.txt")
tg <- data.table::fread("AGEN_lipids_hapmap_tg_m2.txt")
hdl$PHE <- "HDL"; ldl$PHE <- "LDL"; tc$PHE <- "TC"; tg$PHE <- "TG"

cols1 <- c("PHE", "CHR", "MARKERNAME", "POS", "P")
cols2 <- c("PHE", "CHR", "MARKERNAME", "POS", "BETA")

top <- as.data.frame(rbind(hdl[hdl$P < 0.005, ..cols1],
                           ldl[ldl$P < 0.005, ..cols1],
                           tc[tc$P < 0.005, ..cols1],
                           tg[tg$P < 0.005, ..cols1]))
top$P <- -log10(top$P)

bottom <- as.data.frame(rbind(hdl[hdl$P < 0.005, ..cols2],
                              ldl[ldl$P < 0.005, ..cols2],
                              tc[tc$P < 0.005, ..cols2],
                              tg[tg$P < 0.005, ..cols2]))

names(top) <- c("PHE", "CHR", "SNP", "POS", "pvalue")
names(bottom) <- names(top)

ntests <- max(c(nrow(hdl), nrow(ldl), nrow(tc), nrow(tg)))

hudson::phemirror(top = top, 
                  bottom = bottom, 
                  tline = -log10(0.05/ntests), 
                  toptitle = "Spracklen et. al. Lipids (p-values)",
                  bottomtitle = "Spracklen et. al. Lipids (betas)",
                  file = "Figure_1B",
                  chrblocks = FALSE,
                  freey = TRUE,
                  log10 = FALSE,
                  yaxis = c(expression(paste("-log"[10], "(p-value)", sep="")),
                            "beta"))

###### Figure 1C ######
library(hudson)
data(ewas.t)
data(ewas.b)
emirror(top = ewas.t, 
        bottom = ewas.b, 
        annotate_p = c(3e-5, 1e-5),
        highlight_p = c(3e-5, 1e-5),
        highlighter = "cyan",
        toptitle = "Example Data Comparison: 1",
        bottomtitle = "Example Data Comparison: 2",
        file = "Figure_1C")

###### Figure 1D ######
devtools::source_url("https://raw.githubusercontent.com/RitchieLab/utility/master/personal/ana/scripts/phegene_hudson.R")

###### Figure 2A ######
df <- rbind(hdl[hdl$P < 0.005, ],
            ldl[ldl$P < 0.005, ], 
            tc[tc$P < 0.005, ], 
            tg[tg$P < 0.005, ])
df$Link <- paste0("https://www.ncbi.nlm.nih.gov/snp/", df$MARKERNAME)
df$Hover <- paste0("SNP: ", df$MARKERNAME,
                   "\nPOS: ", df$CHR, ":", df$POS,
                   "\nAllele: ", df$EFF_ALLELE, 
                   "\np-value: ", formatC(df$P, format = "e", digits = 3),
                   "\nBeta: ", signif(df$BETA, digits = 2),
                   "\nPhenotype: ", df$PHE,
                   "\nN: ", df$N)
names(df)[names(df)=="MARKERNAME"] <- "SNP"
names(df)[names(df)=="P"] <- "pvalue"
cols <- c("PHE", "SNP", "CHR", "POS", "pvalue", "Link", "Hover")
top  <- as.data.frame(df[df$BETA >= 0, ..cols])
bottom <- as.data.frame(df[df$BETA < 0, ..cols])

hudson::iphemirror(top = top,
                   bottom = bottom,
                   toptitle = "Spracklen et al Lipids: Beta >= 0",
                   bottomtitle = "Spracklen et al Lipids: Beta < 0",
                   file = "Figure_2A")
