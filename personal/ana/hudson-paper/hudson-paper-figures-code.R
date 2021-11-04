###### Figure 1A ######
system("gunzip Okbay_27225129-EduYears_Main_Women.txt.gz")
system("gunzip Okbay_27225129-EduYears_Main_Men.txt.gz")
women <- data.table::fread("Okbay_27225129-EduYears_Main_Women.txt")
men <- data.table::fread("Okbay_27225129-EduYears_Main_Men.txt")
cols <- c("MarkerName", "CHR", "POS", "Pval")

top <- as.data.frame(women[, ..cols])
bottom <- as.data.frame(men[, ..cols])

names(top) <- c("SNP", "CHR", "POS", "pvalue")
names(bottom) <- names(top)
hudson::gmirror(top = top,
                bottom = bottom,
                tline = 0.05/nrow(women),
                bline = 0.05/nrow(men),
                highlight_p = c(0.05/nrow(women),
                                0.05/nrow(men)),
                highlighter = "green",
                annotate_snp = c("rs7613360", "rs11712056"),
                toptitle = "Okbay et al Educational Attainment: Women",
                bottomtitle = "Okbay et al Educational Attanment: Men",
                file="Figure_1A")

###### Figure 1B ######
system("gunzip AGEN_lipids_hapmap_hdl_m2.txt.gz")
system("gunzip AGEN_lipids_hapmap_ldl_m2.txt.gz")
system("gunzip AGEN_lipids_hapmap_tg_m2.txt.gz")
system("gunzip AGEN_lipids_hapmap_tc_m2.txt.gz")
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

top <- rbind(hdl, ldl, tc, tg)
bottom <- rbind(women, men)

set.seed(42); top <- dplyr::sample_n(top[top$pvalue < 0.005, ], size = 100000)
set.seed(42); bottom <- dplyr::sample_n(bottom[top$pvalue < 0.005, ], size = 40000)
phenos <- c("synonymous", "stop lost", "stop gained",
            "splice donor", "frameshift", "missense")

# make some random gene names and assign to df
genes <- c("PPARGC1A","SIRT6","CDKN1A","HMGB1","FGF23","STAT3","RPA1","GSS",
           "GSTP1","MLH1","MAX","ATM","ARHGAP1","EGF","ELN","APOC3","TAF1",
           "RAE1","IL6","APP","EEF2","IGF1","MAPK3","POLB","SIN3A","SOD1",
           "HSP90AA1","STAT5A","TP73","PTK2","TERF1","DBN1","PDGFRA","GSK3B",
           "AIFM1","POU1F1","SIRT7","CDKN2A","HIF1A","IL2RG", "PRDX1","CSNK1E",
           "MIF","INS","GDF11","IFNB1","E2F1","PMCH","SUN1","CACNA1A","POLA1",
           "MAP3K5","GRN","MT1E","MTOR","ARNTL","APOE","MYC","GCLC","VCP",
           "PRKDC","LRP2","PTPN1","RET","TP53BP1","TERF2","CISD2","BRCA2",
           "UBE2I","GRB2","IL7","TERC","TOP2A","HTRA2","HDAC1","BAK1","AGTR1",
           "EEF1A1","HSPA1B","GPX1","SST","NUDT1","GPX4","PDPK1","FOXM1","LEPR",
           "BMI1","SERPINE1","IL7R","NFE2L1","TFAP2A","RAD51","PTPN11","CLU",
           "CYP2C19", "VEGFA","NFE2L2","PDGFRB","TPP2","PCK1")
set.seed(42); top$CHR <- sample(genes, nrow(top), replace = TRUE)
set.seed(42); top$PHE <- sample(phenos, nrow(top), replace = TRUE)
set.seed(42); bottom$CHR <- sample(genes, nrow(bottom), replace = TRUE)
set.seed(42); bottom$PHE <- sample(phenos, nrow(bottom), replace = TRUE)

# subset 
top <- top[, c("PHE", "MARKERNAME", "CHR", "POS", "P")]
bottom <- bottom[, c("PHE", "MarkerName", "CHR", "POS", "Pval")]
names(top)[c(2,5)] <- c("SNP", "pvalue")
names(bottom)[c(2,5)] <- c("SNP", "pvalue")

phegene(top = as.data.frame(top),
        bottom = as.data.frame(bottom),
        levs = genes, 
        title = "Gene Based PheWAS Plot Example",
        file = "Figure_1D")

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
