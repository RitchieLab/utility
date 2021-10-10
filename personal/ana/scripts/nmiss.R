###Get missingness per sample
#df is dataframe where each column is a variable and missing is coded as NA

nmiss <- function(df){
  N <- ncol(df[, -1])
  df$NMISS <- rowSums(!is.na(df[, -1]))
  df$PCTMISS <- (N-df$NMISS)/N
  return(df[, colnames(df) %in% c("ID", "NMISS", "PCTMISS")])
}

