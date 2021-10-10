###Get missingness per variable
#df is dataframe where each column is a variable and missing is coded as NA

vmiss <- function(df){

  getNA <- function(x){
    N <- length(x)
    NMISS <- length(x[is.na(x)])
    PCTMISS <- NMISS/N
    out <- rbind(NMISS=NMISS, PCTMISS=PCTMISS)
    return(out)
  }

  sumNA <- as.data.frame(t(sapply(df[, -1], getNA)))
  sumNA <- cbind(rownames(sumNA), data.frame(sumNA, row.names=NULL))
  names(sumNA) <- c("Variable", "NMISS", "PCTMISS")
  return(sumNA)
}

