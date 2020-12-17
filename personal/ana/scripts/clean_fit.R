
#fit <- lapply(clean_pheno[, c(cov, "PHENO", names(clean_pheno)[5:16])], 
#                         function (x) return(tryCatch(do.call("glm", 
#                                                              list(as.formula(fmla), 
#                                                                   family=as.name("binomial"), 
#                                                                   data=as.name("clean_pheno"))), 
#                                                      error=function(e) NULL)))

clean_fit <- function(fit, phenotype, med_suffix){
  nmco<- fit[!sapply(fit, is.null)]
  sco <- lapply(nmco, function (x) summary(x))
  #Grab sample size, beta, se beta, and pvalue
  rco <- data.frame(t(as.data.frame(sapply(nmco, function(x) as.data.frame(length(x$residuals))))),
                    t(as.data.frame(sapply(nmco, function(x) as.data.frame(x$converged)))),
                    t(as.data.frame(sapply(sco, function(x) as.data.frame(cbind(x$coefficients[2,1],
                                                                                x$coefficients[2,2],
                                                                                x$coefficients[2,4]))))),
                   t(as.data.frame(sapply(nmco, function(x) as.data.frame(exp(cbind("Odds ratio" = coef(x), 
                                                                                    confint.default(x, level = 0.95)))[2,]))))
)
  prco <- data.frame(names = gsub("\\.length.x.residuals.","", row.names(rco)), rco, row.names = NULL)
  names(prco) <- c("Variable", "N", "Converged", "Beta", "SE", "Variable_pvalue","OR", "OR_LB", "OR_UB")
  prco <- as.data.frame(lapply(prco, unlist))
  prco[order(prco$Variable_pvalue), -1]
  prco$Phenotype <- phenotypes
  prco$Med <- gsub(med_suffix, "", prco$Variable)
  return(prco)
}
