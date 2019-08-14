library(MASS)
od600 <- read.csv("orth-2011-supplementary-table-1.csv")
for(ii in colnames(od600)){
       print(ii)
       fit <- fitdistr(od600[,ii], "normal")
       print(fit$estimate)
}