allfluxes <- list()
for (ii in 1:10){
    path <- paste("b_", ii, "/minmax.csv", sep="")
    allfluxes[[ii]] <- read.csv(path, header=FALSE)
}

combinations <- combn(1:10, 2, simplify=FALSE)
for (jj in combinations){
    print(jj)
    print(sum(abs(allfluxes[[jj[1]]] - allfluxes[[jj[2]]])))
}
