allfluxes <- list()
for (ii in 1:10){
    path <- paste("b_", ii, "/minmax.csv", sep="")
    allfluxes[[ii]] <- read.csv(path, header=FALSE)
}

distances <- matrix(1:100, nrow=10, ncol=10)
for (kk in 1:10){
    distances[kk, kk] <- 0
}
combinations <- combn(1:10, 2, simplify=FALSE)
options(digits = 3)
for (jj in combinations){
    print(jj)
    dd = (sum(abs(allfluxes[[jj[1]]] - allfluxes[[jj[2]]])))
    print(dd)
    distances[jj[1], jj[2]] <- dd
    distances[jj[2], jj[1]] <- dd
}
print(distances)
