allfluxes <- list()
for (ii in 1:8){
    path <- paste("c_", ii, "/minmax.csv", sep="")
    allfluxes[[ii]] <- read.csv(path, header=FALSE)
}

distances <- matrix(1:64, nrow=8, ncol=8)
for (kk in 1:8){
    distances[kk, kk] <- 0
}
combinations <- combn(1:8, 2, simplify=FALSE)
options(digits = 3)
for (jj in combinations){
    print(jj)
    dd = (sum(abs(allfluxes[[jj[1]]] - allfluxes[[jj[2]]])))
    print(dd)
    distances[jj[1], jj[2]] <- dd
    distances[jj[2], jj[1]] <- dd
}
print(distances)
