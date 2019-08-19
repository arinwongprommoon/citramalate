library(ggplot2)
library(RColorBrewer)
library(viridis)
library(reshape2)

# Reads files
allfluxes <- list()
for (ii in 1:10){
    path <- paste("b_", ii, "/minmax.csv", sep="")
    allfluxes[[ii]] <- read.csv(path, header=FALSE)
}

# Calculates distance matrix
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

# Plots and colours distance matrix
distances[lower.tri(distances)] <- 0
colnames(distances) <- c('0.1~1.0', '0.1~2.0', '0.1~3.0', '0.1~4.0', '0.1~5.0', '0.1~6.0', '0.1~7.0', '0.1~8.0', '0.1~9.0', '0.1~10.0')
rownames(distances) <- c('0.1~1.0', '0.1~2.0', '0.1~3.0', '0.1~4.0', '0.1~5.0', '0.1~6.0', '0.1~7.0', '0.1~8.0', '0.1~9.0', '0.1~10.0')

longdata <- melt(distances)
longdata <- longdata[longdata$value!=0,]

plot <- ggplot(longdata, aes(x = Var2, y = Var1, label=sprintf("%0.2f", round(value, digits=2)))) +
    geom_raster(aes(fill=value)) +
    scale_fill_viridis() +
    geom_text() +
    labs(x="set", y="set", title="Distance Matrix") +
    theme_bw()
