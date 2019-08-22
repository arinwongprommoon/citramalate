library(reshape2)
library(tidyverse)

# Purpose: creates list of six reactions to vary for each kinetic model reaction in DE, based on one-reaction lists

# Reads files and puts them into a list of data frames
temp = list.files(path = "./VMAXDATA/", pattern = "*.csv", full.names = T)
reaclists = lapply(temp, read.csv, header = FALSE)

# Creates empty data frame: first column is the reaction, and the rest functions as the list of reactions to vary in DE
setdf <- data.frame(reaction = character(),
		    DE1 = character(),
		    DE2 = character(),
		    DE3 = character(),
		    DE4 = character(),
		    DE5 = character(),
		    DE6 = character(),
		    stringsAsFactors = FALSE)
# Populates data frame, each element is each one-reaction list for each kinetic model reaction
for (ii in 1:length(reaclists)) {
	# Gets name of reaction from temp; temp stores paths
	r = substr(temp[ii], 22, nchar(temp[ii])-4)
	# Sort by second column (differences), gets first six reactions (first column)
	l = as.vector(reaclists[[ii]][order(-reaclists[[ii]]$V2),][1:6,1])
	# Adds new row in data frame
	setdf[nrow(setdf)+1,] = c(r,l)
}
# Writes data frame to file
write.csv(setdf, file="fluxdatasets.csv")
