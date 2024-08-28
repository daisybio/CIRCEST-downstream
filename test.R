source("load.R")

circ <- loadCirc()

# circ is a numeric dataframe
# Calculate median of each row
circ_median <- apply(circ, 1, median)

print(median(circ_median))