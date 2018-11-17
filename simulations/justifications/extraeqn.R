# This shows that, even in the absense of an identifiable system,
# the extra equations can approach the correct solution.

library(FurtherNorm2018)
param <- generateRawMeans(ngenes=10000, ncells=200)
counts <- sampleCounts(param$fitted, param$dispersion)

# Adding unbalanced DE to induce composition biases.
# Otherwise the extra equations collapse to library size normalization,
# which would be trivially correct.
mus <- param$fitted
de.cell <- seq_len(50)
goes.up <- seq_len(1000)
mus[goes.up,de.cell] <- mus[goes.up,de.cell] * 10 

counts <- sampleCounts(mus, param$dispersion)

library(scran)
sf <- computeSumFactors(counts, size=seq(20, 100, by=5), min.mean=0.1)
evaluateSizeFactors(sf, param$sf, is.de=de.cell, col=c('black', "red"))

sf <- computeSumFactors(counts, size=seq(21, 101, by=5), min.mean=0.1)
evaluateSizeFactors(sf, param$sf, is.de=de.cell, col=c('black', "red"))

