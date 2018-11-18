# This simulation aims to demonstrate the benefits of pre-scaling
# on the precision of the deconvolution:

set.seed(12345)
mode <- "UMI"
ncells <- 200

library(FurtherNorm2018)
param <- generateRawMeans(ncells=ncells, mode=mode)
truth <- param$sf
mus <- param$fitted
counts <- sampleCounts(mus, param$dispersion)

# WITHOUT:
libsizes <- colSums(counts)
sphere <- scran:::.generateSphere(libsizes)
ave.cell <- scater::calcAverage(counts)
keep <- ave.cell > 0
ave.cell <- ave.cell[keep]
counts <- counts[keep,]

new.sys <- scran:::.create_linear_system(counts, ave.cell, sphere, c(21)) 
final <- Matrix::qr.coef(Matrix::qr(new.sys$design), new.sys$output)
evaluateSizeFactors(final, truth, NULL, plot=FALSE)

# WITH:
new.sys <- scran:::.create_linear_system(t(t(counts)/libsizes), ave.cell, sphere, c(21))
final <- Matrix::qr.coef(Matrix::qr(new.sys$design), new.sys$output)
evaluateSizeFactors(final * libsizes, truth, NULL)

###################################################
# ... as well as a situation where it fails:

set.seed(1000)
mode <- "UMI"
ncells <- 100
genes.up <- 1000
genes.down <- 0
prop <- 0.1
effect <- 100 # huge effect, massive composition bias.

param <- generateRawMeans(ncells=ncells, mode=mode)
truth <- param$sf

mus <- param$fitted
de.cell <- seq_len(ncells*prop)
goes.up <- seq_len(genes.up)
goes.down <- seq_len(genes.down) + genes.up
mus[goes.up,de.cell] <- mus[goes.up,de.cell] * effect
mus[goes.down,de.cell] <- 0
counts <- sampleCounts(mus, param$dispersion)

library(scran)
X <- computeSumFactors(counts, min.mean=0)
evaluateSizeFactors(X, truth, de.cell,col=rep(c("red", "blue"), c(10, 90)))

# We can mitigate this by iterating and using the first-round size factors as t_j.
Y <- computeSumFactors(counts, min.mean=0, scaling=X)
evaluateSizeFactors(Y, truth, de.cell,col=rep(c("red", "blue"), c(10, 90)))
