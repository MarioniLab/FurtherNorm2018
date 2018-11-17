# This calculates the MSE of the size factor estimates from different 'sizes'.

require(scran)
require(FurtherNorm2018)
set.seed(100)
collected.20 <- collected.10 <- collected.2 <- collected.large <- collected.small <- list()      

for (it in 1:10) {
    param <- generateRawMeans(ngenes=10000, ncells=200)
    counts <- sampleCounts(param$fitted, param$dispersion)

    # Size factors with clustering prior to summation:
    sf <- computeSumFactors(counts, sizes=1:5*20+1, min.mean=0.1)
    collected.20[[it]] <- evaluateSizeFactors(sf, param$sf, plot=FALSE) 
    sf <- computeSumFactors(counts, sizes=2:10*10+1, min.mean=0.1)
    collected.10[[it]] <- evaluateSizeFactors(sf, param$sf, plot=FALSE) 
    sf <- computeSumFactors(counts, sizes=10:50*2+1, min.mean=0.1)
    collected.2[[it]] <- evaluateSizeFactors(sf, param$sf, plot=FALSE) 

    # Size factors using only small or large pool sizes.
    # Note, this must be co-prime with 'ncells'!
    sf <- computeSumFactors(counts, size=101, min.mean=0.1)
    collected.large[[it]] <- evaluateSizeFactors(sf, param$sf, plot=FALSE) 
    sf <- computeSumFactors(counts, sizes=21, min.mean=0.1)
    collected.small[[it]] <- evaluateSizeFactors(sf, param$sf, plot=FALSE) 
}

# MSE.
summary(unlist(collected.20))
summary(unlist(collected.10))
summary(unlist(collected.2))
summary(unlist(collected.large))
summary(unlist(collected.small))
