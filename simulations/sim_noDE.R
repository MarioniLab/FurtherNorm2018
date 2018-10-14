# This checks the performance of normalization methods in the presence of no DE genes.
# First we define a function that simulates the data and evaluates each method.

FUN <- function(mode, ncells, dir.out) {
    library(FurtherNorm2018)

    prefix <- sprintf("%s-n%i", mode, ncells)
    stub <- file.path(dir.out, prefix)
    fout <- paste0(stub, ".tsv")

    for (it in 1:10) { 
        param <- generateRawMeans(ncells=ncells, mode=mode)
        truth <- param$sf

        # Adding the requested DE.
        mus <- param$fitted
        counts <- sampleCounts(mus, param$dispersion)
        output <- runAllMethods(counts)
    
        collected <- vector("list", length(output))
        names(collected) <- names(output)
        if (it==1L) {
            # Generating an output plot for this simulation.
            pdf(paste0(stub, ".pdf"))
            par(mar=c(5.1,5.1,4.1,1.1))
            for (x in names(output)) {
                collected[[x]] <- evaluateSizeFactors(output[[x]], truth, main=x)["non-DE"]
            }
            dev.off()
        } else {
            for (x in names(output)) {
                collected[[x]] <- evaluateSizeFactors(output[[x]], truth, plot=FALSE)["non-DE"]
            }
        }

        # Writing summary statistics to table.
        stats <- do.call(cbind, collected)
        write.table(file=fout, data.frame(format(stats, digits=3)),
                append=it!=1, quote=FALSE, sep="\t", col.names=it==1, row.names=FALSE)
    }
    return(NULL)
}

##############################################

library(FurtherNorm2018)
all.scenarios <- spawnScenarios(
    mode=c("UMI", "read"), 
    ncells=c(100, 200, 500, 1000))

dir.out <- "results_noDE"
dir.create(dir.out, showWarnings=FALSE)

library(BiocParallel)
BPPARAM <- createBatchParam(20)
X <- do.call(bpmapply, c(list(FUN=FUN, SIMPLIFY=FALSE, USE.NAMES=FALSE, BPPARAM=BPPARAM, MoreArgs=list(dir.out=dir.out)), all.scenarios))
