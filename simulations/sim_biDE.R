# This checks the brittleness of the summation approach in the presence of DE genes.

library(BiocParallel)
source("functions.R")
dir.out <- "results_biDE"
dir.create(dir.out, showWarnings=FALSE)

FUN <- function(mode, ncells, genes.up, genes.down, prop, effect) {
    # Skipping impossible scenarios.
    if (genes.up==0 && genes.down==0) {
        return(NULL)
    } else if (genes.up==0 && effect!=2) {
        return(NULL)
    }
    
    prefix <- sprintf("%s-n%i-u%i-d%i-p%.1f-e%i", mode, ncells, genes.up, genes.down, prop, effect)
    stub <- file.path(cur.dir, prefix)
    fout <- paste0(stub, ".tsv")

    for (it in 1:10) { 
        param <- generateRawMeans(ncells=ncells, mode=mode)
        truth <- param$sf

        # Adding the requested DE.
        mus <- param$fitted
        de.cell <- seq_len(ncells*prop)
        goes.up <- seq_len(genes.up)
        goes.down <- seq_len(genes.down) + genes.up
    
        mus[goes.up,de.cell] <- mus[goes.up,de.cell] * effect
        mus[goes.down,de.cell] <- 0
    
        counts <- sampleCounts(mus, param$dispersion)
        output <- runAllMethods(counts)
    
        if (it==1L) {
            # Generating an output plot for this simulation.
            collected <- vector("list", length(output))
            names(collected) <- names(output)
            col <- rep("black", ncells)
            col[de.cell] <- "red"
        
            pdf(paste0(stub, ".pdf"))
            par(mar=c(5.1,5.1,4.1,1.1))
            for (x in names(output)) {
                collected[[x]] <- makeSFPlot(output[[x]], truth, de.cell, main=x, col=col)["DE"]
            }
            dev.off()
        }

        # Writing summary statistics to table.
        stats <- do.call(cbind, collected)
        write.table(file=fout, data.frame(format(stats, digits=3)),
                append=it!=1, quote=FALSE, sep="\t", col.names=it==1, row.names=FALSE)
    }
    return(NULL)
}

all.scenarios <- spawnScenarios(
    mode=c("UMI", "read"), 
    ncells=c(100, 200, 500, 1000),
    genes.up=c(0, 1000, 2000), 
    genes.down=c(0, 1000, 2000),
    effects=c(2, 5, 10),
    prop=c(0.1, 0.2, 0.5))

BPPARAM <- createBatchParam(20)
X <- do.call(bpmapply, c(list(FUN=FUN, SIMPLIFY=FALSE, USE.NAMES=FALSE, BPPARAM=BPPARAM), all.scenarios))
