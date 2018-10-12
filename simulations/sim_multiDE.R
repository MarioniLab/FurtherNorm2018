# This checks the performance of normalization methods in the presence of multiple subpopulations. 
# First we define a function that simulates the data and evaluates each method.

FUN <- function(mode, overlap, genes.up, genes.down, prop, effect, dir.out) {
    source("functions.R")

    # Skipping impossible scenarios.
    if (genes.up==0 && genes.down==0) {
        return(NULL)
    } else if (genes.up==0 && effect!=2) {
        return(NULL)
    }
    
    prefix <- sprintf("%s-%s-u%i-d%i-p%.1f-e%i", mode, overlap, genes.up, genes.down, prop, effect)
    stub <- file.path(dir.out, prefix)
    fout <- paste0(stub, ".tsv")

    for (it in 1:10) { 
        ncells <- 1000
        param <- generateRawMeans(ncells=ncells, mode=mode)
        truth <- param$sf

        # Adding the requested DE (second population has reverse profile from the first population).
        mus <- param$fitted
    
        de.cell1 <- seq_len(ncells*prop)
        goes.up1 <- seq_len(genes.up)
        goes.down1 <- seq_len(genes.down) + genes.up
        mus[goes.up1,de.cell1] <- mus[goes.up1,de.cell1] * effect
        mus[goes.down1,de.cell1] <- 0
    
        de.cell2 <- seq_len(ncells*prop) + max(de.cell1)
        goes.up2 <- seq_len(genes.down) 
        if (!overlap) { # Do the DE genes in the second population overlap with the first?
            goes.up2 <- goes.up2 + max(goes.down1)
        }
        goes.down2 <- seq_len(genes.up) + max(goes.up2)
        mus[goes.up2,de.cell2] <- mus[goes.up2,de.cell2] * effect
        mus[goes.down2,de.cell2] <- 0
    
        de.cell <- c(de.cell1, de.cell2)
        counts <- sampleCounts(mus, param$dispersion)
        output <- runAllMethods(counts)

        # Assessing the different methods. 
        collected <- vector("list", length(output))
        names(collected) <- names(output)
        if (it==1L) {
            # Generating an output plot for this simulation.
            col <- rep("black", ncells)
            col[de.cell1] <- "orange"
            col[de.cell2] <- "dodgerblue"

            pdf(paste0(stub, ".pdf"))
            par(mar=c(5.1,5.1,4.1,1.1))
            for (x in names(output)) {
                collected[[x]] <- evaluateSizeFactors(output[[x]], truth, de.cell, main=x, col=col)["DE"]
            }
            dev.off()
        } else {
            for (x in names(output)) {
                collected[[x]] <- evaluateSizeFactors(output[[x]], truth, de.cell, plot=FALSE)["DE"]
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

source("functions.R")
all.scenarios <- spawnScenarios(
    mode=c("UMI", "read"), 
    overlap=c(TRUE, FALSE),
    genes.up=c(0, 1000, 2000), 
    genes.down=c(0, 1000, 2000),
    effect=c(2, 5, 10),
    prop=c(0.1, 0.2, 0.3))

dir.out <- "results_multiDE"
dir.create(dir.out, showWarnings=FALSE)

library(BiocParallel)
BPPARAM <- createBatchParam(20)
X <- do.call(bpmapply, c(list(FUN=FUN, SIMPLIFY=FALSE, USE.NAMES=FALSE, BPPARAM=BPPARAM, MoreArgs=list(dir.out=dir.out)), all.scenarios))
