# Summarizes the simulation results by averaging within scenarios and clustering across scenarios. 

CLUSTERER <- function(indir, k=10) {
    all.files <- list.files(indir, pattern="\\.tsv$", full=TRUE)
    collected <- vector("list", length(all.files))
    for (f in seq_along(all.files)) {
        collected[[f]] <- colMeans(read.table(all.files[f], header=TRUE))
    }

    mat <- do.call(rbind, collected)
    clustered <- kmeans(mat, centers=k)

    # Showing the performance within each cluster.
    pdf(file.path(indir, "summary.pdf"))
    par(mar=c(8.1, 4.1, 4.1, 2.1))
    for (i in seq_len(k)) {
        heights <- clustered$centers[i,]
        X <- barplot(heights, main=paste("Cluster", i), las=2, ylab="Size factor error (%)")

        chosen <- clustered$cluster==i
        sd <- sqrt(matrixStats::colVars(mat[chosen,])) # SD, not SE: we care about the variability.
        segments(X, heights, X, heights + sd)
        segments(X - 0.1, heights + sd, X + 0.1, heights + sd)
    
        highest <- matrixStats::colMaxs(mat[chosen,])
        points(X, highest, pch=4)
    }
    dev.off()

    # Showing the characteristics of each cluster.
    current.files <- sub(".tsv$", "", basename(all.files))
    scenarios <- strsplit(current.files, "-")
    df <- data.frame(do.call(rbind, scenarios))
    colnames(df) <- c("Mode", "Ncells", "GenesUp", "GenesDown", "Prop", "Effect")
    for (i in colnames(df)[-1]) {
        df[[i]] <- factor(as.numeric(sub("^[a-z]+", "", df[[i]])))
    }

    pdf(file.path(indir, "clusters.pdf"), width=5, height=10)
    for (i in seq_len(k)) {
        chosen <- clustered$cluster==i
        cur.df <- df[chosen,]

        # Creating a dot plot of characteristics.
        plot(0,0,type="n", axes=FALSE, xlab="", ylab="", xlim=c(1, ncol(cur.df)), ylim=c(1, nrow(cur.df)), main=paste("Cluster", i))
        for (j in seq_len(ncol(cur.df))) {
            current <- cur.df[[j]]
            if (colnames(cur.df)[j]=="Mode") {
                points(rep(j, length(current)), seq_along(current), pch=ifelse(current=="UMI", 16, 1))
            } else {
                if (any(levels(current)=="0")) { 
                    cols <- c("white", rev(viridis::viridis(nlevels(current)-1)))
                } else {
                    cols <- rev(viridis::viridis(nlevels(current)))
                }
                points(rep(j, length(current)), seq_along(current), pch=21, col="grey80", bg=cols[as.integer(current)])
            }
        }
        axis(1, at=seq_len(ncol(cur.df)), labels=colnames(cur.df), las=2)
    }
    dev.off()
}


##############################################

indir <- "results_biDE"

for (
