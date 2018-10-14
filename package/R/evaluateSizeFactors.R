#' Evaluate normalization error
#'
#' Estimate the mean squared log-error (MSLE) of the size factors compared to their true values.
#' 
#' @param sf Numeric vector, estimated size factors for all cells.
#' @param truth Numeric vector, true size factors (from \code{\link{generateRawMeans}}).
#' @param is.de Logical or integer vector, which cells are from a separate population?
#' @param plot Logical, should the estimates be plotted against the truth.
#' @param main String, title for the plot.
#' @param col String or character vector, point colors for all cells.
#'
#' @details
#' This function will compute the log-fold change of \code{sf} from \code{truth}.
#' It will correct for any shift in location from zero using \code{\link{rlm}}, as these do not affect relative interpretation of the size factors.
#' If \code{is.de=NULL}, all cells are assumed to be from the same population.
#' The MSLE across all cells is then computed using the residuals of the fit.
#'
#' If \code{is.de} is specified, the shift in location is estimated from the cells that are \emph{not} in \code{is.de}.
#' This represents the \dQuote{baseline} cell population for which no DE genes were added.
#' The MSLE of this baseline population is estimated separately from the MSLE across all of the \code{is.de} cells (after correcting for the shift).
#' This enables separate evaluation of size factors in the presence of composition biases introduced by DE genes.
#' 
#' @return
#' A numeric scalar containing the MSLE for size factors of cells with no DE.
#' If \code{is.de} is specified, this is a numeric vector of length two, containing additionally the MSLE for the cells with DE.
#'
#' If \code{plot=TRUE}, a plot of the true and estimated size factors is generated on the current graphics device.
#'
#' @export
#' @importFrom MASS rlm
#' @importFrom stats residuals coef
#'
#' @author Aaron Lun
#'
#' @examples
#' truth <- 2^rnorm(100)
#' sf <- truth * 2^rnorm(100, 0.1)
#' evaluateSizeFactors(sf, truth)
#' 
#' is.de <- 1:20
#' sf[is.de] <- sf[is.de] * 1.1
#' evaluateSizeFactors(sf, truth, is.de)
evaluateSizeFactors <- function(sf, truth, is.de=NULL, plot=TRUE, main="", col="black") {
    logfold <- log2(sf) - log2(truth)

    if (is.null(is.de)) {
        fitted <- rlm(logfold ~ 1, na.action=na.exclude)
        resids <- residuals(fitted)
    } else {
        fitted <- rlm(logfold[-is.de] ~ 1, na.action=na.exclude)
        resids <- logfold - coef(fitted)
    }

    obs <- truth * 2^resids
    if (!is.null(is.de)) {
        non.DE.err <- 2^sqrt(mean(resids[-is.de]^2))-1
        DE.err <- 2^sqrt(mean(resids[is.de]^2))-1
        err.all <- c("non-DE"=non.DE.err, "DE"=DE.err)
    } else {
        non.DE.err <- 2^sqrt(mean(resids^2))-1
        err.all <- c("non-DE"=non.DE.err)
    }

    if (plot) {
        all.range <- range(range(truth), range(obs))
        col <- rep(col, length.out=length(sf))
        shuffle <- sample(length(sf))
 
        plot(truth[shuffle], obs[shuffle], xlim=all.range, ylim=all.range, 
            ylab="Estimated factors", xlab="True factors", log="xy", pch=16, 
            col=col[shuffle], cex.axis=1.5, cex.lab=1.8, main=main, cex.main=1.8)
        abline(0, 1, col="red")
    
        err.formatted <- format(round(err.all*100, 1), nsmall=1)
        legend("topleft", bty="n", cex=1.2,
            legend=paste0(names(err.formatted), " error = ", err.formatted, "%"))
    }
    return(err.all)
}

