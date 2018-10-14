#' Parallelization parameters 
#' 
#' Creates a custom \linkS4class{BatchtoolsParam} object for common use across scripts. 
#'
#' @param N Integer, number of workers.
#' 
#' @details
#' This function assumes a SLURM cluster with the command \code{Rdevel}.
#' It is intended for development environments where installation of new versions on \code{R} itself would compromise release projects.
#'
#' @return
#' A \linkS4class{BatchtoolsParam} object deploying to the specified number of workers.
#'
#' @importFrom BiocParallel BatchtoolsParam
#'
#' @export
#' @author Aaron Lun
createBatchParam <- function(N) {
	template <- system.file("scripts", "slurm.tmpl", package="FurtherNorm2018")
    BatchtoolsParam(N,
        cluster="slurm", template=template,
        logdir="parallel", log=TRUE,
        RNGseed=10000L,
        resources=list(walltime=20000, memory=8000, ncpus=1))
}

