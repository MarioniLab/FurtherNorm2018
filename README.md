# Further scaling normalization

## Overview

This repository contains code for further development of the scaling normalization methods implemented in the [_scran_](https://bioconductor.org/packages/scran) package.
It is based on the code at https://github.com/MarioniLab/Deconvolution2016, which accompanies the paper **Pooling across cells to normalize single-cell RNA sequencing data with many zero counts** by [Lun _et al._ (2016)](https://doi.org/10.1186/s13059-016-0947-7).

## Setting up

First, install the helper package in `package/` with:

```
R CMD INSTALL package/
```

... or some variant thereof.
It is also worth running:

```
BiocManager::install(ask=FALSE, version="devel")
```

... to ensure that the latest versions of all relevant packages are installed.

## Simulations

To run the simulations, enter the `simulations/` directory and run:

- `sim_noDE.R`, which simulates a variety of scenarios involving no DE between populations. 
- `sim_biDE.R`, which simulates a variety of scenarios involving DE between two populations. 
- `sim_multiDE.R`, which simulates a variety of simulations involving DE between three populations.

Note that the parallelization framework assumes a SLURM cluster with the `Rdevel` command (to run a version of R with BioC-devel packages).
This can be changed by modifying `createBatchParam.R` in `package/R` and/or `slurm.tmpl` in `package/inst/scripts`.

To summarize the simulation results, run `summarizer.R` to cluster the simulation scenarios based on the pattern of errors across all methods.

## Real data

To analyze real data, enter the `real` directory and run:

- `zeisel.Rmd`, to compute size factors for the Zeisel brain data set.
- `pbmc4k.Rmd`, to compute size factors for the PBMC 4K data set.
