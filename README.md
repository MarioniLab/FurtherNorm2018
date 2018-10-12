# Further MNN algorithm development

## Overview

This repository contains code for further development of the scaling normalization methods, as implemented in the `computeSumFactors` and `quickSumFactors` functions in the [_scran_](https://bioconductor.org/packages/scran) package.
It is based on the code at https://github.com/MarioniLab/Deconvolution2016, which accompanies the paper **Pooling across cells to normalize single-cell RNA sequencing data with many zero counts** by [Lun _et al._ (2016)](https://doi.org/10.1186/s13059-016-0947-7).

## Simulations

To run the simulations, enter the `simulations/` directory and run:

- `sim_noDE.R`, which simulates a variety of scenarios involving no DE between populations. 
- `sim_biDE.R`, which simulates a variety of scenarios involving DE between two populations. 
- `sim_multiDE.R`, which simulates a variety of simulations involving DE between three populations.

Note that the parallelization framework assumes a SLURM cluster. 
This can be changed by modifying `createBatchParam` in `functions.R`.
