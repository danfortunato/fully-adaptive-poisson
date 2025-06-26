# Reproducibility repository for *A fully adaptive, high-order, fast Poisson solver for complex two-dimensional geometries*

This repository contains the code and data required to reproduce all results from the paper:

- D. Fortunato, D. B. Stein, and A. H. Barnett, *A fully adaptive, high-order, fast Poisson solver for complex two-dimensional geometries*.

Figures 4.1, 4.3&ndash;4.4, 5.1&ndash;5.4 and Table 5.1 in the paper contain numerical results. Each has a corresponding script in this repository which generates a figure or data.

# Requirements

- MATLAB R2021a or newer
- GCC 11 or newer

This repository uses git submodules. To initialize the submodules, either clone this repository with the `--recurse-submodules` flag or run `git submodule update --init --recursive` after cloning.

# Compilation

## Compile the fast coordinate rootfinder

1. Make sure MATLAB's `mex` is accessible from the command line. Edit the file `+util/make.inc` to set `MEX` and `CXX` appropriately.
2. Change to the `+util` directory and type `make` to build the MATLAB MEX interface.

## Compile the k-d tree library

1. Follow the instructions in `external/kdtree` to build the MATLAB MEX interface.

## Compile the FMM library

1. Follow the instructions in `external/fmmlib2d` to build the MATLAB MEX interface. This will likely involve modifying `external/fmmlib2d/matlab/mwrap.inc` to match your specific platform.

## Compile the box code library

1. Follow the instructions in `external/treefun/external/boxcode2d` to build the MATLAB MEX interface.

# Running

First, run the MATLAB script `setup.m` to add the required directories to your MATLAB path. To generate figure X.Y, change to the `figures` directory and run the MATLAB script `figureX_Y.m`. The script may plot a figure in MATLAB and/or write data to a file. Many scripts define a boolean variable `doExport`, which if set to `true` will export the generated figures as images.
