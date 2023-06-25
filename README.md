# PoW model: a tool for reconstructing deep evolutionary history of viruses


### Overview

The Prisoner of War (PoW) model is written in R, allows for the reconstruction of time trees over deep timescales, and is based on a mechanistic evolutionary model that explains the time-dependent changes in the evolutionary rate estimates of viruses over time.

![Prisoner of War model](https://github.com/mg878/PoW_model/blob/main/PoW_sigmoid_TDRP.jpeg)

### Prerequisites and required packages

PoW model requires the posterior rate distribution from dated sequences in the dataset. It also requires an ultrametric distance tree of all sequences in the dataset using a standard HKY substitution model (with no assumed rate heterogenity) along with the posterior distributions of the base frequencies and transition-transverstion rate. All output files should be compatible with the BEAST 1 format. 

* `*.log` includes the posterior rate distribution (meanRate column in the log file) of the dates samples

* `*.log` includes the posterior base frequencies and transition-transversion rate distributions (frequencies1-4 and kappa columns in the log file) of the entire dataset using a standard HKY substitution model

* `*.trees` includes the posterior ultrametric distance trees 


PoW model is compatible with R version 4.0.0 or higher and is tested on 4.0.4 to 4.2.1.  It depends on several libraries:

* ape

* nleqslv

* stringr

### Guide

Ensure all required libraries are installed and the prerequisited files are available on your local directory.
Ensure that in the `PoW_mode.R`, set all the paths to the corresponding directories where the prerequisited files are located -- these can be found in the script lines #9, #13, #14, #95, #98, and #199.

### Developer info

  - Copyright and License: Mahan Ghafari, MIT Licence
  - Reference:
    * [A mechanistic evolutionary model explains the time-dependent pattern of substitution rates in viruses](https://doi.org/10.1016/j.cub.2021.08.020) by Mahan Ghafari, Peter Simmonds, Oliver G Pybus, and Aris Katzourakis. Current Biology (2021).
