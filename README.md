# PoW model: a tool for reconstructing deep evolutionary history of viruses


### Overview

The Prisoner of War (PoW) model is written in R, allows for the reconstruction of time trees over deep timescales, and is based on a mechanistic evolutionary model that explains the time-dependent changes in the evolutionary rate estimates of viruses over time.

![Prisoner of War model](https://github.com/mg878/PoW_model/blob/main/PoW_sigmoid_TDRP.jpeg)

### Prerequisites and required packages

PoW model requires the posterior rate distribution from dated sequences in the dataset. It also requires an ultrametric distance tree of all sequences in the dataset using a standard HKY substitution model (with no assumed rate heterogenity) along with the posterior distributions of the base frequencies and transition-transverstion rate. All output files should be compatible with the BEAST 1 format. 

PoW model is compatible with R version 4.0.0 or higher and is tested on 4.0.4 to 4.2.1.  It depends on several libraries:

* ape

* nleqslv

* stringr

### 

