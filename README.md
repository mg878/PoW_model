# PoW_model

Ghafari et al.: A mechanistic evolutionary model explains the time-dependent pattern of substitution rates in viruses


## Codes

`Evolve_sequences.nb`: Generates simulated sequences under the Wright-Fisher model based on a neutral and/or compensatory evolution. 
Output files are (1) initial_sample.txt (first snapshot from the population after 10N_e generations) (2) final_sample.txt (final snapshot from the population after 10N_e + T generations where T is the time window of rate measurement)

`neutral.nb`: Extracts and plots the inferred root heights and clock rates from BEAST log file for the purely neutral scenario, i.e. no epistatic sites. 

`epistatic.nb`: Extracts and plots the inferred root heights and clock rates from BEAST log file for the epistatic scenario.

`PoW_ultrametric.R`: This is the code to generate PoW-transformed trees from the maximum clade credibility distance trees for FV, HCV, and sarbecovirus datasets.

***stats*** folder includes the log files for parameter inference, along with the inferred mean and median clock rates using BEAST.
The files in this folder are named in the following format: `*_mA_B` letter A varies from 2 to 5 and corresponds to substitution rates ranging from 3x10^-2 to 3x10^-5 and letter B takes the values corresponding to the time window of observation, i.e. T=10, 100, 1000, 10000, 100000, 1000000, 10000000

## Output

`*_&_%_distance.trees`: Contains the inferred trees produced by BEAST (v 1.10.4) for * = FV, HCV, and sarbecovirus datasets using & = JC69 or HKY85 substitution models. This includes % = fixedRate or variedRate whereby the short-term rate and fastest-evolving rate groups are either fixed or varied.

`*_&_%_PoWtransformed.trees`: Contains the PoW-transformed trees (produced by our R package).

`*_logfile.log`: Contains output log files produced by BEAST for * = HCV and sarbecovirus heterochronous datasets (i.e., standard substitution models used to infer the short-term substitution rates).

`HCV_GTR_*.txt`: Contains the inferred MCC time tree, XML and log file for the HCV samples using the GTR+G4 substitution model.

`HCV_alignments.FST`: Contains the HCV alignments.