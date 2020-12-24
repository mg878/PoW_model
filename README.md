# PoW_model

M. Ghafari et al.: Prisoner of War dynamics explains the time-dependent pattern of substitution rates in viruses


## Codes

`Evolve_sequences.nb`: Generates simulated sequences under the Wright-Fisher model based on a neutral and/or compensatory evolution. 
Output files are (1) initial_sample.txt (first snapshot from the population after 10N_e generations) (2) final_sample.txt (final snapshot from the population after 10N_e + T generations where T is the time window of rate measurement)

`neutral.nb`: Extracts and plots the inferred root heights and clock rates from BEAST log file for the purely neutral scenario, i.e. no epistatic sites. 

`epistatic.nb`: Extracts and plots the inferred root heights and clock rates from BEAST log file for the epistatic scenario.

***stats*** folder includes the log files for parameter inference, along with the inferred mean and median clock rates using BEAST.
The files in this folder are named in the following format: `*_mA_B` letter A varies from 2 to 5 and corresponds to substitution rates ranging from 3x10^-2 to 3x10^-5 and letter B takes the values corresponding to the time window of observation, i.e. T=10, 100, 1000, 10000, 100000, 1000000, 10000000

## Nexus

`sarbeco_ultimate_*.txt`: The PoW-transformed branch lengths for the Maximum Clade Credibility (MCC) time tree corresponding to * = median, lower, and upper quantile Sarbecovirus substitution rate estimates

`sarbeco_JC69_*.txt`: Contains the inferred MCC distance tree, XML and log files using for the Sarbecovirus samples using the Jukes-Cantor (JC69) substitution model.

`HCV_ultimate_*.txt`: The PoW-transformed branch lengths for the maximum clade credibility tree corresponding to * = median, lower, and upper quantile HCV substitution rate estimates

`HCV_JC69_*.txt`: Contains the inferred MCC distance tree, XML and log files using for the HCV samples using the JC69 substitution model.

`HCV_GTR_*.txt`: Contains the inferred MCC time tree and log file using for the HCV samples using the GTR+G4 substitution model.

`HCV_alignments.FST`: Contains the HCV alignments.