# PoW_model

M. Ghafari et al.: Prisoner of War dynamics explains the time-dependent pattern of substitution rates in viruses


##Codes

`Evolve_sequences.nb`: Generating neutral and epistatic sites for the simulation analysis. 
Output files are (1) initial_sample.txt (2) final_sample.txt 
`neutral.nb`: Extracting and plotting the inferred root heights and clock rates from `Evolve_sequences.nb` for the purely neutral scenario, i.e. no epistatic sites. 
`epistatic.nb`: Extracting and plotting the inferred root heights and clock rates from `Evolve_sequences.nb` for the epistatic scenario.

*stats* folder includes the log files from BEAST 1.10 parameter inference, along with the inferred mean and median clock rates.
The files in this folder are named in the following format: `*_mA_B` letter A varies from 2 to 5 and corresponds to substitution rates ranging from 3x10^-2 to 3x10^-5 and letter B takes the values correspodning to the time window of observation, i.e. T*=10, 100, 1000, 10000, 100000, 1000000, 10000000

##Nexus

`HCV_ultimate_mean.txt`: The PoW-transformed branch lengths for the maximum clade credibility tree corresponding to mean substitution rate = and maximum rate group = 
`HCV_ultimate_low.txt`: The PoW-transformed branch lengths for the maximum clade credibility tree corresponding to mean substitution rate = and maximum rate group = 
`HCV_ultimate_high.txt`: The PoW-transformed branch lengths for the maximum clade credibility tree corresponding to mean substitution rate = and maximum rate group = 
`sarbeco_ultimate_mean.txt`: The PoW-transformed branch lengths for the maximum clade credibility tree corresponding to mean substitution rate = and maximum rate group = 
`sarbeco_ultimate_low.txt`: The PoW-transformed branch lengths for the maximum clade credibility tree corresponding to mean substitution rate = and maximum rate group = 
`sarbeco_ultimate_high.txt`: The PoW-transformed branch lengths for the maximum clade credibility tree corresponding to mean substitution rate = and maximum rate group = 
