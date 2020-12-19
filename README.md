# PoW_model

M. Ghafari et al.: Prisoner of War dynamics explains the time-dependent pattern of substitution rates in viruses


##Codes

`Evolve_sequences.nb`: Generating neutral and epistatic sites for the simulation analysis. 
Output files are (1) initial_sample.txt (2) final_sample.txt 
`neutral.nb`: Extracting and plotting the inferred root heights and clock rates from `Evolve_sequences.nb` for the purely neutral scenario, i.e. no epistatic sites. 
`epistatic.nb`: Extracting and plotting the inferred root heights and clock rates from `Evolve_sequences.nb` for the epistatic scenario.

*stats* folder includes the log files from BEAST 1.10 parameter inference, along with the inferred mean and median clock rates.
These files are named in the following format `*_mA_B` where letter A varies from 2 to 5 and corresponds to mutation rates ranging from 10^-2 to 10^-5, respectively and letter B

##Fasta files

`Genome_metadata.xlsx`: this contains the metadata from GISAID used for the phylogenetic analysis


