# PoW_model

M. Ghafari et al.: Prisoner of War dynamics explains the time-dependent pattern of substitution rates in viruses


##Codes

`Evolve_sequences.nb`: Generating neutral and epistatic sites for the simulation analysis. Output: (1) initial_sample.txt (2) final_sample.txt 
`neutral.nb`: Extracting and plotting the inferred root heights and clock rates from `Evolve_sequences.nb` for the purely neutral scenario, i.e. no epistatic sites. 
`epistatic.nb`: Extracting and plotting the inferred root heights and clock rates from `Evolve_sequences.nb` for the epistatic scenario.

*stats*

Includes the simulation results for `Evolve_sequences.nb` (inferred mean clock rate), `Evolve_sequences.nb` (inferred median clock rate), and `Evolve_sequences.nb` (log file outputs) from BEAST 1.10
These files are named in the following format `*_mA_B` where letter A varies from 2 to 5 and corresponds to mutation rates ranging from 10^-2 to 10^-5, respectively and letter B

##Fasta files

`Genome_metadata.xlsx`: this contains the metadata from GISAID used for the phylogenetic analysis


