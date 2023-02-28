# ArielNikas-VEHeterogeneity
This project considers the effect of competing heterogeneities on vaccine effectiveness when protection is either constant or waning. The accompanying paper has been submitted to 'Clinical Infectious Diseases'. There are two main parts to the code:

1)Julia code for running each simulation, named for which Figure they produce a panel for. Please note that these example simulations will creat a two files each in your current working directory. To find your current working directory use the pwd() command in Julia.  We provide an example simulation result for the risk-correlate model (Figure 4). Julia files end in .jl . Manifest and Project .toml files give exact package versions used for reproducibility. 

Note: Julia version 1.3 (though it still works through 1.7 with minor differences due to RNG)

2)R code runs the analysis and creates the corresponding figure or figure panel. There are separate labelled chunks for each figure. You need the corresponding csv files created by Julia to run them. An example set of csv files (Example-RiskCorrelated-1.csv &  Example-RiskCorrelated-2.csv) is provided. 


Note: R version 4.2.1, 'survival' package version 3.3-1
