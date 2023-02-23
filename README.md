# ArielNikas-VEHeterogeneity
This project considers the effect of competing heterogeneities on vaccine effectiveness when protection is either constant or waning. The accompanying paper has been submitted to 'Clinical Infectious Diseases'. There are two main parts to the code:

1)Julia code for running the simulation, the code gives a single example simulation for the risk-correlate model. Please note that this example simulation will creat a two files begining with "Example-RiskCorrelate" in your current working directory. To find your current working directory use the pwd() command in Julia. Included is the basic simulation used through most of the paper. Julia files end in .jl 

2)R code runs the analysis and creates the corresponding figure. This code is set up to read in two files and produce one figure that is a panel from the accompanying paper's Figure 4. 
