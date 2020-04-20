# Conidial-heterogeneity

For analyzing CFW and Katushka labeling of single cells, imageJ script 'manual selection of single cells and analyse Katushka mask to CFW signal.ijm' was used on Z-stacks in time. Alternatively, 'analyse CFW mask to CFW signal.ijm' was used on the 'manual selection of single cells and analyse Katushka mask to CFW signal.ijm' processed data.

For analyzing heterogeneity within cell populations the flow cytometry data was first combined of all replicates using 'combine_3_reps_flowjo_cbind_NA.R' and then normalized using 'normalise_3_reps_by_lambda_weighted_mu.R'. These data were then feeding into the heterogeneity analysis script 'Subpopulation analysis and heterogeneity quantification.R'. 
