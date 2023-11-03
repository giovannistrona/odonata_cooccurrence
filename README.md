# odonata_cooccurrence
The repository contains all the data and scripts needed to reproduce the analyses in:
Cerini et al. 2022 "Functional traits predict species co-occurrence patterns in a North American Odonata metacommunity".

First run the hmsc_cooc_vermont.R in R (this will take up severl hours to complete), which takes matrix.csv and sites.csv as input. 
The script will generate various output, including a supported_co-occurrence.csv file which, together with traits.csv serves as input
for the script_post_hmsc_2023.py, which will prepare the data for the Random Forest analyses. This will be conducted with the R script
RF_and_plots_2023.R, which will also generate the plots shown in the article.

The repository includes the intermediate result files, so that the individual scripts can be run independently. 
