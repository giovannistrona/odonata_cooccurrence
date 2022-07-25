# odonata_cooccurrence
The repository contains all the data and scripts needed to reproduce the analyses in:
Cerini et al. 2022 "Functional traits predict species co-occurrence patterns in a North American Odonata metacommunity".

First run the hmsc_cooc_vermont.R in R (this will take up severl hours to complete), which takes matrix.csv and sites.csv as input. 
The script will generate various output, including a supported_co-occurrence.csv file which, together with traits.csv serves as input
for the script.py, which will perform the Random Forest analyses.
