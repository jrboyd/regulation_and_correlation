# regulation_and_correlation
The purpose of this project is to relate long non-coding (lnc) RNAs to mRNA in two different ways.
1) run_clustering_by_correlation.R - calculates a correlation coefficient for each lnc to each mRNA and then performs kmeans clustering in both dimensions.  
    Clusters with the highest correlation are then reported.
2) run_correlate_cis_regulatory.R - looks for differential expressed (DE) mRNA near DE lnc.  
    The default definition of "near" is within 100kb on either strand of DNA.  

Control of both scripts is accomplished by editing a configuration file and passing the filename to parse_config()
