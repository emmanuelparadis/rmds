#Reduced Multidimensional Scaling

The reduced multidimensional scaling (RMDS) is a method to analyse very large data sets that could be not be handled by standard multidimensional scaling (MDS). The method is described in a paper in revision.

This repository contains the following files:
- `common_files.R` with two functions used in the other scripts;
- `RMDS_Euclidean.R` with code to run the MDS with Euclidean distances (typically continuous data);
- `RMDS_DNA.R` with functions to run the MDS with DNA sequence data;
- `rmds.c` C code to be used with the above scripts;
- `code.R` general-purpose functions as described in the paper;
- `application_DNA_HIV_gag.R` a script with an application using ~10k sequences from HIV.


