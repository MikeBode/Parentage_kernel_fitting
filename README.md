# Parentage_kernel_fitting
Methods for fitting dispersal kernels to genetic parentage data

======= DESCRIPTION =======

This project contains (1) MATLAB script files that fit dispersal kernels to genetic parentage assignment datasets, gathered from patchy habitat; and (2) example comma separated value (CSV) data files that contain the necessary information in the correct format.

The article describing this method is:
Bode M, Williamson D, Outram N, Harrison H, Jones G. Estimating dispersal kernels using genetic parentage data. Submitted to Methods in Ecology and Evolution on November 1, 2016.

The files comprise CSV input files:

1.	Example_parentage_matrix.csv
2.	Example_proportion_matrix.csv
3.	Example_reef_size_matrix.csv
4.	Example_sampled_reef_list.csv
5.	Example_reef_location.csv
6.	Example_distance_matrix.csv

And MATLAB scripts:

7.	Fit_Example_Data.m
8.	Kernel_Fitting_Function.m

The structure and content of the CSV input files is as follows:

Example_parentage_matrix.csv
This is a (P+1) x P matrix, where P is the number of sampled patches. Each element contains the number of sampled juveniles from the COLUMN patch whose parents were also sampled on the ROW patch.

Example_proportion_matrix.csv
This is a P x 1 matrix, where each element indicates the estimated proportion of adults on the ROW patch that were sampled and genotyped.

Example_patch_size_matrix.csv
This is an N x 1 matrix, where N is the total number of patches in the metapopulation. Each element indicates the size of the patch (units are unimportant since only relative size matters).

Example_sampled_reef_list.csv
This is a 1 x P matrix, where each element identifies a sampled patch from the N total patches. This data structure links the total number of patches (N) with the sampled patches (P). 

Example_reef_location.csv
This is a N x 2 matrix, where the first column gives the longitude, and the second column gives the latitude, of each patch in the system. This function is not directly used by the project, but (a) it is the basis for the distance matrix, which is essential; and (b) it can be used to visualise the metapopulation and the sampled patches for this example.

Example_distance_matrix.csv
This N x N matrix gives the pairwise Euclidean distance between each patch in the metapopulation.

The purpose of the MATLAB script files is as follows:

Fit_Example_Data.m
This is the central function. It pulls in the input data from the CSV files above. It then defines the candidate kernel shapes as inline Matlab functions. It uses an inbuilt single dimensional optimisation MATLAB routine called FMINBND to search for the best-fit scale parameter for a particular dispersal kernel shape. The performance of a given shape and scale is determined by the function Kernel_Fitting_Function.m, defined below. The script then repeats this search for a number of candidate dispersal kernel shapes. The best fit kernel is the shape and scale that minimises the log likelihood of the parentage assignment dataset. Finally, the method repeatedly subsamples (with replacement) the parentage assignment matrix at the patch scale, and repeats the scale parameter fitting procedure, using the previously identified best-fit dispersal kernel shape. These bootstrap fits are used to construct confidence bounds around the dispersal kernel scale parameter and the resulting estimate of the mean dispersal distance.

Relevant outputs are displayed in the MATLAB command window. All outputs are saved in the file OUTPUTS.MAT.

Kernel_Fitting_Function.m
This is the function that calculates the log likelihood of a particular parentage assignment dataset, given a dispersal kernel shape and scale parameter. 

======= REQUIREMENTS =======
MATLAB, preferably version 9.2.0.556344 (R2017a) or later. 

Refer any questions to:
Email: michael.bode@jcu.edu.au

