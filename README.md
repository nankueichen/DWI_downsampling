# A software procedure for guiding and evaluating diffusion MRI down-sampling schemes

* Here we report 1) a post-processing software procedure for down-sampling diffusion MRI for harmonizing data obtained with two different protocols (with the same b value but different numbers of diffusion gradient directions), and 2) a spatial uniformity index of diffusion encoding directions, aiming to predict the quality of the chosen down-sampling schemes. 

* The details of the methods are discussed in our journal paper authored by Nan-kuei Chen, Ryan P. Bell, and Christina S. Meade, , under revision with "Magnetic Resonance Imaging".

* The matlab code "main.m" under the "figure 1" directory was used to produce Figure 1 of our manuscript

* The matlab code "main.m" under the "bval_bvec_from_University_of_Cardiff" directory was used to process MRI data provided by University of Cardiff. The imaging data are not provided here in github. If you are interested in the imaging data, please contact the authors of https://doi.org/10.1016/j.neuroimage.2019.01.077

* The matlab codes for producing supplemental figure are in the "supplemental_figure" folder

* If you are interested in using the codes for your own data, you are encouraged to convert your bvec and bval files to the same format as the ones in "bval_bvec_from_University_of_Cardiff" directory and use the corresponding "matlab.m" file
