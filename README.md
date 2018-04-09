# ECLIPSE
ECLIPSE is a novel multivariate method to identify and characterize aberrant cells present in individuals out of homeostasis.

## What is ECLIPSE
The ECLIPSE package provides a novel multivariate method to identify and characterize aberrant cells present in individuals out of homeostasis. ECLIPSE combines dimensionality reduction by Simultaneous Component Analysis with Kernel Density Estimates to eliminate cells in responder samples that overlap in marker expression with cells of controls. Subsequent data analyses focus on the immune response-specific cells, leading to more informative and focused models. 
The package includes scripts and functions written for MATLAB, which is needed to apply the method. We suggest running the code in a Windows environment, as some functions could not optimally work otherwise.  ECLIPSE has been thoroughly tested using Windows 10, on MATLAB 2016b.

## Instructions
Instructions on how to use ECLIPSE can be found in the "ECLIPSE tutorial.pdf" file included in the repository.

## Authors
ECLIPSE is developed and maintained by Rita Folcarelli at the Radboud University.
The data was gathered and prepared by Bart Hilvering of the UMC Utrecht.
The GitHub repository is maintained by Roel Bouman at the Radboud University.

## References
ECLIPSE uses a modified version of the two-dimensional kernel density method developed and originally implemented by Botev et al. This method is originally described in ["Kernel Density Estimation via Diffusion"](https://arxiv.org/pdf/1011.2602.pdf).
