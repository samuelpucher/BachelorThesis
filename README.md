# BachelorThesis
This repository contains all the raw code, that I used for my Bachelor thesis. A clone of this repository exists on the GitLab server of the institute for theoretical physics at JKU

## Author:
Samuel Pucher

## Created:
28th October 2024

## Overview:
The repositroy contains all the raw files (.cpp), that can be used to derive all the results from my Bachelor thesis.
In the folder './CorrelatedPolarons', the program is storred, which claculated all my results for the correlated polaron equations.
The mean field approximations can be found in the folder './MeanField' and several important numerical tests are located in './ImportantTests'.
In the folder './CorrelatedPolarons(3D)', the "experimental" implementation of the correlated polaron eqauitons, which is using isotropy, is stored. This implementation did not work unfortunately.

## Miscellaneous:
Some of the programs need libraries linked to compile, for example FFTW. Make sure the libraries are installed on the machine and exessible to the linker of the compiler. 
For these programs, the g++ compiler was used with -O2 optimization.
