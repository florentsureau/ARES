Dictionary Learning Package

Florent Sureau, based on a previous implementation by Simon Beckouche.
florent.sureau@cea.fr
Service d'Astrophysique, CEA Saclay, August 2018


The packages contains the following executables that can be used for dictionary
learning:
+ DL_learning: perform dictionary learning. Sparse Coding is currently only
performed by OMP. The dictionary update can be obtained via MOD, K-SVD, approximate K-SVD
(alternate minimization) (see Rubinstein and Elad papers).
+ sparse_decomp: sparse coding using OMP.

The package should be compiled with cmake. The package is self-contained but it is strongly recommended to install atlas (<http://math-atlas.sourceforge.net>)
and/or armadillo libraries (<http://arma.sourceforge.net/docs.html>).

0) Quick use 
If you already know or do not care about the details, here is a list of commands
to use to run the functions of this package. The image in these examples
are in the same directory as the executables.

+ bin/dl_learning -g 512 -R -n 128 -i100 -s 5. -V -K -k 1 -a approx_approxkSVD1_8x8_5000_it100.fits
patches8x8_5000.fits dicolearn_approxkSVD1_8x8_5000_it100.fits : use approximate KSVD (-K -k 1) with one iteration to
learn a dictionary of 128 atoms (-n 128) using 100 iterations (-i 100), based on fits training set patches8x8_5000.fits. Dictionary
is initialized by taking random examples (-R, with seed 512 : -g 512). For sparse coding, the targeted sparsity is
5 (-s 5), verbosity is on (-V) and we save the learn dictionary in dicolearn_approxkSVD1_8x8_5000_it100.fits.
The approximation obtained is saved in approx_approxkSVD1_8x8_5000_it100.fits

OR

+ bin/dl_learning -g 512  -i100 -s 5. -V -d patches8x8_128.fits patches8x8_5000.fits dicolearn_MOD_8x8_5000_it100.fits : use
MOD to learn a dictionary using 100 iterations (-i 100), based on fits training set patches8x8_5000.fits. Dictionary
is initialized by reading patches8x8_128.fits (128 atoms). For sparse coding, the targeted sparsity is
5 (-s 5), verbosity is on (-V) and we save the learn dictionary in dicolearn_MOD_8x8_5000_it100.fits .


+ bin/sparse_decomp -s 5 dicolearn_MOD_8x8_5000_it100.fits patches8x8_5000.fits code_MOD_8x8_5000.fits : sparse decomposition with sparsity 5 of patches in learned dictionary.


1) dl_learning

This binary learns a sparsifying dictionary from patch data using alternated minimization (-i nbIt, where nbIt is the total
number of iterations). OMP is used as sparse coder by either specifying the target sparsity (-s), and/or the expected RMS level (-e), and either MOD (default), KSVD (-K) or
approximate KSVD (-k x, wiht x the number of sub-iterations) for dictionary update. A initial dictionary can also be specified (-d NameofInitDico.fits) or it can either be a sparse dictionary with targeted density (-r) or taken from random samples in the training set (-R) with a given number of atoms (-n). Input metric in sparse coding (-m) can be provided as full 2D array or only its diagonal.

Usage: bin/DL_learning options nameVectIn.fits nameDicoOut.fits 

   where options =  
      [-e ErrorTarget]
      [-s TargetSparsity]
      [-i Iteration number, default 10]
      [-d InitialDictionary]
      [-m name of fits file containing input metric]
      [-n Natoms if no input dictionary]
      [-p Name of mean patch to save and to subtract
          from training set]
      [-r Init Sparse Density if no input dictionary]
      [-R Random training examples for initialization]
      [-g generator seed - default 0]
      [-c Name for training set code to save]
      [-a Name for training set approx to save]
      [-C center training patches (mean of each patch=0)]
      [-K use KSVD]
      [-k Iterations for approx KSVD (used), default 1]
      [-M Minimal Correlation]
      [-N Minimal Sparsity]
      [-T Timing on]
      [-V verbose]

2) sparse_decomp

This binary performs sparse decomposition (only OMP implemented yet).The target sparsity (-s) and/or the expected RMS level (-e) can be specified. Input metric, used in the dot product, can be specified using -m (either a diagonal, or full 2D metric). If a mean patch is assumed, it can be provided with option (-p), and for centered (mean of the patch=0) patch, one can use -C. Approximation can be saved using -a, with either the mean of the patch added back if -C, or the mean patch added if -p. If the codes are of no interest, the option -N can be used.

Usage: bin/sparse_decomp options nameDicoIn.fits nameVectIn.fits nameCodeOut.fits 

   where options =  
      [-e ErrorTarget, TargetSparsity ignored]
      [-s TargetSparsity]
      [-m name of fits file containing input metric for selection of atoms]
      [-w name of fits file containing input metric for selection and approximation]
      [-p Name of optional mean input to subtract
          from each input]
      [-a Name for training set approx to save]
      [-C center training patches (mean of each patch=0)]
      [-N do not save the final code]
      [-T Timing on]
      [-V verbose]



3) Other executables:

+ DL_test_matrix_decomposition : perform a series of test on matrix
decompositions, with execution time if needed.
