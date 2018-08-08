# ARES

This repository contains two codes, with associated jupyter notebooks as tutorials:

+ DL: dictionary learning package, written by Florent Sureau (CEA, France) and based on a previous implementation by Simon Beckouche (CEA, France). This C++ code should first be compiled with cmake. The package is self-contained, but it is necessary for performance to install atlas (<http://math-atlas.sourceforge.net>) and/or armadillo (<http://arma.sourceforge.net/docs.html>) with dependencies. Adaptation for data living on the sphere (i.e. patch extraction) is coded in python and a notebook illustrate how to use the code for denoising
+ Alpha: alpha-shearlets, including decomposition for data living on the sphere ; the code was mainly written by Felix Voigtlaender (TUB, Germany), the patchwork decomposition by Malte Wust (TUB, Germany). The original alpha shearlet package for euclidean data is described in <https://github.com/dedale-fet/alpha-transform>.

This work was funded by the DEDALE project, contract no. 665044, within the H2020 Framework Program of the European Commission.
