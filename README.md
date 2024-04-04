This repository hosts the code and data of theoretical modeling and computations for the manuscript "Reductive pathways in molten inorganic salts enable III-V semiconductor nanocrystal synthesis". Specifically, we include code and data for the semiempirical pseudopotential parametrization, electronic structure calculation of semiconductor nanocrystals (NC), non-adiabatic electron-phonon coupling calculations, absorption and photoluminescence (PL) spectra simulations, etc. 

Please see the manuscript and SI for details about the data. Please contact the authors (Kailai Lin, Eran Rabani) for questions regarding detailed useage of code in this repository. 

The repository and its sub-directories are structured as follows: 


## code_pseudopotFitting/

Fits semiempirical pseudopotentials to accurately reproduce literature quasi-particle band structures and deformation potentials obtained using GW approximation. 


## code_elecStruct_filterDiag/

Using the filter diagonalization technique, this code is used to obtain quasi-particle energies and eigenstates around the band edge. This code enables computationally tractable calculations for extremely large systems (such as semiconductor nanocrystals of 3000+ atoms and obtain useful optoelectronic data. 


## code_elecStruct_BSE/

Using the Bethe-Salpeter Equation (BSE), this code takes in the eigenstates and energies from the filter diagonalization code, and accounts for electron-hole interactions to obtain accurate energies for exciton states (charge neutral excitations). 


## code_elecPhCoupling/

This code calculates the non-adiabatic coupling between nuclei and electrons. The accompanying python script `compute_vklq.py` converts the electron- and hole- channel coupling with the atomic coordinate motion to the exciton-phonon coupling. 


## code_absEm_spectra/

The provided python script calculates the absorption and emission spectra using the energies and oscillator strengths of individual excitonic transitions. 


## data_pseudopot/

We provide the raw data for the atomistic III-V semiempirical pseudopotentials for In, Ga, As, P atoms (see SI for more details of the pseudopotential function form). We also provide the bulk zbIII-V unit cell configurations, the k-points data, the fitted band structures, and the fitted deformation potentials for these semiconductor systems.


## data_elecStruc/

We provide the data for the electronic structure calculation of spherical GaAs, tetrahedral GaAs, InP, CdSe, and GaP nanocrystals of a wide range of sizes. In each subdirectory of the corresponding nanocrystal calculation, we give the NC geometry (in file conf.par, note the units of atomic coordinates are Bohr radius), the quasi-particle energies (in file BSEeval.par, note the units are Hartree), the exciton energies (in file exciton.dat, with units of Hartree), and the oscillator strength of each excitonic transition (in file OS.dat). 


## data_reorgE_StokesShift/ 

The calculated reorganization energies and Stokes shifts for CdSe, GaAs and InP nanocrystals of various sizes are given. 


## data_calculated_spectra/

The calculated absorption and emission spectra, taking into account the oscillator strength of each exciton transition, as well as their phonon reorganization energies. See the SI for more details. 