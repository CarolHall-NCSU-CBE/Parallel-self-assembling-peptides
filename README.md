# Parallel-self-assembling-peptides
This repository contains the data and python codes associated with the paper "Design of parallel beta-sheet nanofibrils using Monte-Carlo search, coarse-grained simulations and experimental testing" (Sarma, Sudarshan et al.), submitted to Protein Science.

# Requirement and Installation
The analysis codes are in python that can be run using an Anaconda environment. The coordinate files of the peptides can be visualized using VMD (Visual Molecular Dynamics).

# Data for paper

/DMD-PRIME20 simulations/-This directory contains the (1) coordinate files from the last simulation step and (2) the beta-sheet content for all time steps of peptides PP1-PP8 from DMD/PRIME20 simulations.

/coordinate files/-This directory contains the following files:
1. 2omm.pdb -- PDB file of GNNQQNY downloaded available in the Protein Data Bank.
2. Conf-1.pdb, Conf-2.pdb and Conf-3.pdb -- Input PDB files for PepAD in rounds 1 and 2.
3. GNNQQNY-16_peptides_2omm -- This is the PDB file of the biological assembly of 16 GNNQQNY peptides (Class 1 cross beta-spine). Crystal tranformation information is available in PDB ID: 2omm.
4. PP1-PP8.pdb -- Output from PepAD for peptides PP1-PP8.

 /analysis codes/-This directory contains the following python scripts:
 1. eisenberg_transformation.py -- Python script to be able to create a peptide assembly from the unit peptide in PDB ID: 2omm.
 2. beta_sheet_plot_cont.py -- Python script to plot beta-sheet content and compute average beta-sheet content from data available in the /DMD-PRIME20 simulations/ directory.
 3. parallel-antiparallel.py -- Python script to compute total parallel beta-sheet peptides and antiparallel beta-sheet peptides from DMD/PRIME20 simulations data. The input files required for the calculate are the ".bptnr" files that can be made available on request.
 4. fortran.py -- This file has to be loaded along with the parallel-antiparallel.py python script. 
