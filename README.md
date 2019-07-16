# CGautoDipole
Package to calculate and implement coarse-grained (CG) bead charges and Drude oscillator dipole particles from atomistic (AA) structure data files (.gro/.pdb).

The FileReaders.py script contains the class objects and methods used for AA data parsing and CG data manipulation. 

So far we have implemented methods to calculate and modify the static bead charges in a CG topology (.itp/.top) file according to AA data (.gro/.pdb). Furthermore we have implemented a method to calculate the proposed Drude Oscillator dipole moments for each of the CG beads, however implementing these oscillator particles in the CG model requires a considerable amount of further parameterisation and GROMACS file parsing for which the code still needs to be written.


