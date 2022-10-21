# 3D-2PT
Installation instructions
requirements: installed single precision FFTW3 library (www.fftw.org)

To compile the 3D-2PT executables, download the files in this repository and type:

$ make

To get started with 3D-2PT calculations, enter the "test-system" directory
$ cd test-system
The following software is required to run the test simulations:
- gromacs (www.gromacs.org)
- VMD (Visual Molecular Dynamics, http://www.ks.uiuc.edu/Research/vmd)

The test system consists of a sodium and a chloride ion separated by 3.0A in water.
Enter the corresponding sub-directory:

$ cd NaCl-3.00A

Execute the automated scripts in the sequence described in the flowchart (flowchart.pdf).
You can repeat steps-2,3,4 an arbitrary number of times to improve statistics. 
2 iterations are sufficient as a proof of concept that thesoftware works, but the noise-level will not allow for any meaningful interpretations.
50 iterations (or changes in the protocol that provide equivalent statistics) are needed for quantitative results.
Step-5 averages all data generated in step-2,3,4 iterations and step-6 is used for final processing. 
The final results are then stored in the NaCl-3.00A/average folder.

The main results are stored in CUBE files (http://paulbourke.net/dataformats/cube/), an ASCII file format traditionally used to store electron densities from electronic structure calculations (quantum chemistry). However, we use this format to store multiple types of volumetric data as described in the header of each CUBE output file, e.g., water number densities, solvation free energies per water molecule, etc.
Integration over the spatial volume provides the total solvation free energy, enthalpy, entropy etc., which is also done autmoatically in step-6. This data can be used as a reference to ensure that you interpret the results correctly.

Enjoy!
Feel free to send questions and comments to: mheyden1@asu.edu
