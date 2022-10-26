# 3D-2PT
This program can be used to perform a spatially resolved analysis of solvation free energy, enthalpy and entropy calculations based on all-atom molecular dnamics simulations in explicit water.

Applications include the research presented in the following peer-reviewed research articles:
- Persson, Pattni, Singh, Kast & Heyden, J. Chem. Theory Comput. 13, 4467-4481 (2017).
- Pattni, Vasilevskaya, Thiel & Heyden, J. Phys. Chem. B 121, 7431-7442 (2017).
- Heyden, WIREs Comp. Mol. Sci. e1390 (2018).
- Heyden, J. Chem. Phys. 150, 094701 (2019).
- Pattni & Heyden, J. Phys. Chem. B 123, 6014-6022 (2019).
- Conti Nibali, Pezzotti, Sebastiani, Galimberti, Schwaab, Heyden, Gaigeot & Havenith, J. Phys. Chem. Lett. 11, 4809-4816 (2020).
- Fajardo & Heyden, J. Phys. Chem. B. 125, 4634-4644 (2021).
- Waskasi, Lazaric, Heyden, Electrophoresis 42, 2060-2069 (2021).

Installation instructions:

Requirements: installed single precision FFTW3 library (www.fftw.org)

To compile the 3D-2PT executables, download the files in this repository and type:

$ make

Test Case:

To get started with 3D-2PT calculations, enter the "test-system" directory

$ cd test-system

The following software is required to run the test simulations:
- gromacs (www.gromacs.org)
- VMD (Visual Molecular Dynamics, http://www.ks.uiuc.edu/Research/vmd)

The test system consists of a sodium and a chloride ion separated by 3.0A in water.
Enter the corresponding sub-directory:

$ cd NaCl-3.00A

Execute the automated scripts (parent directory: ../step-*.sh) in the sequence described in the provided flowchart (flowchart.pdf).

You can repeat steps-2,3,4 an arbitrary number of times to improve statistics:
- 2 iterations are sufficient as a proof of concept that the software works, but the noise-level will not allow for any meaningful interpretation.
- 40-50 iterations (or modified protocols that provide equivalent statistics) are needed for quantitative results.

In step-5, all data generated in previous step-2,3,4 iterations is averaged and step-6 is used for final processing. 
Final results will be stored in the "NaCl-3.00A/average" folder.

Most results are stored in so-called CUBE files (http://paulbourke.net/dataformats/cube/), an ASCII file format that is traditionally used to store electron densities from electronic structure calculations (quantum chemistry). However, we use this format to store multiple types of volumetric data as described in the header of each CUBE output file, e.g., water number densities, solvation free energies per water molecule, etc.

Integration over the spatial volume provides the total solvation free energy, enthalpy, entropy etc., which is also done autmoatically in step-6 (see file: "NaCl-3.00A/average/results_bulk-distance-to-atom-10.0A_minDensRatio-0.20.dat" after completing all steps and iterations). This data can be used as a reference to ensure that you interpret the results correctly.

Enjoy!
Feel free to send questions and comments to: mheyden1@asu.edu
