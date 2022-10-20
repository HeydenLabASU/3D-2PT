# 3D-2PT
Installation instructions
requirements: installed single precision FFTW3 library (www.fftw.org)

To compile the 3D-2PT executables, download the files in this repository and type:
$ mkdir bin\n
$ make\n

To get started with 3D-2PT calculations, enter the "test-system" directory
$ cd test-system
The following software is required to run the test simulations:
- gromacs (www.gromacs.org)
- VMD (Visual Molecular Dynamics, http://www.ks.uiuc.edu/Research/vmd)

The test system consist of a sodium and a chloride ion separated by 3.0A in water.
Enter the corresponding sub-directory:
adapt path to your installation directory set in variable BIN in the following scripts:
- step-1_prepare.sh
- step-4_average-level-1.sh
- step-5_average-level-2.sh
- step-6_process.sh 
- 3D-2PT-files/step-3x_3D-2PT.sh
- 3D-2PT/getS.sh

Enter the sub-directory with prepared input:
$ cd NaCl-3.00A
Execute the automated scripts in the sequence described in Flowchart.pdf
