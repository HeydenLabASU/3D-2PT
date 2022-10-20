#!/bin/bash

if [ -z $1 ]; then
echo "provide prefix for 3D-2PT files, e.g. grid-96x96x96"
exit
fi

BIN=<your-3D-2PT-installation-directory>/bin

prefix=$1

${BIN}/trans-vdos2entropy_lowMem.exe ${prefix}_water3D-COMvel-VDOS_df-10.450wn_maxf-2090.000wn.dat ${prefix}_water3D-numDens.cube 33.4273 > ${prefix}_water3D-transS.cube
${BIN}/rot-vdos2entropy_lowMem.exe ${prefix}_water3D-AngMomProj1-VDOS_df-10.450wn_maxf-2090.000wn.dat ${prefix}_water3D-AngMomProj2-VDOS_df-10.450wn_maxf-2090.000wn.dat ${prefix}_water3D-AngMomProj3-VDOS_df-10.450wn_maxf-2090.000wn.dat ${prefix}_water3D-numDens.cube 33.4273 > ${prefix}_water3D-rotS.cube
cat << STOP >& tmp
2
${prefix}_water3D-transS.cube
${prefix}_water3D-rotS.cube
STOP
${BIN}/add_cubes.exe tmp > ${prefix}_water3D-totalS.cube
rm tmp
