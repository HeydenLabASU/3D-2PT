#!/bin/bash

if [ ! -d snapshots ]; then
echo "-submit this script from within a run-NPT+posres_x directory"
echo "-exiting"
exit
fi

#first, let's check that all files are in the right place
files=(
../posre_Protein_chain_A.itp
../posre_Protein_chain_B.itp
../complex_Protein_chain_A.itp
../complex_Protein_chain_B.itp
../complex.top
../run-NVE+posres.mdp
../../3D-2PT-files/step-3x_3D-2PT.sh
)
for file in ${files[@]}
do
if [ ! -f ${file} ]; then
echo "-could not find file ${file} relative to current directory"
echo "-you are either starting this in the wrong directory"
echo " or you missed a previous step"
echo "-exiting"
exit
fi
done

#now we count how many snapshots have been extracted from the trajectory
i=0
while [ -f snapshots/state_${i}.gro ]
do
((i+= 1))
done
if [ ${i} -eq 0 ]; then
echo "-no snapshots found"
echo "-exiting"
exit
fi
iMax=`expr $i - 1`
echo "-found ${iMax} snaphots (not counting state_0.gro, which will not be used)"

if [ -z $OMP_NUM_THREADS ]; then
nThreads=$OMP_NUM_THREADS
else
nThreads=4
fi
echo "starting ${iMax} NVE simulations with ${nThreads} threads in directories: run-NVE+posres_1-${iMax}"
echo "(you will want to parallelize this further for real applications)"
i=1
while [ $i -le ${iMax} ]
do
mkdir run-NVE+posres_${i}
cd run-NVE+posres_${i}
cp ../snapshots/state_${i}.gro start.gro
../../../3D-2PT-files/step-3x_3D-2PT.sh ${nThreads}
cd ..
((i+= 1))
done


