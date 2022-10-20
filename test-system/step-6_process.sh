#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 0:02:00                  # wall time (D-HH:MM)

module load gromacs/2018.1
module load vmd/1.9.3

if [ ! -d average ]; then
echo "-could not find directory average"
echo "-you are either starting this in the wrong directory"
echo " or you missed a previous step"
echo "-exiting"
exit
fi

#first, let's check that all files are in the right place
files=(
complex.pdb
../3D-2PT-files/water3D.input
)
for file in ${files[@]}
do
if [ ! -f ${file} ]; then
echo "-could not find file ${file} in current directory"
echo "-you are either starting this in the wrong directory"
echo " or you missed a previous step"
echo "-exiting"
exit
fi
done

cd average

cat << STOP >& tmp.tcl
mol new ../complex.pdb
set sel1 [atomselect top "not water and not hydrogen"]
set sel2 [atomselect top all]
set sel3 [atomselect top "not water"]
set com [measure center \$sel1 weight mass]
set sh [vecscale -1 \$com]
\$sel2 moveby \$sh
animate write pdb prot.pdb sel \$sel3
exit
STOP
vmd -e tmp.tcl -dispdev none
rm tmp.tcl
gmx trjconv -s prot.pdb -f prot.pdb -o prot.gro << STOP >& trjconv.out
0
STOP

#The prefix is part of the file name for all output files and was set in the 3D-2PT input files
prefix=`tail -n 1 ../../3D-2PT-files/water3D.input`

files=(
${prefix}_pot3D-numDens.cube
${prefix}_pot3D-numDens.cube
${prefix}_pot3D-Uxs.cube
${prefix}_pot3D-Uss.cube
${prefix}_water3D-numDens.cube
${prefix}_water3D-totalS.cube
)
for file in ${files[@]}
do
if [ ! -f ${file} ]; then
echo "-could not find file average/${file} in current directory"
echo "-you are either starting this in the wrong directory"
echo " or you missed a previous step"
echo "-exiting"
exit
fi
done

echo "-processing 3D-2PT data in directory average"
cat << STOP >& process.input
prot.gro
${prefix}_pot3D-numDens.cube
10.0
${prefix}_pot3D-numDens.cube
${prefix}_pot3D-Uxs.cube
${prefix}_pot3D-Uss.cube
${prefix}_water3D-numDens.cube
${prefix}_water3D-totalS.cube
300
0.2
STOP
/scratch/mheyden1/MadR2014/grid-dynamics/process-cube.exe < process.input >& process.out
echo "-done"

cd ..

i=1
flag=0
while [ -d run-NPT+posres_${i} ]
do
    cd run-NPT+posres_${i}
    if [ ! -d average ]; then
        echo "-could not find directory run-NPT+posres_${i}/average"
    else
        cd average
        for file in ${files[@]}
        do
            if [ ! -f ${file} ]; then
                echo "-could not find file run-NPT+posres_${i}/average/${file}"
                flag=1
            fi
        done
        if [ ${flag} -eq 0 ]; then
            echo "-processing 3D-2PT data in directory run-NPT+posres_${i}/average"
            cp ../../average/prot.gro .
            cp ../../average/process.input .
            /scratch/mheyden1/MadR2014/grid-dynamics/process-cube.exe < process.input >& process.out
            echo "-done"
        else
            echo "-skipping directory run-NPT+posres_${i}/average"
        fi
        cd ..
        flag=0
    fi
    cd ..
    ((i+= 1))
done
if [ ${i} -eq 1 ]; then
echo "-no sub-driectories run-NPT+posres_x found for processing of level 1 data"
echo "-exiting"
exit
fi
iMax=`expr $i - 1`
