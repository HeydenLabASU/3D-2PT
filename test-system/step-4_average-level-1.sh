#!/bin/bash

BIN=<your-3D-2PT-installation-directory>/bin

if [ ! -d snapshots ]; then
echo "-submit this script from within a run-NPT+posres_x directory"
echo "-exiting"
exit
fi

#first, let's check that all files are in the right place
files=(
../../3D-2PT-files/water3D.input
../../3D-2PT-files/getS.sh
)
for file in ${files[@]}
do
if [ ! -f ${file} ]; then
echo "-could not find file ${file} relative to current directory"
echo "-you are probably starting this in the wrong directory"
echo "-exiting"
exit
fi
done

cubes=( pot3D-numDens water3D-COMmsd-1.600ps water3D-plane2Prot water3D-rot1-1.600ps pot3D-Uss water3D-numDens water3D-pol2Prot water3D-rot2-1.600ps pot3D-Uxs water3D-OH2Prot water3D-pol water3D-rot3-1.600ps )

avCubes=( pot3D-numDens water3D-numDens )

w1avCubes=( pot3D-Uss pot3D-Uxs )

w2avCubes=( water3D-COMmsd-1.600ps water3D-plane2Prot water3D-rot1-1.600ps water3D-pol2Prot water3D-rot2-1.600ps water3D-OH2Prot water3D-pol water3D-rot3-1.600ps )

data=( water3D-AngMomProj1-VDOS_df-10.450wn_maxf-2090.000wn water3D-HydVel-VDOS_df-10.450wn_maxf-2090.000wn water3D-AngMomProj2-VDOS_df-10.450wn_maxf-2090.000wn water3D-AngMomProj3-VDOS_df-10.450wn_maxf-2090.000wn water3D-OxyVel-VDOS_df-10.450wn_maxf-2090.000wn water3D-COMvel-VDOS_df-10.450wn_maxf-2090.000wn )

Scubes=( water3D-rotS water3D-totalS water3D-transS )

#The prefix is part of the file name for all output files and was set in the 3D-2PT input files
prefix=`tail -n 1 ../../3D-2PT-files/water3D.input`

#now we count all run-NVE+posres_x directories and check their content for completeness
i=1
while [ -d run-NVE+posres_${i} ]
do
cd run-NVE+posres_${i}
    for file in ${cubes[@]}
    do
        if [ ! -f ${prefix}_${file}.cube ]; then
            echo "-file run-NVE+posres_${i}/${prefix}_${file}.cube not found"
            echo "-step 3 has not been completed yet or caused an error"
            echo "-exiting"
            exit
        fi
    done
    for file in ${data[@]}
    do
        if [ ! -f ${prefix}_${file}.dat ]; then
            echo "-file run-NVE+posres_${i}/${prefix}_${file}.dat not found"
            echo "-step 3 has not been completed yet or caused an error"
            echo "-exiting"
            exit
        fi
    done
    echo "-directory run-NVE+posres_${i} contains complete 3D-2PT data set"
cd ..
((i+= 1))
done
if [ ${i} -eq 1 ]; then
echo "-no 3D-2PT data found for averaging"
echo "-exiting"
exit
fi
iMax=`expr $i - 1`
echo "-found ${iMax} directories with complete 3D-2PT data sets for averaging"

#We extract the dimensions of the analysis grid directly from one of the output files
grid=`head -n 6 run-NVE+posres_1/${prefix}_water3D-numDens.cube | tail -n 3 | awk '{printf("%d ",$1);};'`

mkdir average
cd average

cubeIdx=0
while [ ${cubeIdx} -lt ${#avCubes[@]} ]
do
i=1
while [ $i -le ${iMax} ]
do
if [ $i -eq 1 ]; then
echo "${iMax}"
fi
echo "../run-NVE+posres_${i}/${prefix}_${avCubes[cubeIdx]}.cube"
if [ $i -eq ${iMax} ]; then
echo "${prefix}_${avCubes[cubeIdx]}.cube"
fi
((i+= 1))
done >& avercube.input
${BIN}/averCubefiles.exe < avercube.input >& avercube.out
((cubeIdx+= 1))
done

cubeIdx=0
while [ ${cubeIdx} -lt ${#w1avCubes[@]} ]
do
i=1
while [ $i -le ${iMax} ]
do
if [ $i -eq 1 ]; then
echo "${iMax}"
fi
echo "../run-NVE+posres_${i}/${prefix}_${w1avCubes[cubeIdx]}.cube"
echo "../run-NVE+posres_${i}/${prefix}_pot3D-numDens.cube"
if [ $i -eq ${iMax} ]; then
echo "${prefix}_${w1avCubes[cubeIdx]}.cube"
fi
((i+= 1))
done >& avercube.input
${BIN}/weightedAverCubeFiles.exe < avercube.input >& avercube.out
((cubeIdx+= 1))
done

cubeIdx=0
while [ ${cubeIdx} -lt ${#w2avCubes[@]} ]
do
i=1
while [ $i -le ${iMax} ]
do
if [ $i -eq 1 ]; then
echo "${iMax}"
fi
echo "../run-NVE+posres_${i}/${prefix}_${w2avCubes[cubeIdx]}.cube"
echo "../run-NVE+posres_${i}/${prefix}_water3D-numDens.cube"
if [ $i -eq ${iMax} ]; then
echo "${prefix}_${w2avCubes[cubeIdx]}.cube"
fi
((i+= 1))
done >& avercube.input
${BIN}/weightedAverCubeFiles.exe < avercube.input >& avercube.out
((cubeIdx+= 1))
done

for file in ${cubes[@]}
do
    if [ ! -f ${prefix}_${file}.cube ]; then
        echo "-file average/${prefix}_${file}.cube not found"
        echo "-averaging seems to have failed"
        echo "-exiting"
        exit
    fi
    size=`ls -rtl ${prefix}_${file}.cube | awk '{printf("%s\n",$5);};'`
    if [ ${size} -lt 100 ]; then
        echo "-file average/${prefix}_${file}.cube is too small (${size} kB)"
        echo "-averaging seems to have failed"
        echo "-exiting"
        exit
    fi
done

dataIdx=0
while [ ${dataIdx} -lt ${#data[@]} ]
do
i=1
while [ $i -le ${iMax} ]
do
if [ $i -eq 1 ]; then
echo "${iMax}"
fi
echo "../run-NVE+posres_${i}/${prefix}_${data[dataIdx]}.dat"
if [ $i -eq ${iMax} ]; then
echo "200"
echo "${grid}"
echo "${prefix}_${data[dataIdx]}.dat"
fi
((i+= 1))
done >& average.input
${BIN}/average.exe average.input >& average.out
((dataIdx+= 1))
done

for file in ${data[@]}
do
    if [ ! -f ${prefix}_${file}.dat ]; then
        echo "-file average/${prefix}_${file}.dat not found"
        echo "-averaging seems to have failed"
        echo "-exiting"
        exit
    fi
    size=`ls -rtl ${prefix}_${file}.dat | awk '{printf("%s\n",$5);};'`
    if [ ${size} -lt 100 ]; then
        echo "-file average/${prefix}_${file}.dat is too small (${size} kB)"
        echo "-averaging seems to have failed"
        echo "-exiting"
        exit
    fi
done

../../../3D-2PT-files/getS.sh ${prefix} >& getS.out

for file in ${Scubes[@]}
do
    if [ ! -f ${prefix}_${file}.cube ]; then
        echo "-file average/${prefix}_${file}.cube not found"
        echo "-entropy calculation seems to have failed"
        echo "-exiting"
        exit
    fi
    size=`ls -rtl ${prefix}_${file}.cube | awk '{printf("%s\n",$5);};'`
    if [ ${size} -lt 100 ]; then
        echo "-file average/${prefix}_${file}.cube is too small (${size} kB)"
        echo "-entropy calculation seems to have failed"
        echo "-exiting"
        exit
    fi
done

cd ..

echo "-averaging 3D-2PT output for ${iMax} data sets complete"
echo "-deleting single data sets to free disk space"
i=1
while [ $i -le ${iMax} ]
do
echo "-deleting directory run-NVE+posres_${i}"
rm -r run-NVE+posres_${i}
((i+= 1))
done
echo "DONE"
