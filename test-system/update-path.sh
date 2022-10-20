#!/bin/bash

instdir=$1

scripts=( 
step-1_prepare 
step-4_average-level-1 
step-5_average-level-2 
step-6_process 
3D-2PT-files/step-3x_3D-2PT 
3D-2PT/getS 
)

for sc in ${scripts[@]}
do
awk -v instdir=${instdir} '/^BIN/ {printf("BIN=%s/bin\n",instdir);}; ! /^BIN/ {printf("%s\n",$0);};' ${sc}.sh > ${sc}-update.sh
mv ${sc}-update.sh ${sc}.sh
done

if [ ! -x 3D-2PT-files/step-3x_3D-2PT.sh ]; then
chmod +x 3D-2PT-files/step-3x_3D-2PT.sh
fi

if [ ! -x 3D-2PT/getS.sh ]; then
chmod +x 3D-2PT/getS.sh
fi
