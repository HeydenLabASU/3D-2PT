#!/bin/bash

instdir=$1
scriptdir=$2

scripts=( 
step-1_prepare 
step-4_average-level-1 
step-5_average-level-2 
step-6_process 
3D-2PT-files/step-3x_3D-2PT 
3D-2PT-files/getS 
)

for sc in ${scripts[@]}
do
awk -v instdir=${instdir} '/^BIN/ {printf("BIN=%s/bin\n",instdir);}; ! /^BIN/ {printf("%s\n",$0);};' ${scriptdir}/${sc}.sh > ${scriptdir}/${sc}-update.sh
mv ${scriptdir}/${sc}-update.sh ${scriptdir}/${sc}.sh
chmod +x ${scriptdir}/${sc}.sh
done

if [ ! -x ${scriptdir}/3D-2PT-files/step-3x_3D-2PT.sh ]; then
chmod +x ${scriptdir}/3D-2PT-files/step-3x_3D-2PT.sh
fi

if [ ! -x ${scriptdir}/3D-2PT-files/getS.sh ]; then
chmod +x ${scriptdir}/3D-2PT-files/getS.sh
fi
