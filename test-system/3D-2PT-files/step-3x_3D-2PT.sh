#!/bin/bash

BIN=../../../../bin

#first, let's check that all files are in the right place
files=(
../../posre_Protein_chain_A.itp
../../posre_Protein_chain_B.itp
../../complex_Protein_chain_A.itp
../../complex_Protein_chain_B.itp
../../complex.top
../../complex.mtop
../../align_ref.gro
../../solv.gro
../../run-NVE+posres.mdp
../../../3D-2PT-files/groups.job
../../../3D-2PT-files/water3D.input
../../../3D-2PT-files/pot3D.input
start.gro
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

gmx grompp -f ../../run-NVE+posres.mdp -c start.gro -p ../../complex.top -r ../../solv.gro -o topol.tpr -maxwarn 1 >& grompp.out
gmx mdrun -v -nt 4 -s topol.tpr -o traj.trr -e ener.edr -g md.log -c confout.gro -cpo state.cpt >& mdrun.out
gmx trjconv -s topol.tpr -f traj.trr -pbc mol -o traj_pbc.trr << STOP >& trjconv.out
0
STOP
rm traj.trr

#openMP parallelization available, but test for each system recommended
#parallelization works best for pot3D
export OMP_NUM_THREADS=1
${BIN}/water3D_noRot.exe ../../../3D-2PT-files/water3D.input >& water3D.out
${BIN}/pot3D.exe ../../../3D-2PT-files/pot3D.input >& pot3D.out

rm traj_pbc.trr
