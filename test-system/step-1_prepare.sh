#!/bin/bash

scriptdir=`awk '{split($0,a,"/");printf("%s\n",a[1]);};' << STOP
$0
STOP`

BIN=bin

if [ ! -f ${BIN}/water3D_noRot.exe ]; then
echo "provide installation directory for 3D-2PT:"
read instdir
BIN=${instdir}/bin
if [ ! -x update-path.sh ]; then
chmod update-path.sh
fi
${scriptdir}/update-path.sh ${instdir}
fi

while [ ! -f ${BIN}/water3D_noRot.exe ]
do
echo "did not find 3D-2PT executables in directory: ${BIN}"
echo "provide correct path:"
read instdir
BIN=${instdir}/bin
if [ ! -x update-path.sh ]; then
chmod update-path.sh
fi
${scriptdir}/update-path.sh ${instdir}
done

#first, let's check that all files are in the right place
files=(
posre_Protein_chain_A.itp
posre_Protein_chain_B.itp
complex_Protein_chain_A.itp
complex_Protein_chain_B.itp
complex.pdb
complex.top
em+posres.mdp
equi-NPT+posres.mdp
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

gmx editconf -f complex.pdb -box 3 3 3 -o box.gro >& editconf.out

gmx solvate -cp box.gro -cs -p complex.top -o solv.gro >& solvate.out

mkdir em+posres
cd em+posres
gmx grompp -f ../em+posres.mdp -c ../solv.gro -r ../solv.gro -p ../complex.top -o topol.tpr -maxwarn 1 >& grompp.out
gmx mdrun -v -nt 1 -s topol.tpr -o traj.trr -e ener.edr -g md.log -c confout.gro -cpo state.cpt >& mdrun.out
cd ..

mkdir equi-NPT+posres
cd equi-NPT+posres
gmx grompp -f ../equi-NPT+posres.mdp -c ../em+posres/confout.gro -p ../complex.top -r ../solv.gro -o topol.tpr -maxwarn 1 >& grompp.out
gmx mdrun -v -nt 16 -s topol.tpr -o traj.trr -e ener.edr -g md.log -c confout.gro -cpo state.cpt >& mdrun.out
cd ..

${BIN}/gmxtop-2.exe << STOP >& convertGMXTOP.out
complex.top
complex_Protein_chain_A.itp
complex_Protein_chain_B.itp
/scratch/mheyden1/ResearchOnline/amber99sb-ildn.ff/tip3p.itp
/scratch/mheyden1/ResearchOnline/amber99sb-ildn.ff/ffnonbonded.itp
complex.mtop
STOP

cat << STOP >& tmp.tcl
mol new box.gro
set sel1 [atomselect top "not water and not hydrogen"]
set sel2 [atomselect top all]
set com [measure center \$sel1 weight mass]
set sh [vecscale -1 \$com]
\$sel2 moveby \$sh
animate write pdb align_ref.pdb sel \$sel1
exit
STOP
vmd -e tmp.tcl -dispdev none
rm tmp.tcl
gmx trjconv -s align_ref.pdb -f align_ref.pdb -o align_ref.gro << STOP >& trjconv.out
0
STOP
