#!/bin/bash

#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 0-24:00                  # wall time (D-HH:MM)
#SBATCH -o step-2.out

module load gromacs/2018.1

#first, let's check that all files are in the right place
files=(
posre_Protein_chain_A.itp
posre_Protein_chain_B.itp
complex_Protein_chain_A.itp
complex_Protein_chain_B.itp
complex.top
solv.gro
equi-NPT+posres/confout.gro
run-NPT+posres.mdp
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

#By default, the simulation will start from the coordinates 
#specificied in file 'equi-NPT+posres.gro'.
#However, if the output files of previous simulations are 
#detected, the simulation will start from the final set of
#coordinates generated by the last simulation.
start=equi-NPT+posres/confout.gro
i=1
#This loop will keep increasing the counter variable 'i' 
#and updat the value of the variable 'start' until 
#the file 'run-NPT+posres_${i}/confout.gro' does NOT exist.
while [ -f run-NPT+posres_${i}/confout.gro ]
do
start=run-NPT+posres_${i}/confout.gro
((i+= 1))
done

#We run the simulation with the usual combination of 
#'grompp' and 'mdrun' steps in a dedicated directory
#created here.
mkdir run-NPT+posres_${i}
cd run-NPT+posres_${i}

gmx grompp -f ../run-NPT+posres.mdp -c ../${start} -p ../complex.top -r ../solv.gro -o topol.tpr -maxwarn 1 >& grompp.out
gmx mdrun -v -nt 4 -s topol.tpr -o traj.trr -e ener.edr -g md.log -c confout.gro -cpo state.cpt >& mdrun.out

#In this post-processing step of the output trajectory, we 
#extract individual 'snapshots' of the system at regular 
#time intervals (specified with '-dt' in picoseconds). 
#These will be used later to start short simulations with 
#a high time resolution that will be used to generate the 
#actual results.
#The snapshots will be stored as separate files in the 'gro' 
#format in a designated directory created here. The files 
#will be numbered and called 
#state_0.gro, state_1.gro, state_2.gro etc.
mkdir snapshots
gmx trjconv -s topol.tpr -f traj.trr -pbc mol -sep -dt 100.0 -ndec 8 -o snapshots/state_.gro << STOP >& trjconv.out
0
STOP

#exit the newly created directory
cd ..