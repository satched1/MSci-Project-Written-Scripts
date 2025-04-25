#PBS -N annealing
#PBS -l select=1:ncpus=8:mem=16gb
#PBS -l walltime=72:00:00

module load lammps/11Aug17

cd $PBS_O_WORKDIR
mpiexec lammps -in input.in
