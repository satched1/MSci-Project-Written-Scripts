#PBS -N annealing
#PBS -l select=1:ncpus=8:mem=16gb
#PBS -l walltime=1:00:00

module load lammps/11Aug17

cd $PBS_O_WORKDIR

# Loop over input files
for i in {0..201..5}; do
    echo "Running LAMMPS for frame_$i"
    mpiexec lammps -in simulation_inputs/input_frame_${i}.in
done
