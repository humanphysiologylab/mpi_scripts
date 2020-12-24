#!/bin/sh
#SBATCH --time 48:00:00
#SBATCH --partition RT_study
#SBATCH -n128 -N8
#SBATCH --job-name maleckar
#SBATCH --comment Atrial_genetic_algorithms

# cd /home/common/pikunov.av/pypoptim/src/model_ctypes/_koivumaki
# make clean && make

cd /home/common/pikunov.av/pypoptim/mpi_scripts/

mpiexec python mpi_script.py ../configs/maleckar/voigt/config_G3C1.json > out_G3C1.log
mpiexec python mpi_script.py ../configs/maleckar/voigt/config_G2C2.json > out_G2C2.log
mpiexec python mpi_script.py ../configs/maleckar/voigt/config_G4C8.json > out_G4C8.log
mpiexec python mpi_script.py ../configs/maleckar/voigt/config_G1C1.json > out_G1C1.log
