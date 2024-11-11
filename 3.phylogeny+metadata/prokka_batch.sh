#!/bin/bash
#SBATCH -A hpc2n2023-117
#SBATCH -n 28
#SBATCH --output=prokka.out
#SBATCH --error=prokka.err
#SBATCH --time=1-00:00:00

module purge > /dev/null 2>&1
module load GCC/12.3.0
module load OpenMPI/4.1.5
module load Biopython/1.83

python prokka.py 7
