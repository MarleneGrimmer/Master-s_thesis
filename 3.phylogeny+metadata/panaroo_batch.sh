#!/bin/bash
#SBATCH -A hpc2n2023-117
#SBATCH -n 16
#SBATCH --output=panaroo.out
#SBATCH --error=panaroo.err
#SBATCH --time=1-00:00:00

module purge > /dev/null 2>&1

apptainer run /proj/nobackup/carroll_hpc2n/containers/panaroo_1.3.4--pyhdfd78af_0.sif panaroo -i /proj/nobackup/carroll_hpc2n/marlene/prokka/prokka_results/*/*.gff -o panaroo_results/ --clean-mode strict -a core --aligner mafft --core_threshold 0.95 -t 4 -f 0.5
