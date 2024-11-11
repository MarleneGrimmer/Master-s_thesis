#SBATCH -A hpc2n2023-117
#SBATCH -n 1
#SBATCH --output=htgcf.out
#SBATCH --error=htgcf.err
#SBATCH --time=05:50:00


module purge > /dev/null 2>&1
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load Biopython/1.79

# install
pip install --user /proj/nobackup/carroll_hpc2n/marlene/htgcf//htgcf_software/htgcf_beta/htgcf/

PATH=/home/m/marlene/.local/bin:$PATH

PATH=/proj/nobackup/carroll_hpc2n/marlene/htgcf/htgcf_software/htgcf_beta/mmseqs/bin:$PATH


htgcf -i /proj/nobackup/carroll_hpc2n/bgc/htgcf/htgcf_mibig.gbk -o results/proG3_rep1.tsv -j $SLURM_CPUS_ON_NODE --workdir /proj/nobackup/carroll_hpc2n/marlene/htgcf/temp --clustering-distance 0.05 --clustering-method "weighted"
