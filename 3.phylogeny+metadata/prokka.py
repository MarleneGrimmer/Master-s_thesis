#!/usr/bin/python

import sys, os, glob
import multiprocessing
from multiprocessing.pool import ThreadPool

def run_prokka(f):
        fname = f.split("/")[-1].strip()
        fname = fname.split(".fna")[0].strip()
        fsl = fname.split("_")[1].strip() #split the filename at _ and take second part (number)
        os.system("apptainer run /proj/nobackup/carroll_hpc2n/containers/prokka_1.14.6--pl5321hdfd78af_5.sif prokka --outdir prokka_results/" + fname + "_prokka_results --kingdom Bacteria --cpus 4 --addgenes --locustag " + fsl + " --centre --compliant --prefix " + fname + " " + f)  #running prokka on files
if __name__ == "__main__": #parallelize execution of prokka on input files
        pool = multiprocessing.Pool(int(sys.argv[1])) #creates multiprocessing pool
        tasks = glob.glob("/proj/nobackup/carroll_hpc2n/marlene/high_quality_fna/*.fna") #tasks: all files .fna
        pool.map(run_prokka, tasks) #distributes tasks to modules
