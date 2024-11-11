nextflow.enable.dsl=2

workflow {
    // define the parameters
    input = Channel
        .fromPath("/proj/nobackup/carroll_hpc2n/marlene/fna/*.fna")
        .map{file -> tuple(file.simpleName, file)}
    // run the process
    antismash(input)
}

process antismash {
  container "/proj/nobackup/carroll_hpc2n/containers/antismash_6.1.1--pyhdfd78af_0.sif"
  publishDir "antiSMASH_results/", mode: "copy"
  input:
    tuple val(accession), path(contigs)
  output:
    tuple val(accession), path("${accession}_antiSMASH")
  script:
  """
  antismash --output-dir "${accession}_antiSMASH" -c ${task.cpus} --genefinding-tool prodigal-m ${contigs}
  """
}
