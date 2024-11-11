nextflow.enable.dsl=2

workflow {
    // define the parameters
    input = Channel
        .fromPath("/proj/nobackup/carroll_hpc2n/marlene/fna/*.fna")
        .map{file -> tuple(file.simpleName, file)}
    // run the process
    gecco(input)
}

process gecco {
  container "/proj/nobackup/carroll_hpc2n/containers/gecco_0.9.8--pyhdfd78af_0.sif"
  publishDir "gecco_results/", mode: "copy"
  input:
    tuple val(accession), path(contigs)
  output:
    tuple val(accession), path("${accession}_gecco")
  script:
  """
  gecco run -g ${contigs} -o "${accession}_gecco" -j ${task.cpus} --mask --threshold 0.3
  """
}
