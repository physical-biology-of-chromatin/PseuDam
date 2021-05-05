version = "1.0.2"
container_url = "lbmc/flexi_splitter:${version}"

params.split = ""
params.split_out = ""


workflow split {
  take:
    reads
    config
  main:
    split_file(reads, config)
    split_file.out.fastq.view()

    

  emit:
    fastq = split_file.out.fastq
}

process split_file {
  // You can get an example of config file here:
  // src/nf_modules/flexi_splitter/marseq_flexi_splitter.yaml
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_prefix"
  if (params.split_out != "") {
    publishDir "results/${params.split_out}", mode: 'copy'
  }

  input:
  tuple val(file_id), path(reads)
  tuple val(config_id), path(config)

  output:
  tuple val(file_id), path("results/*"), emit: fastq

  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }

  if (reads.size() == 2)
  """
  flexi_splitter ${params.split} -n 2 -f ${reads[0]},${reads[1]} -o split -c ${config}
  rm -Rf split/unassigned
  mkdir -p results/
  find split -type "f" | \
    sed -E "s|(.*/split/(.*)/(.*))|\1 \2_\3|g" |
    awk '{system("mv "\$1" results/"\$2)}'
  rm -Rf split
  """
  else
  """
  flexi_splitter ${params.split} -n 1 -f ${reads[0]} -o split -c ${config}
  rm -Rf split/unassigned
  mkdir -p results/
  find split -type "f" | \
    sed -E "s|(.*/split/(.*)/(.*))|\1 \2_\3|g" |
    awk '{system("mv "\$1" results/"\$2)}'
  rm -Rf split
  """
}
