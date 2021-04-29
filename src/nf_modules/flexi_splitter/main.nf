version = "1.0.2"
container_url = "lbmc/flexi_splitter:${version}"

params.split = ""
params.split_out = ""

process split {
  // You can get an example of config file here:
  // src/nf_modules/flexi_splitter/marseq_flexi_splitter.yaml
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_prefix"
  if (params.split_out != "") {
    publishDir "results/${params.split_out}", mode: 'copy'
  }

  input:
  tuple val(file_id), path(reads)
  tuple val(config_id), path(config)

  output:
  tuple val(file_id), path("*"), emit: fastq

  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
  def whitelist_param = ""
  if (whitelist_id != "NO CONFIG"){
    whitelist_param = "-w ${white_list}"
  }

  if (reads.size() == 2)
  """
  flexi_splitter ${params.split} -f ${reads[0]} ${read[1]} -c ${config} -o split
  """
  else
  """
  flexi_splitter ${params.split} -f ${reads[0]} -c ${config} -o split
  """
}
