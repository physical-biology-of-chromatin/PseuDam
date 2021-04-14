# Contributing

When contributing to this repository, please first discuss the change you wish to make via issue,
email, or on the [ENS-Bioinfo channel](https://matrix.to/#/#ens-bioinfo:matrix.org) before making a change. 

## Project organisation

The `LBMC/nextflow` project is structured as follow:
- all the code is in the `src/` folder
- scripts downloading external tools should download them in the `bin/` folder
- all the documentation (including this file) can be found int he `doc/` folder
- the `data` and `results` folders contain the data and results of your piplines and are ignored by `git`

## Code structure

The `src/` folder is where we want to save the pipline (`.nf`) script. This folder also contains:
- the `src/install_nextflow.sh` to install the nextflow executable at the root of the project.
- some pipelines examples (like the one build during the nf_pratical)
- the `src/nextflow.config` global configuration file which contains the `docker`, `singularity`, `psmn` and `ccin2p3` profiles.
- the `src/nf_modules` folder contains per tools `main.nf` modules with predefined process that users can imports in their projects with the [DSL2](https://www.nextflow.io/docs/latest/dsl2.html)

But also some hidden folders that users don't need to see when building their pipeline:
- the `src/.docker_modules` contains the recipies for the `docker` containers used in the `src/nf_modules/<tool_names>/main.nf` files
- the `src/.singularity_in2p3` and `src/.singularity_psmn` are symbolic links to the shared folder where the singularity images are downloaded on the PSMN and CCIN2P3 

# Proposing a new tool

Each tool named `<tool_name>` must have two dedicated folders:

- `src/nf_modules/<tool_name>` where users can find `.nf` files to include
- `src/.docker_modules/<tool_name>/<version_number>` where we have the `.Dockerfile` to construct the container used in the `main.nf` file

## `src/nf_module` guide lines

We are going to take the `fastp`, `nf_module` as an example.

The `src/nf_modules/<tool_name>` should contain a `main.nf` file that describe at least one process using `<tool_name>`

### container informations

The first two lines of `main.nf` should define two variables
```
version = "0.20.1"
container_url = "lbmc/fastp:${version}"
```

we can then use the `container_url` definition in each `process` in the `container` attribute.
In addition to the `container` directive, each `process` should have one of the following `label` attributes (defined in the `src/nextflow.config` file)
- `big_mem_mono_cpus`
- `big_mem_multi_cpus`
- `small_mem_mono_cpus`
- `small_mem_multi_cpus`

```
process fastp {
  container = "${container_url}"
  label = "big_mem_multi_cpus"
  ...
}
```

### process options

Before each process, you shoud declare at least two `params.` variables:
- A `params.<process_name>` defaulting to `""` (empty string) to allow user to add more commmand line option to your process without rewritting the process definition
- A `params.<process_name>_out` defaulting to `""` (empty string) that define the `results/` subfolder where the process output should be copied if the user want to save the process output

```
params.fastp = ""
params.fastp_out = ""
process fastp {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  if (params.fastp_out != "") {
    publishDir "results/${params.fastp_out}", mode: 'copy'
  }
  ...
  script:
"""
fastp --thread ${task.cpus} \
${params.fastp} \
...
"""
}
```

The user can then change the value of these variables:
- from the command line `--fastp "--trim_head1=10"``
- with the `include` command within their pipeline: `include { fastq } from "nf_modules/fastq/main" addParams(fastq_out: "QC/fastq/")
- by defining the variable within their pipeline: `params.fastq_out = "QC/fastq/"

### `input` and `output` format

You should always use `tuple` for input and output channel format with at least:
- a `val` containing variable(s) related to the item
- a `path` for the file(s) that you want to process

for example:
```
process fastp {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.fastp_out != "") {
    publishDir "results/${params.fastp_out}", mode: 'copy'
  }

  input:
  tuple val(file_id), path(reads)

  output:
    tuple val(file_id), path("*.fastq.gz"), emit: fastq
    tuple val(file_id), path("*.html"), emit: html
    tuple val(file_id), path("*.json"), emit: report
...
```

Here `file_id` can be anything from a simple identifier to a list of several variables.
So you have to keep that in mind if you want to use it to define output file names (you can test for that with `file_id instanceof List`).

If you want to use information within the `file_id` to name outputs in your `script` section, you can use the following snipet:

```
  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
```
and use the `file_prefix` variable.

The rational behind taking a `file_id` and emitting the same `file_id` is to facilitate complex channel operations in pipelines without having to rewrite the `process` blocks.

### dealing with paired-end and single-end data

Fastq files opened with `channel.fromFilePairs( params.fastq )` create item of the following shape:

```
[file_id, [read_1_file, read_2_file]]
```

To make this call more generic, we can use the `size: -1` option, and accept arbitrary number of associated fastq file:

```
channel.fromFilePairs( params.fastq, size: -1 )
```

will thus give `[file_id, [read_1_file, read_2_file]]` for paired-end data and `[file_id, [read_1_file]]` for single-end data

```
...
  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
  if (reads.size() == 2)
  """
  fastp --thread ${task.cpus} \
    ${params.fastp} \
    --in1 ${reads[0]} \
    --in2 ${reads[1]} \
    --out1 ${file_prefix}_R1_trim.fastq.gz \
    --out2 ${file_prefix}_R2_trim.fastq.gz \
    --html ${file_prefix}.html \
    --json ${file_prefix}_fastp.json \
    --report_title ${file_prefix}
  """
  else if (reads.size() == 1)
  """
  fastp --thread ${task.cpus} \
    ${params.fastp} \
    --in1 ${reads[0]} \
    --out1 ${file_prefix}_trim.fastq.gz \
    --html ${file_prefix}.html \
    --json ${file_prefix}_fastp.json \
    --report_title ${file_prefix}
  """
...
```

## `src/.docker_modules` guide lines

