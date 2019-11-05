# nextflow pipeline

This repository is a template and a library repository to help you build nextflow pipeline.
You can fork this repository to build your own pipeline.
To get the last commits from this repository into your fork use the following commands:

```sh
git remote add upstream gitlab_lbmc:pipelines/nextflow.git
git pull upstream master
```

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

[you can follow them here.](doc/getting_started.md)

## Available tools

| tool | nf module | docker module | psmn module |
|------|:---------:|:-------------:|:----------:|
BEDtools | ok | ok | ok
BFCtools |**no**  | ok | ok
bioawk |**no**  | ok | ok
Bowtie | ok | ok | **no**
Bowtie2 | ok | ok | ok
BWA | ok | ok | ok
canu | ok | ok | ok
cutadapt | ok | ok | ok
deepTools | ok | ok | ok
fastp | ok | ok | ok
FastQC | ok | ok | ok
file_handle | **no** | ok | ok
GATK | **no** | ok | ok
HISAT2 | ok | ok | ok
HTSeq | ok | ok | ok
Kallisto | ok | ok | ok
MACS2 | ok | ok | ok
MultiQC | ok | ok | ok
MUSIC | ok | ok | ok
picard | **no** | ok | ok
pigz | **no** | ok | ok
RSEM | ok | ok | ok
Salmon | **no** | ok | ok
sambamba | ok | ok | ok
samblaster | ok | ok | ok
SAMtools | ok | ok | ok
SRAtoolkit | ok | ok | ok
subread | **no** | ok | ok
TopHat | **no** | ok | ok
Trimmomatic | **no** | ok | ok
UMItools  | **no** | ok | ok
UrQt | ok | ok | ok


## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow/tags).

## Authors

* **Laurent Modolo** - *Initial work*

See also the list of [contributors](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow/graphs/master) who participated in this project.

## License

This project is licensed under the CeCiLL License- see the [LICENSE](LICENSE) file for details
