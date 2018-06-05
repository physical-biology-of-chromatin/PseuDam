# nextflow pipeline

This repository is a template and a library repository to help you build nextflow pipeline.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

To run nextflow on you computer you need to have java (>= 1.8) installed.

```sh
java --version
```

To be able to easily test tools already implemented for nextflow on your computer (`src/nf_modules/` to see their list). You need to have docker installed.

```sh
docker run hello-world
```

### Installing

To install nextflow on you computer simply run the following command:

```sh
src/install_nextflow.sh
```

Then to initialise a given tools run the following command:

```sh
src/nf_modules/<tool_name>/<tool_version>/docker_init.sh
```

for example to initialise `file_handle` version `0.1.1`, run:

```sh
src/nf_modules/file_handle/0.1.1/docker_init.sh
```

## Running the tests

```sh
cd data
git clone -c http.sslVerify=false https://gitlab.biologie.ens-lyon.fr/LBMC/tiny_dataset.git
cp data/tiny_dataset/fastq/tiny_R1.fastq data/tiny_dataset/fastq/tiny2_R1.fastq
cp data/tiny_dataset/fastq/tiny_R2.fastq data/tiny_dataset/fastq/tiny2_R2.fastq
cp data/tiny_dataset/fastq/tiny_S.fastq data/tiny_dataset/fastq/tiny2_S.fastq
nextflow src/nf_test.nf -c src/nf_test.config --fastq "data/tiny_dataset/fastq/tiny_R{1,2}.fastq"
```

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow/tags). 

## Authors

* **Laurent Modolo** - *Initial work*

See also the list of [contributors](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow/graphs/master) who participated in this project.

## License

This project is licensed under the CeCiLL License- see the [LICENSE](LICENSE) file for details

