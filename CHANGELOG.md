# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.6] - 2018-08-23
### Added
- Added `src/training_dataset.nf` to build a small training dataset from NGS data

### Changed
- the structure of `src/nf_modules`: the `tests` folder was removed

## [0.2.5] - 2018-08-22
### Added
- This fine changelog

### Changed
- the structure of `src/nf_modules`: the `tests` folder was removed


## [0.2.4] - 2018-08-02
### Changed
- add `paired_id` variable in the output of every single-end data processes to match the paired output


## [0.2.3] - 2018-07-25
### Added
- List of tools available as nextflow, docker or sge module to the `README.md`


## [0.2.2] - 2018-07-23
### Added
- SRA module from cigogne/nextflow-master 52b510e48daa1fb7


## [0.2.1] - 2018-07-23
### Added
- List of tools available as nextflow, docker or sge module


## [0.2.0] - 2018-06-18
### Added
- `doc/TP_computational_biologists.md`
- Kallisto/0.44.0

### Changed
- add `paired_id` variable in the output of every paired data processes
- BEDtools: fixes for fasta handling
- UrQt: fix git version in Docker


## [0.1.2] - 2018-06-18
### Added
- `doc/tp_experimental_biologist.md` and Makefile to build the pdf
- tests files for BEDtools

### Changed
- Kallisto: various fixes
- UrQt: improve output and various fixes

### Removed
- `src/nf_test.config` modules have their own `.config`


## [0.1.2] - 2018-06-18
### Added
- `doc/tp_experimental_biologist.md` and Makefile to build the pdf
- tests files for BEDtools

### Changed
- Kallisto: various fixes
- UrQt: improve output and various fixes

### Removed
- `src/nf_test.config` modules have their own `.config`


## [0.1.0] - 2018-05-06
This is the first working version of the repository as a nextflow module repository

