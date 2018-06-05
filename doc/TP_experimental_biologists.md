---
title: "TP for experimental biologists"
author: Laurent Modolo [laurent.modolo@ens-lyon.fr](mailto:laurent.modolo@ens-lyon.fr)
date: 6 Jun 2018
output:
  pdf_document:
    toc: true
    toc_depth: 3
    number_sections: true
    highlight: tango
    latex_engine: xelatex
---

The Goal of this practical is to learn how to build your own pipeline with nextflow and using the tools already *wrapped*.
For this we are going to build a small RNASeq analysis pipeline that should run the following steps:

- remove Illumina adaptors
- trim reads by quality
- build the index of a reference genome
- estimate the number of RNA fragments mapping to the transcript of this genome

# Initialize your own project

You are going to build a pipeline for you or your team. So the first step is to create your own project.

## Forking

Instead of reinventing the wheel, you can use the [pipelines/nextflow](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow) as a template.
To easily do so, go to the [pipelines/nextflow](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow) repository and click on the [**fork**](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow/forks/new) button.

![fork button](img/fork.png)

In git, the [action of forking](https://git-scm.com/book/en/v2/GitHub-Contributing-to-a-Project) means that you are going to make your own private copy of a repository. You can then write modifications in your project, and if they are of interest for the source repository (here [pipelines/nextflow](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow)) create a merge request. Merge request are send to the source repository to ask the maintainers to integrate modifications.

![merge request button](img/merge_request.png)

## Project organisation

This project (and yours) follow the [guide of good practices for the LBMC](http://www.ens-lyon.fr/LBMC/intranet/services-communs/pole-bioinformatique/ressources/good_practice_LBMC)

You are now on the main page of your fork of the [pipelines/nextflow](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow). You can explore this project, all the code in it is under the CeCILL lience (in the [LICENCE](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow/blob/master/LICENSE) file).

The [README.md](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow/blob/master/README.md) file contains instructions to run your pipeline and test it's installation.

The [CONTRIBUTING.md](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow/blob/master/CONTRIBUTING.md) file contains guidelines to follow if you want to contribute to the [pipelines/nextflow](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow) (making a merge request for example).

The [data](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow/tree/master/data) folder will be the place were you store the raw data for your analysis.
The [results](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow/tree/master/results) folder will be the place were you store the results of your analysis.
Note that the content of these two folders should never be saved on git.

The [doc](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow/tree/master/doc) folder contains the documentation of this practical course.

And most interestingly for you, the [src](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow/tree/master/src) contains code to wrapp tools. This folder contains two subdirectory. A `nf_modules` and a `sge_modules` folder. 
 






