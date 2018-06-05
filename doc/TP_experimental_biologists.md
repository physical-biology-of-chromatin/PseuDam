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

# TP for experimental biologists

The Goal of this practical is to learn how to build your own pipeline with nextflow and using the tools already *wrapped* in Docker and SGE.
For this we are going to build a small RNASeq analysis pipeline that should run the following steps:

- remove Illumina adaptors
- trim reads by quality
- build the index of a reference genome
- estimate the number of RNA fragments mapping to the transcript of this genome

## Initialize your own project

You are going to build a pipeline for you or your team. So the first step is to create your own project.

Instead of reinventing the wheel, you can use the [pipelines/nextflow](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow) as a template.
To easily do so, go to the [pipelines/nextflow](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow) repository and click on the [**fork**](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow/forks/new) button.

In git, the [action of forking](https://git-scm.com/book/en/v2/GitHub-Contributing-to-a-Project) means that you are going to make your own private copy of a repository. You can then write modifications in your project, and if they are of interest for the source repository (here [pipelines/nextflow](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow)) create a merge request. Merge request are send to the source repository to ask the maintainers to integrate modifications.


