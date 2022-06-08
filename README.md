# DamID pipeline

* This DamID pipeline is designed to treat single-end and paired-end fastq files in a novel way by using the DamID experiments property of having all reads in windows between 2 GATC sites. This experimental particularity allows us to treat the same way as RNAseq datas. We use Kallisto to pseudoalign the reads to a fasta file containing all the GATC fragments in a given reference genome.

* The pipeline executes the following steps :

* - Quality control and adapter trimming using fastp (ref)
* - Creation of a bed file of the GATC fragments with gatc_finder
* - Creation of a fasta files from the bed file with fasta_from_bed from bedtools
* - Indexing of the fragments contained in the fasta files by Kallisto
* - Pseudomapping of the reads to the fragments by Kallisto
* - Generation of a global report with multiqc ()


* requirements to launch the pipeline :
* - Docker installed

* The pipeline is launched using the following command :

* ./nextflow src/Dam_ID_analysis.nf -profile docker \
* -fastq "path/to/reads.fastq" \
* -fasta "path/to/genome.fasta" \
* -length int : mean fragments' length* (If you do not have this information use the average length of the reads)
* -std int : standard deviation of the fragments' length* (If you do not have this information use les than 10% of the mean length)

* * : nescessary only if you are processing single end datas. Kallisto computes it by itself if with paired end datas






## Authors

* **Laurent Modolo** - *Initial work*

* **Nathan Lecouvreur** - *DamID pipeline*


## License

This project is licensed under the CeCiLL License- see the [LICENSE](LICENSE) file for details
