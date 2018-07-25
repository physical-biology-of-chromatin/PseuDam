/*
* multiqc :
* Imputs : report files
* Output :  multiqc report
*/

/*                      MultiQC                                     */

process multiqc {
  tag "$repport.baseName"
  publishDir "results/fastq/multiqc/", mode: 'copy'
  cpus = 1

  input:
    file repport from fastqc_repport.collect()

  output:
    file "*multiqc_*" into multiqc_report

  script:
"""
multiqc -f .
"""
}

