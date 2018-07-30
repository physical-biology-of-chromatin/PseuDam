/*
* multiqc :
* Imputs : report files
* Output :  multiqc report
*/

/*                      MultiQC                                     */

process multiqc {
  tag "$report.baseName"
  publishDir "results/fastq/multiqc/", mode: 'copy'
  cpus = 1

  input:
    file report from fastqc_report.collect()

  output:
    file "*multiqc_*" into multiqc_report

  script:
"""
multiqc -f .
"""
}

