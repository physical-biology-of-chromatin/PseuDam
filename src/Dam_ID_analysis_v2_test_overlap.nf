nextflow.enable.dsl=2

include { fastp                             } from "./nf_modules/fastp/main.nf"

include { gatc_finder as gatc_no_overlap    } from "./nf_modules/gatc_finder/main.nf" addParams(gatc_finder_out: "sites/test/no/", gatc_finder: "--overlap_size 0")
include { gatc_finder as gact_overlap       } from "./nf_modules/gatc_finder/main.nf" addParams(gatc_finder_out: "sites/test/yes/", gatc_finder: "--overlap")
include { gatc_finder as gatc_overlap_10    } from "./nf_modules/gatc_finder/main.nf" addParams(gatc_finder_out: "sites/test/10/", gatc_finder: "--overlap_size 14")
include { gatc_finder as gatc_overlap_50    } from "./nf_modules/gatc_finder/main.nf" addParams(gatc_finder_out: "sites/test/50/", gatc_finder: "--overlap_size 54")

include { index_fasta as index_test_no ; mapping_fastq as mapping_test_no       } from "./nf_modules/kallisto/main.nf"    addParams(mapping_fastq_out: "pseudo/test/no/")
include { index_fasta as index_test_yes ; mapping_fastq as mapping_test_yes     } from "./nf_modules/kallisto/main.nf"    addParams(mapping_fastq_out: "pseudo/test/yes/")
include { index_fasta as index_test_10 ; mapping_fastq  as mapping_test_10      } from "./nf_modules/kallisto/main.nf"    addParams(mapping_fastq_out: "pseudo/test/10/")
include { index_fasta as index_test_50 ; mapping_fastq  as mapping_test_50      } from "./nf_modules/kallisto/main.nf"    addParams(mapping_fastq_out: "pseudo/test/50/")

include { multiqc     as multiqc_test_no    } from "./nf_modules/multiqc/main.nf"     addParams(multiqc_out: "mapping/test/no/")
include { multiqc     as multiqc_test_yes   } from "./nf_modules/multiqc/main.nf"     addParams(multiqc_out: "mapping/test/yes/")
include { multiqc     as multiqc_test_10    } from "./nf_modules/multiqc/main.nf"     addParams(multiqc_out: "mapping/test/10/")
include { multiqc     as multiqc_test_50    } from "./nf_modules/multiqc/main.nf"     addParams(multiqc_out: "mapping/test/50/")

include { fasta_from_bed as fasta_test_no   } from "./nf_modules/bedtools/main.nf"
include { fasta_from_bed as fasta_test_yes  } from "./nf_modules/bedtools/main.nf"
include { fasta_from_bed as fasta_test_10   } from "./nf_modules/bedtools/main.nf"
include { fasta_from_bed as fasta_test_50   } from "./nf_modules/bedtools/main.nf"


params.fasta = "data/genome/S288C_reference_sequence_R64-3-1_20210421.fsa"
params.fastq = "data/reads/Dam_ID/*_{1,2}.fq"


channel
    .fromPath(params.fasta)
    .ifEmpty { error "Cannot find any fasta files matching: ${params.fasta}" }
    .map { it -> [it.simpleName, it]}
    .set {fasta_files}

channel
    .fromFilePairs(params.fastq, size: -1)
    .set {fastq_files}


workflow {
    
    gatc_no_overlap(fasta_files)
    gact_overlap(fasta_files)
    gatc_overlap_10(fasta_files)
    gatc_overlap_50(fasta_files)


    fastp(fastq_files)


    fasta_test_no(fasta_files,
                  gatc_no_overlap.out.bed)
    fasta_test_yes(fasta_files,
                   gact_overlap.out.bed)
    fasta_test_10(fasta_files,
                  gatc_overlap_10.out.bed)
    fasta_test_50(fasta_files,
                  gatc_overlap_50.out.bed)


    index_test_no(fasta_test_no.out.fasta)
    index_test_yes(fasta_test_yes.out.fasta)
    index_test_10(fasta_test_10.out.fasta)
    index_test_50(fasta_test_50.out.fasta)


    mapping_test_no(index_test_no.out.index.collect(),
                    fastp.out.fastq)
    mapping_test_yes(index_test_yes.out.index.collect(),
                     fastp.out.fastq)
    mapping_test_10(index_test_10.out.index.collect(),
                    fastp.out.fastq)
    mapping_test_50(index_test_50.out.index.collect(),
                    fastp.out.fastq)


    multiqc_test_no(mapping_test_no.out.report)
    multiqc_test_yes(mapping_test_yes.out.report)
    multiqc_test_10(mapping_test_10.out.report)
    multiqc_test_50(mapping_test_50.out.report)    
}