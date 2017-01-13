
BAM_FILE=/home/vincent/Desktop/Sequencing/Methods_Compare/Total/Rsubread_test/CoCo/MED24.sam
#GTF=/home/vincent/Desktop/Results/Seq/hg38/SNORD124.gtf
#GTF=/home/vincent/Desktop/Results/Seq/hg38/FabDupSan_human_ensembl_83.gtf
GTF=/home/vincent/Desktop/Sequencing/Methods_Compare/Total/Rsubread_test/CoCo/SNORD124.simplified_gap.gtf
OUTPUT_PREFIX=/home/vincent/Desktop/Sequencing/Methods_Compare/Total/Rsubread_test/test.exon.out.$1

#samtools sort -n $BAM_FILE $BAM_FILE.sorted
#BAM_FILE=/home/vincent/Desktop/Sequencing/Methods_Compare/Total/Rsubread_test/test.1000.bam.sorted.bam

Rscript /home/vincent/Desktop/Sequencing/Methods_Compare/Total/Rsubread_test/Rsubread.featureCounts_exon.R \
$BAM_FILE $GTF $OUTPUT_PREFIX

