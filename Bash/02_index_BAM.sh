#create a sorted and indexed bam file necessary for use with bioconductor
samtools sort Processed_Data/all_flu_reads.bam -o Processed_Data/all_flu_reads_sorted.bam && \
samtools index Processed_Data/all_flu_reads_sorted.bam