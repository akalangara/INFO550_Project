# load needed packages
library(BiocManager)
#BiocManager::install("Rsubread")
library(Rsubread)


here::i_am('R/01_make_alignment.R')

## build reference needed for alignment
buildindex(basename="flu",
           reference= here::here("Processed_Data","flu_ref.fa"))
## align sequences so that positions match up
align.stat <- align(index="flu",
                    readfile1= here::here("Processed_Data","all_flu_reads.fa"),
                    output_file= here::here("Processed_Data","all_flu_reads.bam"),
                    type='dna')
