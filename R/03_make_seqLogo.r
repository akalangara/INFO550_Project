# load needed packages
library(BiocManager)
#BiocManager::install("Rsamtools")
library(Rsamtools)
#install.packages('devtools')
#devtools::install_github("omarwagih/ggseqlogo")
library(ggseqlogo)
library(ggplot2)

here::i_am('R/03_make_seqLogo.R')

# load aligned sequences using bioconductor for efficiency
bf <- BamFile(here::here("Processed_Data","all_flu_reads_sorted.bam"))
seqinfo(bf)
sl <- seqlengths(bf)
all_gr <- GRanges("IND_cds:AVZ00479",IRanges(1, sl["IND_cds:AVZ00479"]))
all_reads <- scanBam(bf, param=ScanBamParam(what=c("qwidth","qname","seq"), which=all_gr))

# create position frequency matrix
pfm_all <-consensusMatrix(all_reads[[1]]$seq)

# use position freq matrix to create position probability matrix
# a ppm is the input for seqLogo
num_seq = length(all_reads[[1]]$seq)
convert2PPM <- function(pfm) (pfm / num_seq)
ppm_all <- as.data.frame(convert2PPM(pfm_all[1:4,])) #1:4 so just looking at A,G,T,C

# subset to positions of high variation (max prob for a nt is <= 90%)
maxV <- sapply(ppm_all, max)
i_moreThan90 <- which(maxV <= 0.90)
highV <- pfm_all[1:4,i_moreThan90]

# max probabilities will be used in PCA analysis later
write.table(maxV, here::here("Processed_Data","maxV.txt"), sep="\t")

# generate and export seqLogo plot
colnames(highV) <- i_moreThan90
seq_fig<-ggplot() + geom_logo(highV, method='prob')
seq_fig$scales$scales[[1]] <- scale_x_continuous(name ="Position", 
                                              breaks= seq(1,26,by=1), 
                                              labels=i_moreThan90)
png(here::here("Figs", "seq_logo.png"),width=1000)
seq_fig
dev.off()