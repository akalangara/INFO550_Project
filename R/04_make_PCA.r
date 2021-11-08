#load needed packages
library(BiocManager)
library(Rsamtools)
library(stringr)
library(ggplot2)
library(caret)

here::i_am('R/04_make_PCA.R')

# load aligned sequences
bf <- BamFile(here::here("Processed_Data","all_flu_reads_sorted.bam"))
seqinfo(bf)
sl <- seqlengths(bf)
all_gr <- GRanges("IND_cds:AVZ00479",IRanges(1, sl["IND_cds:AVZ00479"]))
all_reads <- scanBam(bf, param=ScanBamParam(what=c("qwidth","qname","seq"), which=all_gr))

# load max probabilities at each position
maxV <- read.table(here::here("Processed_Data","maxV.txt"), header=TRUE, sep="\t")

# determine number of seqs for each country 
# (some were not aligned so not the same as query)
seq_names = all_reads[["IND_cds:AVZ00479:1-1701"]][["qname"]]
num_can <- table(str_count(seq_names, "CAN_"))[2]
num_aus <- table(str_count(seq_names, "AUS_"))[2]
num_ind <- table(str_count(seq_names, "IND_"))[2]

# create crude data matrix 
# row and columns were filled this way for efficiency
# row: 1 for each position in the HA sequence
# col: 1 for every sequence, identified by country and number
# Ex: row 8 contains the 8th basepair in all sequences
# EX: column labeled CAN 2 contains the sequence for the 2nd Canadian sequence
nts <- matrix(all_reads[["IND_cds:AVZ00479:1-1701"]][["seq"]])
seq_list <- nts[,1]
seq_list <-paste(seq_list,collapse ="")
seq_list <- unlist(strsplit(seq_list, split=""))
data.matrix <- matrix(seq_list, nrow=1209, ncol=644,byrow=FALSE)
colnames(data.matrix) <- c(
     paste("CAN", 1:num_can, sep=""),
     paste("AUS", 1:num_aus, sep=""),
     paste("IND", 1:num_ind, sep=""))
rownames(data.matrix) <- paste(1:1209, sep="")

# transpose the crude matrix (rows and columns have switched)
# proper orientation for analysis
t.data.matrix <- t(data.matrix)

# use maxV to drop positions with no variability
# PCA doesn't run if there are positions (variables) that are constant
i_max <- which(maxV == 1)
dm.df <- as.data.frame(t.data.matrix)
dropped_df <- dm.df[-i_max]

# create and fill dummy variables
# PCA is for continuous data but we are trying to run it on A, G, T, and C
# For each non varying position, there is a binary variable for A, G, T, and C (# of columns has 3x)
dmy <- dummyVars(" ~ .", data = dropped_df)
dmy_df <- data.frame(predict(dmy, newdata = dropped_df))
dmy_matrix <-as.matrix(dmy_df)

# conduct PCA and determine PC variation
pca <- prcomp(dmy_matrix, scale=TRUE)
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)


# create dataframe for mapping (sample=country, X=PC1, Y=PC2)
pca.data <- data.frame(Sample=substring(rownames(pca$x), 1, 3),
                       X=pca$x[,1],
                       Y=pca$x[,2])

# generate and export PCA plot
pca_plot <- ggplot(data=pca.data, aes(x=X, y=Y, color=Sample)) +
     geom_point() +
     xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
     ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
     theme_bw()

png(here::here("Figs", "pca_plot.png"))
pca_plot
dev.off()

# get top 30 most influential positions
loading_scores <- pca$rotation[,1]
pos_scores <- abs(loading_scores) ## get the magnitudes (ignoring +/-)
pos_score_ranked <- sort(pos_scores, decreasing=TRUE)
top_30_pos <- names(pos_score_ranked[1:30])

top_30_score <- pca$rotation[top_30_pos,1]
top_30_pos <- substr(top_30_pos, 3, nchar(top_30_pos)) #gets rid of automatically added X. in front of name

top_30 <- data.frame(top_30_pos, top_30_score)

# export top 30 and their scores
write.csv(top_30,here::here("Processed_Data","top_30.csv"), row.names = FALSE)