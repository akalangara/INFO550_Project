# rule for making final report
FinalReport.pdf: Figs/seq_logo.png Figs/pca_plot.png Processed_Data/top_30.csv FinalReport.Rmd
	Rscript -e "rmarkdown::render('FinalReport.Rmd', quiet = TRUE)"

# rule for making the seq logo plot
Figs/seq_logo.png: Processed_Data/all_flu_reads_sorted.bam R/03_make_seqLogo.R
	Rscript R/03_make_seqLogo.R

# rule for making the pca plot
Figs/pca_plot.png: Processed_Data/all_flu_reads_sorted.bam Processed_Data/maxV.txt R/04_make_PCA.R
	Rscript R/04_make_PCA.R

# rule for making the top 30 most influential pos
Processed_Data/top_30.csv: Processed_Data/all_flu_reads_sorted.bam Processed_Data/maxV.txt R/04_make_PCA.R
	Rscript R/04_make_PCA.R

# rule for creating sorted and indexed bam file after alignment
Processed_Data/all_flu_reads_sorted.bam: Processed_Data/all_flu_reads.bam Bash/02_index_BAM.sh
	chmod +x Bash/02_index_BAM.sh && \
	Bash/02_index_BAM.sh

# rule for creating bam file of all reads
Processed_Data/all_flu_reads.bam: Processed_Data/flu_ref.fa Processed_Data/all_flu_reads.fa R/01_make_alignment.R
	Rscript R/01_make_alignment.R && \
	mv flu* Processed_Data

# rule for creating fasta file of all reads
Processed_Data/all_flu_reads.fa: Raw_Data/IND_flu_reads.fa Raw_Data/AUS_flu_reads.fa Raw_Data/CAN_flu_reads.fa Bash/00_prep_fasta.sh
	chmod +x Bash/00_prep_fasta.sh && \
	Bash/00_prep_fasta.sh

# rule for creating reference sequence
Processed_Data/flu_ref.fa: Raw_Data/IND_flu_reads.fa Raw_Data/AUS_flu_reads.fa Raw_Data/CAN_flu_reads.fa Bash/00_prep_fasta.sh
	chmod +x Bash/00_prep_fasta.sh && \
	Bash/00_prep_fasta.sh

# rule to reset to just raw data
clean:
	rm -r Processed_Data/* && \
	rm Figs/seq_logo.png Figs/pca_plot.png && \
	rm FinalReport.pdf