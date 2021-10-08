## Project Description

For this project I am analyzing the GC content and examining the locations of CpG islands on the 1st chromosome of the human genome.

To analyze the data all that is necessary is to have R installed. The Rmd file should automatically download any packages that are necessary to conduct the analysis, as well as download the necessary data from online.
In case, an error should occur, the following code can be run in order to manually install the packages using R commands. Update all packages if prompted.

```{r}
#install bioconductor
if (!requireNamespace("BiocManager")) 
     install.packages("BiocManager") 

#install bioconductor packages
BiocManager::install()
BiocManager::install(c("GenomicRanges", "GenomicFeatures", "BSgenome", "BSgenome.Hsapiens.UCSC.hg19", "Biostrings"))

#install R packages
installed_pkgs <- row.names(installed.packages())
pkgs <- c("tinytex", "kableExtra", "tidyverse","devtools")
for(p in pkgs){
	if(!(p %in% install_pkgs)){
		install.packages(p, repos = "http://cran.us.r-project.org")
	}
}
```
   
## Execute the analysis

To execute the analysis, from the project folder you can run 

``` bash
Rscript -e "rmarkdown::render('markdownHW.Rmd')"
```

This will create a file called `markdownHW.pdf` output in your directory that contains a report of the findings. Additionally, it will also create 2 folders necessary for caches and figure formatting.
Please note the process takes about 2-3 minuts as the analysis involves large files.
