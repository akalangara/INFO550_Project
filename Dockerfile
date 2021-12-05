FROM rocker/verse:3.6.3-ubuntu18.04

# install dependencies and updates for samtools
RUN sudo apt-get update -y
RUN sudo apt-get install -y gcc
RUN sudo apt-get install -y make
RUN sudo apt-get install -y libbz2-dev
RUN sudo apt-get install -y zlib1g-dev
RUN sudo apt-get install -y libncurses5-dev 
RUN sudo apt-get install -y libncursesw5-dev
RUN sudo apt-get install -y liblzma-dev

# install HTSlib
RUN cd /usr/bin && \
wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
tar -vxjf htslib-1.9.tar.bz2 && \
cd htslib-1.9 && \
make && \
make install

# install samtools
RUN cd .. && \
wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
tar -vxjf samtools-1.9.tar.bz2 && \
cd samtools-1.9 && \
make && \
make install

# install bcftools
RUN cd .. && \
wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && \
tar -vxjf bcftools-1.9.tar.bz2 && \
cd bcftools-1.9 && \
make && \
make install

# add to path
ENV PATH=/usr/bin/bcftools-1.9:$PATH
ENV PATH=/usr/bin/samtools-1.9:$PATH
ENV PATH=/usr/bin/htslib-1.9:$PATH

# install tiny tex for publishing
RUN wget -qO- "https://yihui.org/tinytex/install-bin-unix.sh" | sh

# make a project directory
RUN mkdir /project

# set workdir to project folder
WORKDIR /project

# copy contents of local folder to project folder in container
COPY ./ /project/

# install renv, and restore package environment
ENV RENV_VERSION 0.14.0
RUN Rscript -e "Sys.setenv(R_INSTALL_STAGED = FALSE);\
install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'));\
remotes::install_github('rstudio/renv@${RENV_VERSION}');\
renv::restore()"

# make R files executable
RUN chmod +x /project/R/*.R

# make container entry point creating report
CMD make