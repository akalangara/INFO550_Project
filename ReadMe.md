# Project Description

Influenza shows different seasonal patterns depending on the location. Tropical regions may experience a year-long flu season, while temperate regions experience a limited flu season that occurs in cold and dry months. For this project, I aim to explore the differences between tropical and temperate influenza strains by analyzing the nucleotide positional variations between the two types of strains.

## Analysis Execution Method 1: Using Docker
<details>
  <summary>Instructions</summary>
  
  ### Pulling Needed Image
  In order to be able to run the necessary container for this analysis, we need to first pull the appropriate image. You can use the following code to do so.
  ```bash
  docker pull akalangara/info550_project
  ```
  This is a chunky image so it takes around 10 minutes.
  
  ### Creating Output Folder
  When the container is run, the analysis will automatically execute to produce a report of the findings in the container. In order to be able to retrieve the report, an `Output` directory needs to be mounted. The `Output` directory does not need to be in a specific location, so long as the location is noted. 
  
  An example of how to create an `Output` directory in the Downloads directory is given below. You can create your folder anywhere and call it whatever you want so long as you remember what it is called, and its filepath (location). This can be done by navigating inside the directory and using the `pwd` command.
  ```bash
  cd ~/Downloads
  mkdir Output
  ```
  The filepath for this folder is `~/Downloads/Output`.

  ### Running Container to Execute Analysis
  To run the container and be able to access the report you can use the following pseudocode. You will need to fill in the filepath and output folder name.
  ```bash
  docker run -v ~/path/to/your/output_folder_name:/project/Output akalangara/info550_project
  ```
  For the example given previously, this is what the command would look like.
  ```bash
  docker run -v ~/Downloads/Output:/project/Output akalangara/info550_project
  ```
  Please note that generating the report takes around 10 minutes as the analysis involves many packages and large files. Once this code is run, the container will close and the report will be in your local `Output` directory. As the code is running, a number of warnings regarding packages will appear. These can be ignored.

  ### Running Container Interactively
  If you do not want to automatically run the analysis, you can use the following psuedocode to overide the entry point.
  ```bash
  docker run -v ~/path/to/your/output_folder_name:/project/Output -it akalangara/info550_project /bin/bash
  ```
  If you want to produce the full report, you can run `make` from the `/project` directory. You can also reset to raw data as noted below.
</details>

## Analysis Execution Method 2: Using Only Renv
<details>
  <summary>Instructions</summary>
  
  ### Necessary Programs
  For this analysis it is necessary to have R and samtools installed.
  You can use the following to install samtools
  ```bash
  sudo apt install samtools
  ```
  ### Restore Package Environment
  Prior to executing the analysis, navigate to the `INFO550_Project` directory and start an R session
  ```bash
  R
  ```
  In the R session, run the following to restore the package environment
  ```R
  renv::restore()
  ```
  There will be a lot of output. Please note any errors or warnings of uninstalled packages that come up.
  Quit the R session once this step is completed.
  ```R
  q()
  ```
  ### Execute Analysis
  To execute the analysis, navigate to the `INFO550_Project` directory. You can then run 
  ``` bash
  make
  ```
  This will create in the `INFO550_Project/Output` directory, a file called `FinalReport.pdf` which contains a report of the findings.
  Please note the process takes about 2-3 minutes as the analysis involves many packages and large files.
</details>

## Reset to Raw Data
  To remove all created files and return to just having the raw data, you can run
  ``` bash
  make clean
  ```
  If the entire analysis did not execute prior to cleaning, this may produce some errors. Ultimately, this should not impact the organization of files.

