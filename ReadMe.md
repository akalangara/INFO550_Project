## Project Description

Influenza shows different seasonal patterns depending on the location. Tropical regions may experience a year-long flu season, while temperate regions experience a limited flu season that occurs in cold and dry months. For this project, I aim to explore the differences between tropical and temperate influenza strains by analyzing the nucleotide positional variations between the two types of strains.

## Necessary Programs

For this analysis it is necessary to have R and samtools installed.

You can use the following to install samtools
```bash
sudo apt install samtools
```

## Restore Package Environment
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

## Execute Analysis
To execute the analysis, navigate to the `INFO550_Project` directory. You can then run 

``` bash
make
```

This will create in the `INFO550_Project` directory, a file called `FinalReport.pdf` which contains a report of the findings.

Please note the process takes about 2-3 minutes as the analysis involves many packages and large files.


## Reset to Raw Data

To remove all created files and return to just having the raw data, you can run

``` bash
make clean
```