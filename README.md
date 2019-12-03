# Saiz *et al* (2020) 

This repository contains the data and code associated with the article by Saiz *et al.*, *Coupling between lineage size and cell fate specification underlies developmental robustness* soon to be deposited in [bioRXiv](https://www.biorxiv.org).

The repository tries to follow a standard structure, mostly based on the [Cookiecutter Data Science](https://drivendata.github.io/cookiecutter-data-science/#cookiecutter-data-science) project and other best practices recommendations. This repository was made post-hoc, in 2019, based on rawer code written between 2017 and 2019 (sorry - trying to be better).

## Repository structure:  

* README.md - you're reading it
* .Rprofile - contains global settings and options for RStudio.
* regdev_project.Rproj - RStudio project. If open with RStudio in the main directory where all files are, it should allow to reproducibly run all scripts, analyses and create graphs. If you clone the repository, it should work out of the box. 
* /packrat - contains logs and package library for this project, as created with the [packrat](http://rstudio.github.io/packrat/) R package.
* /data
   * /uncorfiles - original data files generated by MINS segmentation
   * /corfiles - manually curated files, corrected for over and undersegmentation, as described in [Morgani *et al.*, (2018) *Dev Bio*](https://doi.org/10.1016/j.ydbio.2018.06.017).
   * /raw - unprocessed files resulting from running ```*_read.R``` scripts on corresponding /corfiles  (combine corfiles and count cell numbers)
   * /processed - processed files where data has been transformed and cells classified
   * /interim - intermediate tables to be read in by other scripts
* /figures - contains plots generated during analyis
* /results - contains tables and results generated during analysis
* /notebooks - contains interactive Jupyter notebooks detailing the main the data transformations
   * Z-correction compares raw with transformed data after correcting for fluorescence decay along teh Z-axis
   * nanogata-tx describes the conversion of NANOG and GATA6 intensitiy levels obtained with different antibodies to a common scale
* /references - contains the experimental reference files (metadata), generated manually as logs over the course of experimentation
   * ```*_exp_ref.csv``` are files with experimental data for each experiment, such as experimental date, experimental group, etc
   * ```*_if.csv``` are files containing details about the antibodies used to stain the corresponding embryos in each experiment
* /src - contains source code to use in the project
   * /data - scripts that read in and transform data tables
      - ```*_read.R``` read in files in ```./data/corfiles```, combine them, clean them and count cells. Generate ```*-raw.csv``` files stored in ```./data/raw```.
      - ```*_tx.R``` apply transformations to data generated by ```*_read.R```. Generate ```*-tx.csv``` files stored in ```./data/interim```.
      - ```*_classify.R``` classify cells in ```*_tx.csv``` data tables into populations as necessary. Generate ```*-processed.csv``` files stored in ```./data/processed```.
      - ```*_counter.R``` count number of cells in each lineage in ```*-processed.csv``` data tables and writes them to ```./data/processed```. 
   * /functions - scripts defining functions used throughout analysis.
   * / models - scripts that fit models to transform some of the data.

