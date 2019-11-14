# Saiz *et al* (2020) 

This repository contains the data and code associated with the article by Saiz *et al.*, *INSERT CURRENT TITLE VERSION HERE* **JOURNAL** ISSUE, PPs (2020) [LINKTO ARTICLE](WEBSITE HERE).

The repository tries to follow a standard structure, broadly based on the [Cookiecutter Data Science](https://drivendata.github.io/cookiecutter-data-science/#cookiecutter-data-science) project and other best practices recommendations. This repository was made post-hoc, in 2019, based on rawer code written between 2017 and 2019 (sorry - trying to be better).

## Repository structure:  

* README.md - you're reading it
* /packrat - contains logs and package library for this project, as created with the [packrat](http://rstudio.github.io/packrat/) R package.
* /data - contains raw data files and manually curated files (corfiles) with immunofluorescence quantification data for all embryos analyzed in the project
   * /corfiles
   * /rawfiles
* /figures - contains plots generated during analyis
* /results - contains tables and results generated during analysis
* /notebooks - contains interactive Jupyter notebooks taking you through my data analysis process
* /references - contains the experimental reference files (metadata), generated manually as logs over the course of experimentation
   * \*_exp_ref.csv are files with experimental data for each experiment, such as experimental date, experimental group, etc
   * \*_if.csv are files containing details about the antibodies used to stain the corresponding embryos in each experiment
* /src - contains source code to use in the project

