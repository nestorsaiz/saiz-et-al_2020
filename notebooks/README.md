This folder contains computational notebooks detaling some of the data transformation steps performed in the study. These steps are described in detail in the Methods section of the manuscript, but the notebooks show the code used for the correction steps and their effect on the data:

* ```_data-processing-wf.pdf``` is a flowchart describing the processing steps for each subset of the data in the study.
* ```GFP-classification.ipynb``` details the automatic classificaton of GFP+ and GFP- cells in embryo-embryo aggregation chimeras.
* ```Z-correction.ipynb``` details the correction for signal decay along the Z axis. This procedure builds on the steps described in [Saiz *et al* (2016)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4828170), [Saiz *et al* (2016)](https://www.nature.com/articles/ncomms13463) and [Morgani *et al* (2018)](https://doi.org/10.1016/j.ydbio.2018.06.017).
* ```nanogata-tx.ipynb``` describes the transformation of fluorescence values obtained with anti-NANOG(rat) and anti-GATA6(rb) to values equivalent to those obtained with anti-NANOG(rb) and anti-GATA6(gt) antibodies, as described in the manuscript.
