

[![DOI](https://zenodo.org/badge/217027606.svg)](https://zenodo.org/badge/latestdoi/217027606)


# Mathematical reconstruction of land carbon models from their numerical output: computing soil radiocarbon from <sup>12</sup>C dynamics

This repository contains all files related to the manuscript

```{yaml}
Title: Mathematical reconstruction of land carbon models from their numerical output: computing soil radiocarbon from 12C dynamics
Authors: Holger Metzler, Qing Zhu, William Riley, Alison Hoyt, Markus Müller, Carlos A. Sierra
Journal: Journal of Advances in Modeling the Earth System
```

It contains all manuscript files as well as the code to reproduce all results presented in the text. A description of the main folders is given below.

## Manuscript
This folder contains all original LaTeX files used to prepare the publication, including a BibTeX file with references and a folder with all figures.

## Data
We provide two datasets to reproduce all results from the manuscript. The file `C14Atm_NH.csv` contains atmospheric radiocarbon data for the northern hemisphere. These data are used to inform the model about the inputs of radiocarbon to the ecosystem. These atmospheric radiocarbon data are consistent with the data used to run the original ELM model.

The second dataset is the output from the ELM model used for the reconstruction. Given its size, the data was compressed and split in four different files: `JAMES.nc.7z.001` to `JAMES.nc.7z.004`. The compression was done with [7zip](https://www.7-zip.org/). On the command line, the compressed data can be extracted and merged with the command

```7z e JAMES.nc.7z.001```

## Code
There are two main files that run the code that reproduce all graphs and tables in the manuscript: `simple_models.py` and `ELM_discrete.py`. There are a number of preparation steps that are required to run the code:
* Decompress the data as explained in the previous section.
* Create a python virtual environment where all package requirements will be stored.
* Install the full set of python packages `testinfrastructure`, `LAPM`, `CompartmentalSystem`, `bgc-md`. The description of how to install these packages can be found in <https://github.com/MPIBGC-TEE/bgc-md>.
* Using pip3, install in the virtual environment the additional packages: `pandas`, `xarray`, and `netcdf4`.
* Now you can run the examples as `python simple_models.py` and as `python ELM_discrete.py`. The ouput figures will be stored in the folder `Code/Output/`.

### Updates
In case of updates, they will be available in the main repository at <https://github.com/MPIBGC-TEE/modelReconstruction>.
