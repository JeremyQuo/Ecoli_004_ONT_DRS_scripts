# Ecoli_004_ONT_DRS_scripts
This GitHub repository contains various scripts,
including plotting scripts for all main figures, some custom scripts for the QC part, and the commands used in our paper. Additionally, it serves as a comprehensive resource for replicating our analyses and further exploring our methodologies.

## Project Overview

This project focuses on the analysis of ONT DRS data for E. coli. The provided scripts help in the visualization, quality control, and overall data processing required for our study.

In the QC section, we provide three custom scripts to calculate mapped information from BAM files and FASTQ files, namely 
[count_base_num_each_gene.py](1_QC%2Fpython_script%2Fcount_base_num_each_gene.py), 
[count_read_num_each_gene.py](1_QC%2Fpython_script%2Fcount_read_num_each_gene.py) and 
[count_map_unmap.py](1_QC%2Fpython_script%2Fcount_map_unmap.py)

## Installation and Dependencies

To run the scripts in this repository, you will need the following software and libraries:

- Python 3.8+
  - plotnine
  - pysam
  - tqdm
  - pandas
  - numpy
  - sklearn
  - matplotlib
- R
  - ggplot2
  - waffle
 
## nanoSundial

For instructions on using nanoSundial, please refer to this [documentation](https://github.com/lrslab/nanoSundial).
