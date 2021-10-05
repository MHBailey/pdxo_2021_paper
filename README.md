# pdxo_2021_paper
Code availability and reproducibility for Guillen et al. 2021

Over all this project highlights two snakemake workflows to recreate the genomic and drug screening dataset for the following manuscript, "A breast cancer patient-derived xenograft and organoid platform for drug discovery and precision oncology" by Guillen et al. 2021. I'm one of the co-authors of this project and we made it to publication so I'm providing the code that I used to generate the foundation figures for these two sets of analysis. Many of these figures were "touched up" in Adobe illustrator but the majority of the condent is there. 

In each of the directories I provide a snakemake workflow that executes the commands. Data for this project will not be stored in this GitHub Repo. Your job is to find those files and pull them in here. Some of the publically available data that I've trimmed (say driver mutatation data from TCGA files, or other gene lists are provided here) your job is to get permission from dbGap to pull in the relevant data and then it should run.... (We'll see)

## STEP 1: Packages to install 
Due to odd dependence with the DESeq2 libraries I've set up two different conda envs stats and rnaseq2 (both available in the landing directory) 

#### STEP 1.1 Create the stats env 
conda env create -f stats.yml

#### STEP 1.2 Create the rnaseq2 env 
conda env create -f rnaseq2.yml 


## STEP 2: Building the genomics figures 
This code is organized using [snakemake](https://snakemake.readthedocs.io/en/stable/). After data is acquired through dbGaP the each of the base figures (before tidying them up using Adobe Illustrator) can be generated one at a time or all together. 

I was a little sloppy having two conda environments to get this off the ground but I found that bioconda's DESeq2 didn't play well with other packages like ComplexHeatmap. What I'm trying to say is sometimes you'll have to switch environments to get the analysis working correctly (All rules that require the rule normalize_RNA) will have to be run with source activate rnaseq2 and stats.   

#### STEP 2.1
source activate stats

#### STEP 2.2 (use -np for a dry-run) 
snakemake -p Figure_1a\
snakemake -p Figure_1b\
snakemake -p Figure_4b\
snakemake -p Figure_4c\
snakemake -p Figure_57b\
snakemake -p Figure_66c

## Building the screening figures 


