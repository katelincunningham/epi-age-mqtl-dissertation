# Age and Genetic Impact on Age-associated CpGs

## Project Overview
These analysis were part of submission for MSc in Bioinformatics at Queen Mary University of London.

This repository contains code and scripts used to analyse DNA methylation at the CpG site cg10501210, with a focus on its association with ageing across multiple tissues and the potential genetic regulation by the mQTL rs41317993.

The analysis integrates array-based methylation datasets and Oxford Nanopore long-read sequencing data to:

* Test whether cg10501210 is consistently hypomethylated with age across tissues.

* Investigate age and genotype effects on cg10501210 methylation in long-read datasets.

* Explore haplotype-specific methylation patterns.

## Script Usage

The `DNAm_array_data_alaysis` folder conatins scripts used to investiagte the pan-tissue effect of cg10501210 DNA methylation.

The `Age_DNAm_analysis_LRS` folder included scripts used to evaluate effect of age on cg10501210 DNA methylation in ONT LRS data.

The `methyl_AGT_plot.sh` script along with the required R scripts `00_Rlibs.R` , `01_Combine_bed.R` , `02_Combine_Data.R` and `03__LRS_plots.R` can be used to visualise the genotypic effect of an mQTL on DNA methylation at a CpG site.

It requres:
* BAM and phased VCF files
* BED file of CpG location
* Chromosome location of mQTL



