# da_analysis_smk

======================================================

This Snakemake workflow is focused on the analysis of pre-processed whole genome sequencing data in bcf/vcf format of the water flea *Daphnia*.
More information on the [pre-processing steps](https://github.com/tholtzem/daphnia_snakemake_pbs/blob/main/rules/hts.smk) from raw sequencing data to mapped bam files, and on [generating all-sites bcf/vcf from bams](https://github.com/tholtzem/daphnia_snakemake_pbs/blob/main/rules/bcftools.smk) is available in the [daphnia_snakemake_pbs workflow](https://github.com/tholtzem/daphnia_snakemake_pbs).

For more information on [snakemake](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html). 


======================================================

Conda/Mamba and other [dependencies](github)

### Mamba (https://github.com/mamba-org/mamba) is a reimplementation of the conda package manager in C++.

```
# Create environment from yaml file (in envs/):
conda init bash
mamba env create -f envs/s21.yaml

# Activate the environment
conda activate da_analysis

# if you've added new software to install to the conda environment, then you can update:
mamba env update --name da_analysis --file envs/s21.yaml

```

#### python2 env

```
mamba env create -f envs/py2.yaml

```

## Plotting stats and hard-filter all-sites bcf for average site-level using bcftools
here: INFO/DP (2x meanDP or HengLi max depth treshold) and MQ (mapping quality)

see rules/hardfilter_allSites.smk

## Filter genotypes of all-sites bcf using bcftools
### Filter variant and invariant sites separately
see rules/genotypefilter_allSites.smk


## Get SNPs
### with and without mafs
see rules/get_SNPs.smk

## Infer and evaluate admixture proportions using ADMIXTURE and evalAdmix
### using SNPs with maf-filter
see rules/admixture.smk

## Infer structure using PCA implemented in scikit-allel
### using SNPs with maf-filter
note: this was run in a jupyternotebook, not snakemake
see [FastPCA tutorial](http://alimanfoo.github.io/2015/09/28/fast-pca.html)

## Calculate summary stats with pixy
### on filtered all-sites bcf
see rules/pixy.smk

## Infering species and lineage trees with SVDQ in paup
### using SNPs without maf-filter
note: the svdq-step requires some manual adjustments of input files!!!
see rules/svdq.smk
