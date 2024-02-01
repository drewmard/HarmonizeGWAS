# Harmonize GWAS files!

This repo contains a basic pipeline for downloading and harmonizing GWAS summary statistics into a common format.

```git clone https://github.com/drewmard/HarmonizeGWAS```

## (0) Pre-requisites.

#### (a) Download liftOver chain file (hg19-->hg38).

```cd bin; wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz```

#### (b) Install any packages necessary.

I use a conda/mamba environment called `r`, where `jsonlite` and `data.table` packages are installed in R.

The liftOver script uses a conda/mamba environment called `kent-tools` for running liftOver.

## (1) Download GWAS.

Within the following script, manually edit in the directory where 

```scripts/download_gwas/immune_gwas_download.sh```

## (2)


Rscript /oak/stanford/groups/smontgom/amarder/HarmonizeGWAS/scripts/munge_gwas_sumstats.R $config
