# Harmonize GWAS files!

This repo contains a basic pipeline for downloading and harmonizing GWAS summary statistics into a common format.

```git clone https://github.com/drewmard/HarmonizeGWAS```

## (0) Pre-requisites.

#### (a) Download liftOver chain file (hg19-->hg38).

```cd bin; wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz```

#### (b) Install any packages necessary.

I use a conda/mamba environment called `r`, where `jsonlite` and `data.table` packages are installed in R.

I use a conda/mamba environment called `kent-tools` for running liftOver.

## (1) Download GWAS.

Within the following script:

```scripts/download_gwas/immune_gwas_download.sh```

manually edit `dir=/path/to/HarmonizeGWAS` to your own path.

Next, files can be downloaded by specifying `URL` and `direc`. This will download the file from `URL` (e.g. the GWAS catalog), and save using the `direc` information. For example:

```
URL="http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003156/bentham_2015_26502338_sle_efo0002690_1_gwas.sumstats.tsv.gz"
direc=Lupus_Bentham_2015
mkdir -p $dir/gwas/$direc
wget -c $URL --directory $dir/gwas/$direc
mv $dir/gwas/$direc/* $dir/gwas/$direc/$direc.txt.gz
```

This will download the SLE GWAS by Bentham et al. 2015 from the GWAS catalog (the specified URL) to `gwas/Lupus_Bentham_2015/Lupus_Bentham_2015.txt.gz`.

## (2) Harmonize GWAS.


```
config=config/immune.config
Rscript /oak/stanford/groups/smontgom/amarder/HarmonizeGWAS/scripts/munge_gwas_sumstats.R $config
```