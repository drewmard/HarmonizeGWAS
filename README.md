# Harmonize GWAS files!

This repo contains a basic pipeline for downloading and harmonizing GWAS summary statistics into a common format.

```git clone https://github.com/drewmard/HarmonizeGWAS```


## (0) Pre-requisites.

#### (a) Download liftOver chain file (hg19-->hg38).

```cd bin; wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz; cd ..```

#### (b) Install any packages necessary.

I use a conda/mamba environment called `r`, where `jsonlite` and `data.table` packages are installed in R.

I use a conda/mamba environment called `kent-tools` for running `liftOver`.

And finally, I use a conda/mamba environment called `tabix` for running `tabix`.

#### (c) Download and pre-process dbSNP files.

Note: this is a large download (12G x 2).

```
conda activate tabix
mkdir -p bin/dbsnp/hg19; cd bin/dbsnp/hg19
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp151.txt.gz
zcat snp151.txt.gz | cut -f2,3,4,5,10 | bgzip -c > snp151.v2.txt.gz
tabix -p bed snp151.v2.txt.gz
cd ../../..
mkdir -p bin/dbsnp/hg38; cd bin/dbsnp/hg38
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp151.txt.gz
zcat snp151.txt.gz | cut -f2,3,4,5,10 | bgzip -c > snp151.v2.txt.gz
tabix -p bed snp151.v2.txt.gz
cd ../../..
```

or working with only common variants:

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp151Common.txt.gz
zcat snp151Common.txt.gz | cut -f2,3,4,5,10 | bgzip -c > snp151Common.v2.txt.gz


## (1) Download GWAS.

One way to download GWAS to use the following script:

```scripts/download_gwas/immune_gwas_download.sh```

First, manually edit `dir=/path/to/HarmonizeGWAS` to your own path.

Next, files can *usually* be downloaded by specifying `URL` and `direc` and running the following commands. This will download the file from `URL` (e.g. the GWAS catalog), and save using the `direc` information. For example:

```
URL="http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003156/bentham_2015_26502338_sle_efo0002690_1_gwas.sumstats.tsv.gz"
direc=Lupus_Bentham_2015
mkdir -p $dir/gwas/$direc
wget -c $URL --directory $dir/gwas/$direc
mv $dir/gwas/$direc/* $dir/gwas/$direc/$direc.txt.gz
```

This will download the SLE GWAS by Bentham et al. 2015 from the GWAS catalog (the specified URL) to `gwas/Lupus_Bentham_2015/Lupus_Bentham_2015.txt.gz`. However, for some files (particularly those not in GWAS catalog), there may be additional preprocessing necessary if the file is in a weird format in order to get it to be read into R.





## (2) Modify the config file.

The config file is used to specify GWAS summary statistics information.



## (3) Harmonize GWAS.

This script will reformat the GWAS file into a common format.

```
config=config/immune.config
Rscript /oak/stanford/groups/smontgom/amarder/HarmonizeGWAS/scripts/munge_gwas_sumstats.R $config
```