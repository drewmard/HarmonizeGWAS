# Harmonize GWAS files!

This repo contains a basic pipeline for downloading and harmonizing GWAS summary statistics into a common format.

```git clone https://github.com/drewmard/HarmonizeGWAS```

## (0) Pre-requisites.

#### Download liftOver chain file (hg19-->hg38).

```cd bin; wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz```

## (1) Download GWAS.

## (2)
/oak/stanford/groups/smontgom/amarder/HarmonizeGWAS/scripts/munge_gwas_sumstats.R $config
