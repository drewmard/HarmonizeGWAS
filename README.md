# Harmonize GWAS files!

This repo contains a basic pipeline for downloading and harmonizing GWAS summary statistics into a common format.

```git clone https://github.com/drewmard/HarmonizeGWAS```


## (0) Pre-requisites.

#### (a) Download liftOver chain file (hg19-->hg38).

```cd bin; wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz; cd ..```

#### (b) Install any software necessary.

These scripts use `r`, `liftOver`, and `tabix`.


## (1) Download GWAS.

One way to download GWAS is by using the following script:

```scripts/download_gwas/immune_gwas_download.sh```

First, manually edit `dir=/path/to/HarmonizeGWAS` to your own path in the above script.

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

The config file is used to specify GWAS summary statistics information. See syntax in config/immune.config.

*What do the config arguments mean?*

Global arguments:
```
"input_base_dir" | where gwas files are downloaded to
"output_base_dir" | where harmonized gwas files are moved to
"HEADDIR" | directory to github repo
"TMPDIR" | tmp directory
"hg19ToHg38chain" | liftOver chain file
"liftOver" | path to liftOver software
```

Trait-specific arguments: (each number refers to a specific column in the file; for "study_info", I typically format as Trait_Author_Year)

```
"study_info": "Lupus_Bentham_2015",
"delimiter": "\t",
"rsid_index": "3",
"chr_index": "1",
"snp_pos_index": "2",
"pvalue_index": "6",
"effect_allele_index": "5",
"non_effect_allele_index": "4",
"effect_index": "7",
"se_index": "8",
"direction_index": "7",
"source_build": "hg19"
```


## (3) Harmonize GWAS.

This script will reformat the GWAS file into a common format, and liftOver to hg38 if necessary.

```
config=config/immune.config
Rscript /oak/stanford/groups/smontgom/amarder/HarmonizeGWAS/scripts/munge_gwas_sumstats.R $config
```


## (4) Future to-do.

Code:
* Map SNPs from author-provided effect/non-effect alleles to reference/alternate alleles (I might have the code to do this in CONVERT_TO_REF portion of the script.)
* Map SNPs to RSID. (I don't have a computationally efficient approach for doing this yet.)
* Map SNPs to SPDI. (I pretty much have code for this last part (spdi.R), but this just takes a while to run.) 

Application:
* Does this script work with the following GWAS studies?

```
Psoriasis
https://pubmed.ncbi.nlm.nih.gov/28537254/

Crohn’s/Ulcerative Colitis/IBD
https://pubmed.ncbi.nlm.nih.gov/36038634/
https://pubmed.ncbi.nlm.nih.gov/28067908/

Type 1 Diabetes
https://pubmed.ncbi.nlm.nih.gov/34012112/

Multiple Sclerosis
https://pubmed.ncbi.nlm.nih.gov/31604244/

Allergy
https://pubmed.ncbi.nlm.nih.gov/30013184/
```