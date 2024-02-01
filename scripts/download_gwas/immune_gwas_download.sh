#!/usr/bin/sh

# EDIT MANUALLY:
dir=/oak/stanford/groups/smontgom/amarder/LDSC_pipeline

####################################################################

URL="http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003156/bentham_2015_26502338_sle_efo0002690_1_gwas.sumstats.tsv.gz"
direc=Lupus_Bentham_2015
mkdir -p $dir/gwas/$direc
wget -c $URL --directory $dir/gwas/$direc
mv $dir/gwas/$direc/* $dir/gwas/$direc/$direc.txt.gz

URL="http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132222/GCST90132222_buildGRCh37.tsv.gz"
direc=RheumatoidArthritisTRANS_Ishigaki_2022
mkdir -p $dir/gwas/$direc
wget -c $URL --directory $dir/gwas/$direc
mv $dir/gwas/$direc/* $dir/gwas/$direc/$direc.txt.gz

URL="http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132223/GCST90132223_buildGRCh37.tsv.gz"
direc=RheumatoidArthritisEUR_Ishigaki_2022
mkdir -p $dir/gwas/$direc
wget -c $URL --directory $dir/gwas/$direc
mv $dir/gwas/$direc/* $dir/gwas/$direc/$direc.txt.gz

URL="http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132224/GCST90132224_buildGRCh37.tsv.gz"
direc=RheumatoidArthritisEAS_Ishigaki_2022
mkdir -p $dir/gwas/$direc
wget -c $URL --directory $dir/gwas/$direc
mv $dir/gwas/$direc/* $dir/gwas/$direc/$direc.txt.gz

