# module load R
# conda init bash
# source ~/.bashrc

trait=$1
TMPDIR=$2
HEADDIR=$3 # HarmonizeGWAS directory
hg19ToHg38chain=$4
liftOver=$5

# hg19ToHg38chain=/oak/stanford/groups/smontgom/amarder/bin/chain/hg19ToHg38.over.chain.gz
# liftOver=/home/amarder/micromamba/envs/kent-tools/bin/liftOver
path_to_downloaded_gwas=$HEADDIR/out/gwas/munge

echo "Running for $trait ..."

mkdir -p $TMPDIR
mkdir -p $path_to_downloaded_gwas/hg38/$trait

echo "Script: Performing hg19 to hg38 liftOver for sumstats..."

# reformat:
echo "Reformatting..."
sumstatpath=$path_to_downloaded_gwas/hg19/$trait/$trait.txt.gz
zcat $sumstatpath | awk -f $HEADDIR/scripts/t.awk -v cols=chr,snp_pos,rsid | awk -v OFS='\t' '{print "chr"$1,$2-1,$2,$3}' > $TMPDIR/tmp_hg19.bed

# liftover: hg19 to hg38 ################################################
echo "Running liftOver..."
bedFile=$TMPDIR/tmp_hg19.bed
$liftOver $bedFile $hg19ToHg38chain $TMPDIR/tmp_hg38.bed $TMPDIR/tmp_unmapp.bed

echo "Preparing to write hg38 GWAS sumstat file..."
Rscript $HEADDIR/scripts/create_hg38_sumstat.R $trait $TMPDIR $path_to_downloaded_gwas

echo "Tabix the output file..."
rm $path_to_downloaded_gwas/hg38/$trait/$trait.txt.gz
rm $path_to_downloaded_gwas/hg38/$trait/$trait.txt.gz.tbi
bgzip -f $path_to_downloaded_gwas/hg38/$trait/$trait.txt
tabix -f -s 2 -b 3 -e 3 -S 1 $path_to_downloaded_gwas/hg38/$trait/$trait.txt.gz
echo "hg38 transformation complete."
echo "File saved to $path_to_downloaded_gwas/hg38/$trait/$trait.txt.gz"

# fi
