# module load R

conda init bash
source ~/.bashrc

trait=$1
TMPDIR=$2
HEADDIR=$3
path_to_downloaded_gwas=$4

# hgBuild=`awk -F, -v pat=$trait '$1~pat {print $7}' /oak/stanford/groups/smontgom/amarder/neuro-variants/scripts/snps/helper_func/colocal.csv`

# if [[ "$hgBuild" == *19* ]]; 
# then

echo "Running for $trait ..."

mkdir -p $TMPDIR

# path_to_downloaded_gwas=/oak/stanford/groups/smontgom/amarder/LDSC_pipeline/gwas/munge
mkdir -p $path_to_downloaded_gwas/hg38/$trait

echo "Script: Performing hg19 to hg38 liftOver for sumstats..."

# reformat:
echo "Reformatting..."
sumstatpath=$path_to_downloaded_gwas/hg19/$trait/$trait.txt.gz
zcat $sumstatpath | awk -f $HEADDIR/t.awk -v cols=chr,snp_pos,rsid | awk -v OFS='\t' '{print "chr"$1,$2-1,$2,$3}' > $TMPDIR/tmp_hg19.bed

# liftover: hg19 to hg38 ################################################
echo "Running liftOver..."
bedFile=$TMPDIR/tmp_hg19.bed
hg19ToHg38chain=/oak/stanford/groups/smontgom/amarder/bin/chain/hg19ToHg38.over.chain.gz
source /oak/stanford/groups/smontgom/amarder/bin/conda_init.sh
conda activate kent-tools
liftOver $bedFile $hg19ToHg38chain $TMPDIR/tmp_hg38.bed $TMPDIR/tmp_unmapp.bed
conda deactivate 

# echo "If not converted in liftOver, then looking for dbsnp file..."
# # unmappedRs="rs6541019|rs7543597"
# # ugly scrappy combo of sed/awk code, courtesy of Rachel & Andrew's sheer geniusness
# unmappedRs=`grep "^[^#]" $TMPDIR/tmp_unmapp.bed | head | cut -f4  | sed -z -r 's/\n/|/g;' | awk -F'|' 'OFS="|" {printf $1; for(i=2;i<NF;i++){printf "|"$i}; print ""}'`
# dbsnphg38=/oak/stanford/groups/smontgom/mgloud/projects/gwas-download/gwas-download/munge/dbsnp/sorted_1kg_matched_hg38_snp150.txt.gz
# zcat $dbsnphg38 | grep -E "$unmappedRs$"
# 
# grep "^[^#]" $TMPDIR/tmp_unmapp.bed | head | cut -f4 > $TMPDIR/rs_unmapp.bed
# grep -F -f $unmappedRs $(zcat $dbsnphg38)

echo "Preparing to write hg38 GWAS sumstat file..."
conda activate r
Rscript $HEADDIR/scripts/create_hg38_sumstat.R $trait $TMPDIR $path_to_downloaded_gwas
conda deactivate

echo "Tabix the output file..."
conda activate tabix
bgzip -f $path_to_downloaded_gwas/hg38/$trait/$trait.txt
tabix -f -s 2 -b 3 -e 3 -S 1 $path_to_downloaded_gwas/hg38/$trait/$trait.txt.gz
conda deactivate
echo "hg38 transformation complete."
echo "File saved to $path_to_downloaded_gwas/hg38/$trait/$trait.txt.gz"

# fi
