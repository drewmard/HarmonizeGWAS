library(data.table)

args <- commandArgs(trailingOnly = TRUE)
trait <- args[1]
TMPDIR <- args[2]
path_to_downloaded_gwas <- args[3]

f.bed_hg38 <- paste0(TMPDIR,"/tmp_hg38.bed")
f.gwas_hg19 <- paste0(path_to_downloaded_gwas,"hg19/",trait,"/",trait,".txt.gz")

print("Reading files...")
gwas_hg19 <- fread(f.gwas_hg19,data.table = F,stringsAsFactors = F)
bed_hg38 <- fread(f.bed_hg38,data.table = F,stringsAsFactors = F,header = F)

print("Matching...")
gwas_hg38 <- gwas_hg19
i=match(gwas_hg38$rsid,bed_hg38$V4)
gwas_hg38$chr <- substring(bed_hg38$V1[i],4)
gwas_hg38$snp_pos <- bed_hg38$V3[i]

print("Removing unmatched positions...")
gwas_hg38 <- gwas_hg38[!is.na(gwas_hg38$snp_pos),]

print("Reordering incase hg19-to-hg38 conversion swapped the order...")
gwas_hg38 <- gwas_hg38[order(gwas_hg38$chr,gwas_hg38$snp_pos),]

print("Saving hg38 file...")
system(paste0('mkdir -p ',path_to_downloaded_gwas,"hg38/",trait))
f.gwas_hg38 <- paste0(path_to_downloaded_gwas,"hg38/",trait,"/",trait,".txt")
print(paste0('file: ',f.gwas_hg38))
fwrite(x=gwas_hg38,file=f.gwas_hg38,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)



