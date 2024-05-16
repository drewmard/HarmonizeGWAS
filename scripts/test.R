library(jsonlite)
library(data.table)
library(parallel)
library(dplyr)
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(plyranges)

configFileName = "/oak/stanford/groups/smontgom/amarder/HarmonizeGWAS/config/munge.config"
config = fromJSON(configFileName)

# Initialize
# TMPDIR = "/oak/stanford/groups/smontgom/amarder/tmp"
# HEADDIR = "/oak/stanford/groups/smontgom/amarder/HarmonizeGWAS"
TMPDIR = config$TMPDIR
HEADDIR = config$HEADDIR
AUTOSOMAL_ONLY = TRUE
CONVERT_TO_REF = TRUE
source(paste0(HEADDIR,"/scripts/stringSplitter.R"))

num_cores <- detectCores()
studies = config$studies$study_info

i=12

# initialize:
trait = studies[i]
print(paste0("Running trait: ",trait,"..."))

df = fread(paste0(config$output_base_dir,"hg38","/",trait,"/",trait,".txt.gz"),data.table = F,stringsAsFactors = F)
source(paste0(HEADDIR,"/scripts/convert_to_ref_alt_1kg.R"))

df.lst = readRDS("/oak/stanford/groups/smontgom/amarder/HarmonizeGWAS/out/gwas/munge/hg38/Alzheimers_Bellenguez_2022/tmp.rds")
df.sub = df.lst[[8]]

df.lst = list()
for (chrUse in 22:1) {
  print(paste0("Running chr ",chrUse,"..."))
  df.lst[[chrUse]] = convert_to_ref_alt_1kg(df,chrNum=chrUse)
}
# saveRDS(df.lst,"/oak/stanford/groups/smontgom/amarder/HarmonizeGWAS/out/gwas/munge/hg38/Alzheimers_Bellenguez_2022/tmp.rds")
df = as.data.frame(do.call(rbind,df.lst))
fwrite(df,paste0(config$output_base_dir,"hg38","/",trait,"/",trait,".v2.txt.gz"))


df.sub = subset(df,is.na(ReferenceAllele) & nchar(A0)==1 & nchar(A1)==1)
df.sub[order(df.sub$pvalue),][1:5,]
