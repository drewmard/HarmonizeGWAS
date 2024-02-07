library(BSgenome.Hsapiens.UCSC.hg38)
library(plyranges)
library(data.table)
library(dplyr)

df = fread(paste0(config$output_base_dir,"hg38","/",trait,"/",trait,".txt.gz"),data.table = F,stringsAsFactors = F)

if (AUTOSOMAL_ONLY) {
  df = subset(df,chr %in% c(1:22))
} else {
  df = subset(df,chr %in% c(1:22,"X","Y"))
}

refallele = data.frame(seqnames=paste0("chr",df$chr),start=df$snp_pos,width=1) |> 
  # refallele = data.frame(seqnames=paste0("chr",chrNum),start=tmp$pos,width=1) |> 
  # head(1) |>
  as_granges() |> 
  getSeq(x=Hsapiens) %>% as.character()
tmp2 = data.frame(chr=df$chr,snp_pos=df$snp_pos,ReferenceAllele=refallele)
rm(refallele)

##############

sumstat_expanded = expand_sumstats_via_flip_reverse(df)
colnames(sumstat_expanded)[colnames(sumstat_expanded)=="non_effect_allele"] = "ReferenceAllele"

df2 = merge(sumstat_expanded,tmp2,by=c('chr',"snp_pos","ReferenceAllele"))
df2 = df2[!duplicated(df2),]
colnames(df2)[colnames(df2)=="effect_allele"] = "AltAllele"

# these are SNPs that removed (they dont match the ref genome)
#  subset(df,!(rsid %in% df2$rsid))[1,]
# rsid chr  snp_pos non_effect_allele effect_allele         z pvalue direction
# 60052 1_12891414_G_C   1 12831560                 G             C 0.7893628 0.4299         +
#    subset(df,snp_pos==12831560)
# rsid chr  snp_pos non_effect_allele effect_allele         z pvalue direction
# 60052 1_12891414_G_C   1 12831560                 G             C 0.7893628 0.4299         +
#   60053 1_12891414_T_C   1 12831560                 T             C 0.7735448 0.4392         +
#    subset(df,!(rsid %in% df2$rsid))[2,]
# rsid chr  snp_pos non_effect_allele effect_allele         z    pvalue direction
# 60327 1_12911429_A_T   1 12851576                 A             T 0.5165904 0.6054421         +
#    subset(df,snp_pos==12851576)
# rsid chr  snp_pos non_effect_allele effect_allele          z    pvalue direction
# 60325 1_12911429_C_A   1 12851576                 C             A -0.8806671 0.3784980         -
#   60326 1_12911429_C_T   1 12851576                 C             T  1.8248743 0.0680200         +
#   60327 1_12911429_A_T   1 12851576                 A             T  0.5165904 0.6054421         +

# var = head(df2,10)
# var = subset(df2,ReferenceAllele!="<CN0>" & AlternativeAllele!="<CN0>")
# subset(var,pos==115686847)
# removes structural variants

# var[var$pos %in% df2$pos,][1:5,]

fwrite(df2[,c("chr","snp_pos","ReferenceAllele","AltAllele","rsid")],paste0(TMPDIR,"/pre_var_to_spdi.txt"),quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

cmd = paste0("Rscript ",HEADDIR,"/scripts/var_to_spdi.R ",paste0(TMPDIR,"/pre_var_to_spdi.txt")," ",paste0(TMPDIR,"/post_var_to_spdi.txt"))
system(cmd)

outfile=paste0(TMPDIR,"/post_var_to_spdi.txt")
# write.table(res, file = outfile, row.names = F, col.names = F, quote = F, sep = "\t")
df2.out = fread(outfile,data.table = F,stringsAsFactors = F,header=F)

colnames(df2)[5]= "SPDI"
df2.final = merge(df2,df2.out,by.x="SPDI",by.y="V1") %>%
  select(-SPDI) %>% rename(SPDI=V2)
# df2.final = merge(df2,df2.out,by.x="SPDI",by.y="V1") %>%
#   select(-SPDI) %>% rename(SPDI=V2)
# head(df2.final)
# df2.final %>% rename(SPDI=V2)
out.file = paste0(config$output_base_dir,"hg38","/",trait,"/",trait,".spdi.txt.gz")
fwrite(df2.final,out.file,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T,compress = "gzip")
# fwrite(df2.final,"~/Downloads/Marderstein_Lupus_ClinVar.txt",quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

