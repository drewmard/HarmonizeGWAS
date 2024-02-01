# GWAS Munge Pipeline!
# this is a modified pipeline from Mike Gloudesman
# the issue with Mike's python pipeline was that if SNPs gained a new rsid but were retained between hg19 --> hg38
# they would be dropped because the rsid would be different
# however, liftover should be done on the coordinates
# and then new rsid's assigned!

library(jsonlite)
library(data.table)

# Config path:
# config <- fromJSON("~/Documents/Research/munge_test.config")
# config <- fromJSON("~/Downloads/munge.config")
# config = fromJSON("/oak/stanford/groups/smontgom/amarder/LDSC_pipeline/config/munge_test.config")
args = commandArgs(trailingOnly=TRUE)
configFileName = args[1]
configFileName = "/oak/stanford/groups/smontgom/amarder/HarmonizeGWAS/config/munge_test.config"
config = fromJSON(configFileName)

studies = config$studies$study_info
for (i in 1:length(studies)) {
  
  # initialize:
  trait = studies[i]
  print(paste0("Running trait: ",trait,"..."))
  
  # Pull out columns specified (not NULL):
  study_info = subset(config$studies,study_info==trait)
  
  if (study_info$rsid_index==-1) {
    print("Note: no rsid info in the file!")
    # rsid = identify_rsid(chrom,snp_pos,effect_allele_index,non_effect_allele_index)
  }
  
  tmp = study_info[,grep("index",colnames(study_info))]
  cols_to_use = colnames(tmp)[!is.na(tmp) & tmp!=-1]
  
  # Need to do: Switch this such that create a pvalue column from log10 p-value column!
  if ("pvalue_index" %in% cols_to_use) {
    cols_to_use = cols_to_use[!(cols_to_use %in% c("log_pvalue_index","neg_log_pvalue_index"))]
  } else if ("log_pvalue_index" %in% cols_to_use)  {
    cols_to_use = cols_to_use[!(cols_to_use %in% c("neg_log_pvalue_index"))]
  } else if ("neg_log_pvalue_index" %in% cols_to_use)  {
    cols_to_use = cols_to_use
  }

  # Rearrange dataframe in desired order:
  desired_order = c("rsid_index","chr_index","snp_pos_index",
                    "non_effect_allele_index","effect_allele_index",
                    "effect_index","or_index","se_index","zscore_index",
                    "pvalue_index","log_pvalue_index","neg_log_pvalue_index",
                    "direction_index")
  cols_to_use_ordered = desired_order[desired_order %in% cols_to_use]
  idx = as.numeric(unlist(study_info[,cols_to_use_ordered]))
  
  # Read in sum stats:
  f = paste0(config$input_base_dir,"/",trait,"/",trait,".txt.gz")
  df = fread(f,data.table = F,stringsAsFactors = F)
  print(paste0("Summary statistics read into memory..."))
  df = df[,idx]
  colnames(df) = sub("_index","",cols_to_use_ordered) # rename! (remove _index)
  df$chr = sub("chr","",df$chr) # remove chr!
  df$direction = ifelse(sign(df$effect)>0,"+","-")
  colnames(df)[colnames(df)=="effect"] = "beta" # rename if necessary!
  print(paste0("Summary statistics reformatted..."))
  
  # Save dataframe to the specified source build (e.g. hg19 or hg38):
  dir.create(paste0(config$output_base_dir,study_info$source_build,"/",trait),showWarnings = FALSE)
  f.out = paste0(config$output_base_dir,study_info$source_build,"/",trait,"/",trait,".txt.gz")
  fwrite(df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T,compress = "gzip")
  print(paste0("Data frame saved to: ",f.out," ..."))
  
  if (study_info$source_build=="hg19") {
    
  }
}

