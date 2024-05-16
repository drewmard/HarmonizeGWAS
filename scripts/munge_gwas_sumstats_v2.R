# GWAS Munge Pipeline!
# this is a modified pipeline from Mike Gloudesman
# the issue with Mike's python pipeline was that if SNPs gained a new rsid but were retained between hg19 --> hg38
# they would be dropped because the rsid would be different
# however, liftover should be done on the coordinates
# and then new rsid's assigned!

library(jsonlite)
library(data.table)
library(parallel)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(plyranges)


# Config path:
args = commandArgs(trailingOnly=TRUE)
configFileName = args[1]
# configFileName = "/oak/stanford/groups/smontgom/amarder/HarmonizeGWAS/config/immune.config"
# configFileName = "/oak/stanford/groups/smontgom/amarder/HarmonizeGWAS/config/munge.config"
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
for (i in 1:length(studies)) {
  
  # initialize:
  trait = studies[i]
  print(paste0("Running trait: ",trait,"..."))
  
  # Pull out columns specified (not NULL):
  study_info = subset(config$studies,study_info==trait)
  
  if (study_info$rsid_index==-1) {
    print("Note: no rsid info in the file!")
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
  
  # Reformat:
  df = df[,idx]
  colnames(df) = sub("_index","",cols_to_use_ordered) # rename! (remove _index)
  df$chr = sub("chr","",df$chr) # remove chr!
  
  # add direction column
  if ("effect" %in% colnames(df)) {
    df$direction = ifelse(sign(df$effect)>0,"+","-")
  } else if ("zscore" %in% colnames(df)) {
    df$direction = ifelse(sign(df$zscore)>0,"+","-")
  } else if ("or" %in% colnames(df)) {
    df$direction = ifelse(sign(df$or)>1,"+","-")
  } else {
    stop("Missing effect column!")
  }
  colnames(df)[colnames(df)=="effect"] = "beta" # rename if necessary!
  colnames(df)[colnames(df)=="zscore"] = "z" # rename if necessary!
  if ("non_effect_allele" %in% colnames(df)) {df[,"non_effect_allele"] = toupper(df[,"non_effect_allele"])}
  if ("effect_allele" %in% colnames(df)) {df[,"effect_allele"] = toupper(df[,"effect_allele"])}
  print(paste0("Summary statistics reformatted..."))
  
  # # Add rsid if needed:
  # build = ifelse(is.na(study_info$source_build),config$genome_build,study_info$source_build)
  # if (!("rsid" %in% colnames(df))) {
  #   # Specify rsid as chr, snp_pos, non_effect_allele, effect_allele
  #   df$rsid = unlist(mclapply(mc.cores=8,1:nrow(df),function(i) {
  #     # if(i%%10000 == 0){print(i)} 
  #     return(paste(df[i,c("chr","snp_pos","non_effect_allele","effect_allele")],collapse = '_'))
  #   }))
  #   df = df[,c("rsid",colnames(df)[colnames(df)!="rsid"])]
  # }
  # 
  # # Search chr and pos if needed
  # if (!("chr" %in% colnames(df))) {
  #   test=1 # not yet implemented
  #   # df = dbsnpQuery(data_input=df,trait=trait,rsid_col="rsid",tmpdir=TMPDIR,a1_col="non_effect_allele",a2_col="effect_allele")
  #   # df = df[,c("rsid",colnames(df)[colnames(df)!="rsid"])]
  # }
  
  
  # Save dataframe to the specified source build (e.g. hg19 or hg38):
  dir.create(paste0(config$output_base_dir),showWarnings = FALSE)
  dir.create(paste0(config$output_base_dir,build),showWarnings = FALSE)
  dir.create(paste0(config$output_base_dir,build,"/",trait),showWarnings = FALSE)
  f.out = paste0(config$output_base_dir,build,"/",trait,"/",trait,".txt.gz")
  fwrite(df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T,compress = "gzip")
  print(paste0("Data frame saved to: ",f.out," ..."))
  
  # If build == hg19, also save an hg38 version:
  if (build=="hg19") {
    # cmd = paste0(HEADDIR,"/scripts/liftover_hg19_to_hg38.sh ",trait," ",TMPDIR," ",HEADDIR," ",config$output_base_dir)
    cmd = paste0(HEADDIR,"/scripts/liftover_hg19_to_hg38.sh ",trait," ",TMPDIR," ",HEADDIR," ",config$output_base_dir," ",config$hg19ToHg38chain," ",config$liftOver)
    system(cmd)
  }

  # if (AUTOSOMAL_ONLY) {
  #   df = subset(df,chr %in% c(1:22))
  # } else {
  #   df = subset(df,chr %in% c(1:22,"X","Y"))
  # }
  
  if (CONVERT_TO_REF) {
    
    df = fread(paste0(config$output_base_dir,"hg38","/",trait,"/",trait,".txt.gz"),data.table = F,stringsAsFactors = F)
    source(paste0(HEADDIR,"/scripts/convert_to_ref_alt_1kg.R"))
    
    df.lst = list()
    for (chrUse in 22:1) {
      print(paste0("Running chr ",chrUse,"..."))
      df.lst[[chrUse]] = convert_to_ref_alt_1kg(df,chrNum=chrUse)
    }
    # saveRDS(df.lst,"/oak/stanford/groups/smontgom/amarder/HarmonizeGWAS/out/gwas/munge/hg38/Alzheimers_Bellenguez_2022/tmp.rds")
    df = as.data.frame(do.call(rbind,df.lst))
    f.out = paste0(config$output_base_dir,"hg38","/",trait,"/",trait,".v2.txt.gz")
    fwrite(df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T,compress = "gzip")
  } else {
    df$snp_id = unlist(
      mclapply(
        1:nrow(df),
        function(i) paste(
          df$chr[i],df$snp_pos[i],df$non_effect_allele[i],df$effect_allele[i],sep = "_"
        ),
        mc.cores = 8
      )
    )
  }
}

