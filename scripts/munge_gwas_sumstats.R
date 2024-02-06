# GWAS Munge Pipeline!
# this is a modified pipeline from Mike Gloudesman
# the issue with Mike's python pipeline was that if SNPs gained a new rsid but were retained between hg19 --> hg38
# they would be dropped because the rsid would be different
# however, liftover should be done on the coordinates
# and then new rsid's assigned!

library(jsonlite)
library(data.table)

# Config path:
args = commandArgs(trailingOnly=TRUE)
configFileName = args[1]
configFileName = "/oak/stanford/groups/smontgom/amarder/HarmonizeGWAS/config/immune.config"
configFileName = "/oak/stanford/groups/smontgom/amarder/HarmonizeGWAS/config/munge.config"
config = fromJSON(configFileName)

# Initialize
path_to_dbsnp = "/oak/stanford/groups/smontgom/amarder/data/dbsnp" # this could also be: /oak/stanford/groups/smontgom/amarder/HarmonizeGWAS/bin/dbsnp/hg19
TMPDIR = "/oak/stanford/groups/smontgom/amarder/tmp"
HEADDIR = "/oak/stanford/groups/smontgom/amarder/HarmonizeGWAS"
path_to_dbsnp = config$path_to_dbsnp # this could also be: /oak/stanford/groups/smontgom/amarder/HarmonizeGWAS/bin/dbsnp/hg19
TMPDIR = config$TMPDIR
HEADDIR = config$HEADDIR

# Source function:
source(paste0(HEADDIR,"/scripts/dbsnpQuery.R"))

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
  print(paste0("Summary statistics reformatted..."))
  
  # Add rsid if needed:
  build = ifelse(is.na(study_info$source_build),config$genome_build,study_info$source_build)
  if (!("rsid" %in% colnames(df))) {
    df = dbsnpQuery(data_input=df,
                    type="chr_pos",
                    tmpdir=TMPDIR,
                    trait=trait,
                    SNPFILE=paste0(path_to_dbsnp,"/",build,"/snp151Common.v2.txt.gz")
                    )
  }
  
  # Save dataframe to the specified source build (e.g. hg19 or hg38):
  dir.create(paste0(config$output_base_dir),showWarnings = FALSE)
  dir.create(paste0(config$output_base_dir,build),showWarnings = FALSE)
  dir.create(paste0(config$output_base_dir,build,"/",trait),showWarnings = FALSE)
  f.out = paste0(config$output_base_dir,build,"/",trait,"/",trait,".txt.gz")
  fwrite(df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T,compress = "gzip")
  print(paste0("Data frame saved to: ",f.out," ..."))
  
  # If build == hg19, also save an hg38 version:
  if (build=="hg19") {
    cmd = paste0(HEADDIR,"/scripts/liftover_hg19_to_hg38.sh ",trait," ",TMPDIR," ",HEADDIR," ",config$output_base_dir)
    system(cmd)
  }
}

