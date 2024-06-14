#
# GWAS Munge Pipeline!
# this is a modified pipeline from Mike Gloudesman
# the issue with Mike's python pipeline was that if SNPs gained a new rsid but were retained between hg19 --> hg38
# they would be dropped because the rsid would be different
# however, liftover should be done on the coordinates
# and then new rsid's assigned!

# srun --account=default --partition=interactive --time=24:00:00 --mem=64G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash 
# conda activate /oak/stanford/groups/smontgom/amarder/micromamba/envs/HarmonizeGWAS

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
configFileName = "/oak/stanford/groups/smontgom/amarder/HarmonizeGWAS/config/munge2.config"
config = fromJSON(configFileName)

# Initialize
# TMPDIR = "/oak/stanford/groups/smontgom/amarder/tmp"
# HEADDIR = "/oak/stanford/groups/smontgom/amarder/HarmonizeGWAS"
TMPDIR = config$TMPDIR
HEADDIR = config$HEADDIR
AUTOSOMAL_ONLY = TRUE
CONVERT_TO_REF = TRUE
# source(paste0(HEADDIR,"/scripts/stringSplitter.R"))

num_cores <- detectCores()
studies = config$studies$study_info
for (i in 11:length(studies)) { # need to fix 7:Major_Depression_Howard_2019... and 10: BMI_Yengo_2018
  
  # initialize:
  trait = studies[i]
  print(paste0("Running trait",i,"/",length(studies),": ",trait,"..."))
  
  # Pull out columns specified (not NULL):
  study_info = subset(config$studies,study_info==trait)
  
  if (is.na(study_info$rsid_index) | study_info$rsid_index==-1) {
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
                    "freq_index",
                    "direction_index")
  cols_to_use_ordered = desired_order[desired_order %in% cols_to_use]
  idx = as.numeric(unlist(study_info[,cols_to_use_ordered]))
  
  # Read in sum stats:
  f = paste0(config$input_base_dir,trait,"/",trait,".txt.gz")
  # f = paste0(config$input_base_dir,"/",trait,"/",trait,".txt.gz")
  # df = fread(f,data.table = F,stringsAsFactors = F,fill = TRUE,sep = study_info$delimiter)
  # df = fread(f,data.table = F,stringsAsFactors = F,sep = study_info$delimiter)
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
  if (!("rsid" %in% colnames(df))) {df[,"rsid"] = paste0("snp",1:nrow(df))}
  df = df %>% select(rsid,everything())
  
  # Specify build:
  build = ifelse(is.na(study_info$source_build),config$genome_build,study_info$source_build)
  
  # Save dataframe to the specified source build (e.g. hg19 or hg38):
  dir.create(paste0(config$output_base_dir),showWarnings = FALSE)
  dir.create(paste0(config$output_base_dir,build),showWarnings = FALSE)
  dir.create(paste0(config$output_base_dir,build,"/",trait),showWarnings = FALSE)
  f.out = paste0(config$output_base_dir,build,"/",trait,"/",trait,".txt.gz")
  options(scipen = 999)
  fwrite(df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T,compress = "gzip")
  options(scipen = 0)
  print(paste0("Data frame saved to: ",f.out," ..."))
  
  # dbsnp = fread("/oak/stanford/groups/smontgom/mgloud/projects/gwas-download/gwas-download/munge/dbsnp/sorted_1kg_matched_hg38_snp150.txt.gz",data.table = F,stringsAsFactors = F)
  # y=paste(dbsnp$V1,dbsnp$V2)
  # dbsnp = dbsnp[!duplicated(y),]
  # dbsnp$V1 = substring(dbsnp$V1,4)
  # fwrite(dbsnp,"/oak/stanford/groups/smontgom/amarder/data/dbsnp/hg38/sorted_1kg_matched_hg38_snp150.no_dup.txt.gz",compress = "gzip",quote = F,sep = '\t',row.names = F,col.names = F,na = "NA")
  
  # # Function to process a batch of SNPs
  # process_batch <- function(df_batch, tbx) {
  #   param <- GRanges(df_batch$chr, IRanges(df_batch$snp_pos - 1, df_batch$snp_pos))
  #   result <- scanTabix(tbx, param=param)
  #   
  #   match_data <- lapply(seq_along(result), function(i) {
  #     matches <- result[[i]]
  #     if (length(matches) > 0) {
  #       match_list <- strsplit(matches, "\t")
  #       do.call(rbind, lapply(match_list, function(match) {
  #         data.frame(
  #           chr = match[1],
  #           pos = as.integer(match[2]),
  #           a1 = match[4],
  #           a2 = match[5],
  #           rsid = match[3],
  #           stringsAsFactors = FALSE
  #         )
  #       }))
  #     } else {
  #       NULL
  #     }
  #   })
  #   do.call(rbind, match_data)
  # }
  # 
  # library(Rsamtools)
  # # Specify the path to the Tabix file
  # tabix_file <- paste0("/oak/stanford/groups/akundaje/soumyak/refs/dbsnp_hg38/","hg38",".dbsnp156.gz")
  # open(tbx)
  # 
  # # Split the data frame into batches
  # batch_size <- 100
  # df_batches <- split(df[1:5000,], as.factor(ceiling(seq_len(1000) / batch_size)))
  # match_data_list <- lapply(df_batches, process_batch, tbx = tbx)
  # match_data <- do.call(rbind, match_data_list)
  # 
  # # Close the Tabix file
  # close(tbx)
  
  # If build == hg19, also save an hg38 version:
  if (build=="hg19") {
    # cmd = paste0(HEADDIR,"/scripts/liftover_hg19_to_hg38.sh ",trait," ",TMPDIR," ",HEADDIR," ",config$output_base_dir)
    cmd = paste0(HEADDIR,"/scripts/liftover_hg19_to_hg38.sh ",trait," ",TMPDIR," ",HEADDIR," ",config$output_base_dir," ",config$hg19ToHg38chain," ",config$liftOver)
    system(cmd)
  }
  
  # Read in hg38 dataset:
  df = fread(paste0(config$output_base_dir,"hg38","/",trait,"/",trait,".txt.gz"),data.table = F,stringsAsFactors = F)
  
  if (!("rsid_index" %in% cols_to_use)) {
    print("Crudely adding RSIDs... (not matching on alleles because it takes too long...)")
    dbsnp = fread("/oak/stanford/groups/smontgom/amarder/data/dbsnp/hg38/sorted_1kg_matched_hg38_snp150.no_dup.txt.gz",data.table = F,stringsAsFactors = F)
    df = merge(df[,colnames(df)!="V3"],dbsnp,by=c("chr","snp_pos"),by.y=c("V1","V2"),all.x = TRUE)
    colnames(df)[colnames(df)=="V3"] = "rsid"
  }
  
  if (AUTOSOMAL_ONLY) {
    df = subset(df,chr %in% c(1:22))
  } else {
    df = subset(df,chr %in% c(1:22,"X","Y"))
  }
  
  if (CONVERT_TO_REF) {
    
    source(paste0(HEADDIR,"/scripts/convert_to_ref_alt_1kg.R"))
    
    df.lst = list()
    # for (chrUse in 22:1) {
    for (chrUse in 22:1) {
      print(paste0("Running chr ",chrUse,"..."))
      attempt = 0
      repeat {
        # attempt = attempt + 1
        result = convert_to_ref_alt_1kg(df,chrNum=chrUse,parallel=FALSE)
        # if (attempt >= 2) {result = convert_to_ref_alt_1kg(df,chrNum=chrUse,parallel=FALSE)}
        # else {
        #   result <- tryCatch({
        #     convert_to_ref_alt_1kg(df,chrNum=chrUse,parallel=TRUE)
        #   }, error = function(e) {
        #     NULL
        #   })
        # }
        if (!is.null(result)) break
      }
      df.lst[[chrUse]] = result
    }
    df = as.data.frame(do.call(rbind,df.lst))
    f.out = paste0(config$output_base_dir,"hg38","/",trait,"/",trait,".v2.txt.gz")
    fwrite(df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T,compress = "gzip")
  } else {
    df$snp_id = unlist(
      mclapply(
        1:nrow(df),
        function(i) paste(
          paste0("chr",df$chr[i]),df$snp_pos[i],df$non_effect_allele[i],df$effect_allele[i],sep = "_"
        ),
        mc.cores = 8
      )
    )
  }
}

