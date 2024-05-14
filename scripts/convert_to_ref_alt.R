convert_to_ref_alt = function(df) {
  
  # This function aligns GWAS sumstat file to the reference genome.
  # Therefore, ref/alt alleles are used instead of effect allele notation.
  # Betas are convert to refer to the alt allele.
  
  # Preserve order:
  df$i = 1:nrow(df)
  
  # SNP identifier - will be overwritten in the end:
  num_cores <- detectCores()
  df$snp_id = unlist(
    mclapply(
      1:nrow(df),
      function(i) paste(
        df$chr[i],df$snp_pos[i],df$non_effect_allele[i],df$effect_allele[i],sep = "_"
      ),
      mc.cores = num_cores - 1
    )
  )
  
  # Fetch reference alleles:
  refallele = data.frame(seqnames=paste0("chr",df$chr),start=df$snp_pos,width=1) |> 
    as_granges() |> 
    getSeq(x=Hsapiens) %>% as.character()
  # consider modifying to submit X SNPs at a time:
  # ???????????
  tmp2 = data.frame(chr=df$chr,snp_pos=df$snp_pos,ReferenceAllele=refallele)

  # Flip & reverse sumstats, and then merge w/ ref allele matching:
  source(paste0(HEADDIR,"/scripts/expand_sumstats_via_flip_reverse.R"))
  idx = nchar(df$non_effect_allele)==1 & nchar(df$effect_allele)==1
  sumstat_expanded = expand_sumstats_via_flip_reverse(df[idx,])
  
  # Merge with reference alleles:
  tmp2 = merge(sumstat_expanded,tmp2,by.x=c("chr","snp_pos","non_effect_allele"),by.y=c('chr',"snp_pos","ReferenceAllele"))
  
  # Remove duplicates (from expand_sumstats_via_flip_reverse() function):
  tmp2 = tmp2[!duplicated(tmp2$snp_id),]
  rm(sumstat_expanded)
  
  # Append indels:
  indel = !idx
  if (sum(indel) > 0) {
    df = rbind(tmp2,df[indel,])
  }
  rm(tmp2)
  
  # Reformat:
  colnames(df)[colnames(df)=="non_effect_allele"] = "ReferenceAllele"
  colnames(df)[colnames(df)=="effect_allele"] = "AltAllele"
  df$snp_id = unlist(
    mclapply(
      1:nrow(df),
      function(i) paste(
        df$chr[i],df$snp_pos[i],df$ReferenceAllele[i],df$AltAllele[i],sep = "_"
      ),
      mc.cores = num_cores
    )
  )
  df = df[order(df$i),]
  
  return(df)
}

