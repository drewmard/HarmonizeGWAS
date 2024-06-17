convert_to_ref_alt_1kg = function(df,chrNum,parallel=FALSE) {
  
  # This function aligns GWAS sumstat file to the reference genome.
  # Therefore, ref/alt alleles are used instead of effect allele notation.
  # Betas are convert to refer to the alt allele.
  
  # Fetch reference alleles:
  # f.geno=paste0("/oak/stanford/groups/smontgom/amarder/bin/1kg/30x/plink_indel/1kg.all.hg38_30x.chr",chrNum,".bim")
  f.geno=paste0("/oak/stanford/groups/akundaje/soumyak/refs/plink_eur_1kg_30x_hg38_snvs_indels/chr",chrNum,".eur.filtered.bim")
  geno <- fread(f.geno,data.table = F,stringsAsFactors = F)#, verbose = FALSE)
  colnames(geno)[c(1,4,6,5)] = c("chr","pos","a0","a1")
  
  # Subset gwas
  sumstats = subset(df,chr==chrNum)
  print(paste0("Processing ",nrow(sumstats)," variants..."))
  
  # Reformat:
  sumstats$chr = as.numeric(sumstats$chr)
  colnames(sumstats)[colnames(sumstats)=="snp_pos"] = "pos"
  sumstats$A0 = sumstats$non_effect_allele
  sumstats$A1 = sumstats$effect_allele
  
  # Preserve order:
  sumstats$iden = 1:nrow(sumstats)
  
  # Flip & reverse sumstats, and then merge w/ allele matching:
  HEADDIR = "/oak/stanford/groups/smontgom/amarder/HarmonizeGWAS"
  source(paste0(HEADDIR,"/scripts/expand_sumstats_via_flip_reverse.R"))
  idx = nchar(sumstats$non_effect_allele)==1 & nchar(sumstats$effect_allele)==1
  sumstat_expanded = expand_sumstats_via_flip_reverse(sumstats[idx,])
  
  # Merge with alleles:
  ind <- as.data.table(geno)[, .(chr, pos)][as.data.table(sumstats)[, .(chr, pos)], on = .(chr, pos),
                                            which = TRUE]
  geno.sub <- as.data.table(geno)[unique(ind), .(chr, pos, a0, a1)]
  tmp2 = merge(sumstat_expanded,geno.sub,by.x=c("chr","pos","non_effect_allele","effect_allele"),by.y=c("chr","pos","a0","a1"))
  rm(sumstat_expanded)
  
  # Remove duplicates (from expand_sumstats_via_flip_reverse() function):
  tmp2 = tmp2[!duplicated(tmp2$iden),]
  
  # Rename:
  colnames(tmp2)[colnames(tmp2)=="non_effect_allele"] = "ReferenceAllele"
  colnames(tmp2)[colnames(tmp2)=="effect_allele"] = "AltAllele"
  colnames(tmp2)[colnames(tmp2)=="A0"] = "non_effect_allele"
  colnames(tmp2)[colnames(tmp2)=="A1"] = "effect_allele"
  tmp2$A0 = tmp2$ReferenceAllele
  tmp2$A1 = tmp2$AltAllele
  
  tmp3 = subset(sumstats,!(iden %in% tmp2$iden))
  if (nrow(tmp3)==0) { # if the dataframe were only SNVs
    sumstats = tmp2
  } else {
    tmp3$ReferenceAllele = NA
    tmp3$AltAllele = NA
    tmp3$A0 = tmp3$non_effect_allele
    tmp3$A1 = tmp3$effect_allele
    
    # Bind:
    sumstats = rbind(tmp2,tmp3)
  }
  rm(tmp2);rm(tmp3)
  
  # SNP ids for identity:
  if (parallel) {
    num_cores <- detectCores()/2
    sumstats$snp_id = unlist(
      mclapply(
        1:nrow(sumstats),
        function(i) paste(
          # sumstats$chr[i],sumstats$pos[i],sumstats$A0[i],sumstats$A1[i],sep = "_"
          paste0("chr",sumstats$chr[i]),sumstats$pos[i],sumstats$A0[i],sumstats$A1[i],sep = ":"
        ),
        mc.cores = num_cores
      )
    )
  } else {
    sumstats$snp_id = unlist(
      lapply(
        1:nrow(sumstats),
        function(i) paste(
          # sumstats$chr[i],sumstats$pos[i],sumstats$A0[i],sumstats$A1[i],sep = "_"
          paste0("chr",sumstats$chr[i]),sumstats$pos[i],sumstats$A0[i],sumstats$A1[i],sep = ":"
        )
      )
    )
  }
  
  # order:
  sumstats = sumstats[order(sumstats$iden),]
  colnames(sumstats)[colnames(sumstats)=="pos"] <- "snp_pos"
  
  
  return(sumstats)
}

