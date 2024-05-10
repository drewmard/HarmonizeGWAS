expand_sumstats_via_flip_reverse <- function(ss) {
  # Next, account for strand flip + reverse alleles:
  # When comparing sumstats (GWAS) to genotype ref (1000G),
  # A/C in sumstats and C/A in the ref is fine,
  # the sumstats just need to be reversed (multiplying Z-score by -1).
  # However, A/C and A/G is not fine.
  # A/G and T/C is the same, just needs to be flipped according to strand.
  # Same with A/G and C/T, but need to flip AND reverse.
  # This function handles all the above.
  # However:
  # We are going to keep SNPs on ambiguous strands (e.g. A/T, C/G),
  # but maybe these should be removed?
  
  flip_strand <- function(allele) {
    dplyr::case_when(
      allele == "A" ~ "T",
      allele == "C" ~ "G",
      allele == "T" ~ "A",
      allele == "G" ~ "C",
      TRUE ~ NA_character_
    )
  }
  
  # ss is original sumstats
  
  # ss3 is reversed sumstats
  ss3 <- ss
  ss3$effect_allele <- ss$non_effect_allele
  ss3$non_effect_allele <- ss$effect_allele
  if ("z" %in% colnames(ss3)) {
    ss3$z <- -ss$z
  }
  if ("beta" %in% colnames(ss3)) {
    ss3$beta <- -ss$beta
  } 
  if ("or" %in% colnames(ss3)) {
    ss3$or <- exp(-log(ss$or))
  } 
  
  
  #ss3 becomes original + reversed sumstats
  ss3 <- rbind(ss, ss3) #####
  
  # for non-ambiguous snps, do strand flip:
  ss4 <- ss3[!((ss3$non_effect_allele == "A" & ss3$effect_allele == "T") |
                 (ss3$non_effect_allele == "T" & ss3$effect_allele == "A") |
                 (ss3$non_effect_allele == "G" & ss3$effect_allele == "C") |
                 (ss3$non_effect_allele == "C" & ss3$effect_allele == "G")),]
  ss4$effect_allele <- flip_strand(ss4$effect_allele)
  ss4$non_effect_allele <- flip_strand(ss4$non_effect_allele)
  
  # for 
  
  #ss4 becomes original + reversed + strand flip original + strand flip reverse sumstats
  ss4 <- rbind(ss3, ss4) ######
  return(ss4)
}
