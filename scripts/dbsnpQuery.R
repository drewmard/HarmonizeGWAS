# zcat snp151.txt.gz | cut -f2,4,5,10 > snp151.v2.txt
# gzip -c snp151.v2.txt > snp151.v2.txt.gz

# zcat snp151.txt.gz | head -10 | cut -f2,3,4,5,10 | bgzip -c > snp151.v2.txt.bgz
# tabix -p bed snp151.v2.txt.bgz
# 
# drop 2 from snp id for start, this will capture too much but im having trouble w/ deletions


# data_input=df2
# type="rsid"
# tmpdir="~/tmp"
# SNPFILE="/oak/stanford/groups/smontgom/amarder/data/dbsnp/snp151.v2.txt.gz"


####################################################

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
  
  #ss4 becomes original + reversed + strand flip original + strand flip reverse sumstats
  ss4 <- rbind(ss3, ss4) ######
  return(ss4)
}


dbsnpQuery = function(data_input,
                      type="rsid",
                      tmpdir="~/tmp",
                      trait="tmp",
                      SNPFILE="/oak/stanford/groups/smontgom/amarder/data/dbsnp/snp151.v2.txt.gz",
                      ignoreAlleles=FALSE) {
  if (type=="rsid") {
    
    print("Querying dbSNP using chr + pos...")
    
    y = data_input$RSID
    RSIDFILE=paste0(tmpdir,"/rsid.txt")
    RSIDTMPFILE=paste0(tmpdir,"/rsidtmp.txt")
    MAPFILE=paste0(tmpdir,"/rsid_mapping.txt")
    fwrite(as.data.frame(data_input$RSID),RSIDFILE,quote = F,na = "NA",sep = '\t',row.names = F,col.names = F)
    
    cmd = paste0("awk '{print $0 \"\\t\"}' ",RSIDFILE, " > ",RSIDTMPFILE)
    system(cmd)
    
    # note: this is not an exact match, as rs1010989 will capture rs10109893 and rs101098934 and rs1010989343
    # this is fine because we are going to merge!
    # system(paste0("zgrep -f ",RSIDFILE," ",SNPFILE," > ",MAPFILE))
    cmd = paste0("zcat ",SNPFILE," | grep -Ff ",RSIDTMPFILE," > ",MAPFILE)
    system(cmd)
    # cmd = paste0("zcat ",SNPFILE," | grep -Ff <(awk '{print $0 \"\\t\"}' ",RSIDTMPFILE,") > ",MAPFILE);
    
    # Read in RSID info and QC:
    RSIDINFO = fread(MAPFILE,data.table = F,stringsAsFactors = F)
    colnames(RSIDINFO) = c("CHR","POS","RSID","A1/A2")
    # colnames(RSIDINFO) = c("CHR","POS","RSID")
    RSIDINFO$CHR = sub("chr","",RSIDINFO$CHR)
    RSIDINFO = stringSplitter(RSIDINFO,
                              column_to_use = "A1/A2",
                              divider="/",
                              removeCol = TRUE,idCol = FALSE,
                              idx_to_use = c(1,2),cols_to_use = c("A1","A2"))
    RSIDINFO = subset(RSIDINFO,CHR %in% c(1:22,"X","Y"))
    RSIDINFO = subset(RSIDINFO,!duplicated(RSID))
    if (ignoreAlleles) {
      RSIDINFO = unique(RSIDINFO[,c("CHR","POS","RSID")])
    }
    
    # Merge data input and dbsnp
    if (!ignoreAlleles) {
      data_input[,"tmp"] = 1:nrow(data_input)
      data_input1a = merge(data_input[,!(colnames(data_input)%in%c("CHR","POS"))],RSIDINFO,by=c("RSID","A1","A2"))
      data_input1b = merge(data_input[,!(colnames(data_input)%in%c("CHR","POS"))],RSIDINFO,by.x=c("RSID","A1","A2"),by.y=c("RSID","A2","A1"))
      data_input1 = as.data.frame(rbind(data_input1a,data_input1b))
      data_input1c = merge(subset(data_input,!(RSID %in% data_input1$RSID))[,!(colnames(data_input)%in%c("CHR","POS"))],RSIDINFO[,!(colnames(RSIDINFO)%in%c("A1","A2"))],by="RSID")
      data_input1 = as.data.frame(rbind(data_input1,data_input1c))
      data_input1d = subset(data_input,!(tmp %in% data_input1$tmp))     # data_input1d = subset(data_input,!(RSID %in% data_input1$RSID))
      data_input1 = as.data.frame(rbind(data_input1,data_input1d))
      data_input = data_input[,!(colnames(data_input) %in% c("tmp"))]; data_input1 = data_input1[,!(colnames(data_input1) %in% c("tmp"))]
    } else {
      data_input1 = merge(data_input[,!(colnames(data_input)%in%c("CHR","POS"))],RSIDINFO,by="RSID",all.x=TRUE)
    }
    
    # Final formatting
    data_input1 = data_input1[,colnames(data_input)]
    data_input1$hg = "hg38"
    data_input1$CHRPOS = paste0(data_input1$CHR,"_",data_input1$POS)
    
    # Remove temp files
    system(paste0("rm ",RSIDFILE))
    system(paste0("rm ",RSIDTMPFILE))
    system(paste0("rm ",MAPFILE))
  }
  
  if (type=="chr_pos") {
    
    print("Querying dbSNP using chr + pos...")
    
    f = SNPFILE
    dbsnp = fread(SNPFILE,data.table = F,stringsAsFactors = F)
    
    print("---Preprocess dbsnp file...")
    y = dbsnp$V5
    y[grep("/",y,invert = T)] = "?/?"
    y = strsplit(y,"/")
    chrNum = sub("chr","",dbsnp$V1)
    dbsnp$id = paste(
      chrNum, # remove chr!
      dbsnp$V3,
      unlist(lapply(y,function(x)x[[1]])),
      unlist(lapply(y,function(x)x[[2]])),
      sep = '_'
    )
    dbsnp = dbsnp[!duplicated(dbsnp$id),]
    colnames(dbsnp) = c("chr","start","snp_pos","rsid","alleles","id")
    
    # merge based on chr, pos, a1, a2
    print("---Merge based on chr, pos, a1, a2...")
    data_input$orig_id = paste0(data_input$chr,"_",data_input$snp_pos,"_",data_input$non_effect_allele,"_",data_input$effect_allele)
    sumstat_expanded = expand_sumstats_via_flip_reverse(data_input)
    sumstat_expanded$new_id = paste0(sumstat_expanded$chr,"_",sumstat_expanded$snp_pos,"_",sumstat_expanded$non_effect_allele,"_",sumstat_expanded$effect_allele)
    data_out <- merge(sumstat_expanded, dbsnp[,c("id","rsid")], 
                      by.x=c("new_id"),by.y=c("id"), all = FALSE)
    data_out = data_out[,colnames(data_out)!="new_id"]
    
    # merge based on chr, pos for other snps
    # need to then add sumstat_expanded, which didnt merge with dbsnp:
    # tmp = those that didnt get merged
    print("---Merge based on chr + pos...")
    tmp = subset(data_input,!(orig_id %in% data_out$orig_id))
    if (length(tmp) > 0) { 
      data_out2 = merge(tmp,dbsnp[,c("rsid","chr","snp_pos")],
                        by=c("chr","snp_pos"))
      data_out2 = data_out2[!duplicated(data_out2$rsid),]
      data_out = as.data.frame(rbind(data_out,data_out2))
    }
    
    tmp = subset(data_input,!(orig_id %in% data_out$orig_id))
    if (length(tmp)>0) {
      tmp$rsid = tmp$orig_id
      data_out = rbind(data_out,tmp)
    }
    # sort by chr and position, and remove orig_id column
    print("---Sorting...")
    data_out = data_out[order(data_out$chr,data_out$snp_pos),]
    
    # Reformat
    data_out = data_out[!duplicated(data_out$orig_id),]
    data_out = data_out[,colnames(data_out)!="orig_id"]
    data_out = data_out[,c("chr","snp_pos","non_effect_allele","effect_allele","rsid",
                           colnames(data_out)[!(colnames(data_out) %in% c("chr","snp_pos","non_effect_allele","effect_allele","rsid"))])]

    # return:
    data_input1 = data_out
    
  }
  print("dbsnpQuery complete.")
  return(data_input1)
}
