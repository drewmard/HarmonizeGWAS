# Assuming your data is in a data frame named 'data_input'
# Replace this with your actual data import logic

# Create a basic VCF header
vcf_header <- c(
  '##fileformat=VCFv4.2',
  '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',
  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
  '##contig=<ID=1>',
  '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE'
)

# Create a data frame with placeholder values
data_input.tmp = data_input[1:10,]
vcf_data <- data.frame(
  CHROM = 1,
  POS = data_input.tmp$snp_pos,
  ID = '.',
  REF = data_input.tmp$non_effect_allele,
  ALT = data_input.tmp$effect_allele,
  QUAL = '.',
  FILTER = '.',
  INFO = 'AF=0.5',
  FORMAT = 'GT',
  SAMPLE = '0/1'
)

# Combine the header and data
vcf_content <- c(vcf_header, apply(vcf_data, 1, function(row) paste(row, collapse = '\t')))

# Write the content to a VCF file
writeLines(vcf_content, '~/tmp/tmp.vcf')

inputVCF=/home/amarder/tmp/tmp.vcf
outputVCF=/home/amarder/tmp/tmp.out.vcf
annot1=9
echtvar anno -e $annot1 $inputVCF $outputVCF




