library(ieugwasr) # For clumping
library(TwoSampleMR) # The MR-Base R package
library(gwasvcf)
library(parallel)
library(biomaRt)
library(foreach)
library(dplyr)
library("readxl")
biomartCacheClear()

#load files
out <- read.table("~/decode_analysis/data/all_format.txt", header=TRUE)

protein_files <- list.files(path = "vcf",
                            pattern = "\\.vcf.bgz$",
                            full.names = T)
mart <- useEnsembl("ensembl","hsapiens_gene_ensembl", host = "https://www.ensembl.org")

for (i in protein_files) {
  gene <- strsplit(i, "_")[[1]][3]
  start <- lapply(i, function(x) as.character(strsplit(x, "/", fixed = T)[[1]][13]))
  check <- lapply(start, function(x) as.character(strsplit(x, ".", fixed = T)[[1]][1]))
  print(i)
  vcf <- readVcf(i)
  vcfsubset <- query_gwas(vcf, rsid = out$SNP) #add all glioma snps
                  
  if (dim(vcfsubset)[1] > 0){
  exp_dat <- vcfsubset %>% gwasglue::gwasvcf_to_TwoSampleMR()
  exp_dat['phenotype'] <- gene
  exp_clump <- subset(exp_dat, pval.exposure < 5e-8)
    
  if (dim(exp_clump)[1] > 0){
   table <- try(getBM(attributes = c("hgnc_id", "hgnc_symbol","start_position","end_position", "chromosome_name"), values = unique(exp_clump$SNP), mart = mart)) #cis SNPs for hg38
   table <- table[c("hgnc_id", "hgnc_symbol", "start_position", "end_position", "chromosome_name")]
   table <- subset(table, hgnc_symbol == exp_clump$phenotype)
    table <- table[1,]
   merged <- merge(exp_clump, table, by.x = "phenotype", by.y = "hgnc_symbol", all.x = T, all.y = F, no.dups = T)
   merged$pos.exposure<- as.numeric(merged$pos.exposure)
   merged$start_position <- as.numeric(merged$start_position)
   merged$end_position <- as.numeric(merged$end_position)
   merged$cistrans <- "trans"
   merged$cistrans <- with(merged,
                        ifelse(pos.exposure >= start_position - 500000 & pos.exposure <= end_position + 500000,
                               "cis", "trans"))
  subset <- subset(merged, cistrans == "cis") #find cis SNPs
    
if (dim(subset)[1] > 0){
subset$chromosome_name <- paste0("chr", subset$chromosome_name)
check <- subset(subset, chromosome_name == chr.exposure)

if(dim(check)[1]>0){
  names(check)[names(check) == "pval.exposure"] <- "pval"
  check$chr.exposure <- gsub("chr*", "\\1", check$chr.exposure)
  check$chrpos<-paste0(check$chr.exposure,":",check$pos.exposure)
  check$rsid <- check$SNP 
  check$id<-check$id.exposure
  exp <- check
  print(i)
  
if (dim(exp)[1]>0){
ids <- unique(exp$phenotype)

clump<- ld_clump( #clump exposure data
    dplyr::tibble(rsid=exp$rsid, pval=exp$pval, id=exp$phenotype),
    plink_bin = genetics.binaRies::get_plink_binary(),
    bfile = "EUR"
)

print(clump)
exp_dat <- merge(clump, exp, by="rsid")
  
dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = out) #harmonise data
  
res <- mr(dat) #carry out MR
  
start <- lapply(i, function(x) as.character(strsplit(x, "/", fixed = T)[[1]][13]))
check <- lapply(start, function(x) as.character(strsplit(x, ".", fixed = T)[[1]][1])) 
out <- paste0(check, "_all_sig_mr_cis.txt")
write.table(res, file=out)
}}}}}}
