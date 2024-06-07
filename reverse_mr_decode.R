library(remotes)
library(gwasvcf)
library(GenomicFiles)
library(ieugwasr) # For clumping
library(TwoSampleMR) # The MR-Base R package
library(MRInstruments) 
library(tidyverse)
library(mrpipeline)
library(pegas)
library(R.utils)
library(vcfR)
library(parallel)

gwasvcf::set_bcftools("/sw/apps/bcftools-1.9/bin/bcftools")

#load files

all <- read.table("all_clumped_p1.txt", header=TRUE)

protein_files <- list.files(path = "vcf",
                            pattern = "\\.vcf.bgz$",
                            full.names = T)

for (i in protein_files){
  start <- lapply(i, function(x) as.character(strsplit(x, "/", fixed = T)[[1]][13]))
  check <- lapply(start, function(x) as.character(strsplit(x, ".", fixed = T)[[1]][1]))
  mr_files <- paste("result", check, "_all_sig_revmr.txt", sep='')
                  
  if (!file.exists(mr_files)){
  vcf <- readVcf(i)
  vcfsubset <- query_gwas(vcf, rsid = all$SNP) #add all glioma snps
    
  if (dim(vcfsubset)[1] > 0){
  out_dat <- vcfsubset %>% gwasglue::gwasvcf_to_TwoSampleMR() #extract data from vcf
  names(out_dat)[names(out_dat) == "chr.exposure"] <- "chr.outcome"
  names(out_dat)[names(out_dat) == "pos.exposure"] <- "pos.outcome"
  names(out_dat)[names(out_dat) == "other_allele.exposure"] <- "other_allele.outcome"
  names(out_dat)[names(out_dat) == "effect_allele.exposure"] <- "effect_allele.outcome"
  names(out_dat)[names(out_dat) == "beta.exposure"] <- "beta.outcome"
  names(out_dat)[names(out_dat) == "se.exposure"] <- "se.outcome"
  names(out_dat)[names(out_dat) == "pval.exposure"] <- "pval.outcome"
  names(out_dat)[names(out_dat) == "eaf.exposure"] <- "eaf.outcome"
  names(out_dat)[names(out_dat) == "samplesize.exposure"] <- "samplesize.outcome"
  names(out_dat)[names(out_dat) == "ncase.exposure"] <- "ncase.outcome"
  names(out_dat)[names(out_dat) == "ncontrol.exposure"] <- "ncontrol.outcome"
  names(out_dat)[names(out_dat) == "exposure"] <- "outcome"
  names(out_dat)[names(out_dat) == "mr_keep.exposure"] <- "mr_keep.outcome"
  names(out_dat)[names(out_dat) == "pval_origin.exposure"] <- "pval_origin.outcome"
  names(out_dat)[names(out_dat) == "id.exposure"] <- "id.outcome"
    
  dat <- harmonise_data(exposure_dat = all, outcome_dat = out_dat) #harmonise data
    
  res <- mr(dat) #perform MR
    
  outpath <- "result"
  start <- lapply(i, function(x) as.character(strsplit(x, "/", fixed = T)[[1]][13]))
  check <- lapply(start, function(x) as.character(strsplit(x, ".", fixed = T)[[1]][1])) 
  out <- paste0(check, "_all_sig_revmr.txt")
  write.table(res, file=paste0(outpath, out))
}}}
