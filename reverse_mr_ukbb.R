library(data.table)
library(dplyr)
library(glue)
library(tidyr)
library(gwasvcf)
library(genetics.binaRies)
library(stringr)
library(readr)
library(parallel)
library(biomaRt)

#upload UKBB data from Synapse https://www.synapse.org/Synapse:syn51364943/wiki/622119

source("synapse_install.py")
files <- list.files(path = "ukbb_pqtl", pattern = "\\.tar$", full.names = T) #list of tar files from UKBB protein data

all <- read.table("all_format.txt", header=T) #all glioma data formatted for forward MR

untar <- glue("tar -xvf {i}") #untar protein files
system(untar)
split <- str_split(i, "/")[[1]][[6]]
split_final <- str_sub(split, end = -5)
combine <- paste0("ukbb_pqtl", split_final)
zip_files <- list.files(path = combine, pattern = "\\.gz$", full.names = T)

for (x in zip_files){
unzip <- glue("gzip -d {x}") #unzip protein files
system(unzip)
print("unzipped")
   chr_file <- str_sub(x, end = -4)
prot <- read.table(chr_file, header=T)
print("file read")
    prot$pval <- 10^(-prot$LOG10P) #calculate p-value
    marker <- prot$ID
    position <- data.frame(do.call("rbind", strsplit(as.character(marker), ":", fixed = TRUE)))
    subset <- position[, 1:4]
    prot$CHROM <- subset$X1 #add chr 
    prot$GENPOS <- subset$X2 #add pos 
    prot$A0 <- subset$X3
    prot$A1 <- subset$X4
    prot$pos.outcome <- prot$GENPOS
    df <- merge(prot, all, by="pos.outcome")
    if (dim(df)[1]>0){
    new_df <- df %>% filter(CHROM == chr.outcome)
    if (dim(new_df)[1]>0){
    all(new_df$CHROM == new_df$chr.outcome) #include SNPs which map to all glioma file
    file <- str_split(chr_file, "/")[[1]][[7]]
    out <- paste0("mapped_", file, ".txt")
    new_df$gene <- i
    write.table(new_df, file=out, row.names=F, quote = F, sep = "\t") #save mapped file
}}}

txt_files <- list.files(path = "ukbb_pqtl", pattern = "\\.txt$",full.names = T) #text files which have been mapped

for (t in txt_files){
print(t)
  t_split <- str_split(t, "/")[[1]][[7]]
t_split_end <- str_split(t_split, "_")[[1]][[4]]
empty_list = rbind(empty_list, t_split_end)}
list <- unique(empty_list[,1])
print(list)

for (x in list) {
num <- grep(x, txt_files)
files <- txt_files[num]
df <- data.frame(matrix(ncol = 29, nrow = 0))
column_title <- c("pos.outcome", "CHROM",   "GENPOS",  "ID",  "ALLELE0", "ALLELE1", "A1FREQ",  "INFO",    "N",   "TEST",    "BETA",    "SE",  "CHISQ",   "LOG10P", "EXTRA",   "pval",    "SNP", "chr.outcome", "effect_allele.outcome",   "other_allele.outcome",    "eaf.outcome", "beta.outcome",    "se.outcome",  "pval.outcome",    "outcome", "mr_keep.outcome", "pval_origin.outcome", "id.outcome",  "gene")
colnames(df) <- column_title
  
for (i in files){
df_temp = read.table(i,fill=T, header=T, sep = "\t")
df = rbind(df, df_temp)}

outcome_dat <- format_data( #read protein data as outcome
    df,
    type = "outcome",
    snps = NULL,
    header = TRUE,
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    eaf_col = "A1FREQ",
    effect_allele_col = "ALLELE1",
    other_allele_col = "ALLELE0",
    pval_col = "pval",
    samplesize_col = "N",
    chr_col = "CHROM",
    pos_col = "GENPOS", phenotype="gene"
)

exp <- read.table("~/decode_analysis/data/all_clumped_p1.txt", header=TRUE) #glioma data formatted for p1 threshold(r2=0.01, p-value=1x10-5)

dat <- harmonise_data(exposure_dat = exp, outcome_dat = outcome_dat)

#heterogeneity analysis
if (dim(dat)[1]>0){
het <- mr_heterogeneity(dat)
print(het)
loo <- mr_leaveoneout(dat)
out <- paste0(x, "_loo_all_p1.txt")
write.table(loo, file = out, row.names = F, quote = F, sep = "\t")

if (dim(dat)[1] > 0){
    res <- generate_odds_ratios(mr(dat))
    out <- paste0("results",x, "_all_p1_rev_mr.txt")
    write.table(res, file = out, row.names = F, quote = F, sep = "\t")
}}}

