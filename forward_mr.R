library(data.table)
library(dplyr)
library(glue)
library(tidyr)
library(gwasvcf)
library(genetics.binaRies)
library(stringr)
library(readr)
library(parallel)

#upload UKBB data from Synapse https://www.synapse.org/Synapse:syn51364943/wiki/622119

files <- list.files(path = "/user/work/ry20205/ukbb_data", pattern = "\\.tar$", full.names = T)
source("synapse_install.py")

df <- data.frame(matrix(ncol = 17, nrow = 0))
colnames(df) <- c("GENPOS",  "CHROM.x ID",  "ALLELE0", "ALLELE1", "A1FREQ",  "INFO",    "N",   "TEST",    "BETA",    "SE",  "CHISQ",   "LOG10P",  "EXTRA",   "pval",    "CHROM.y", "SNP", "gene")

gene<- df$gene[1]
split <- str_split(gene, "/")[[1]][[6]]
split_end <- str_split(split, "_")[[1]][[1]]
df$gene <- split_end
exp <- format_data(
    df,
    type = "exposure",
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
    chr_col = "CHROM.x",
    pos_col = "GENPOS", phenotype="gene"
)
exp_clump <- subset(exp, pval.exposure < 5e-8)
if (dim(exp_clump)[1]>0){
mart <- useEnsembl("ensembl","hsapiens_gene_ensembl", host = "https://www.ensembl.org")
table <- try(getBM(attributes = c("hgnc_id", "hgnc_symbol","start_position","end_position", "chromosome_name"), 
               values = unique(exp_clump$SNP), 
               mart = mart)) #cis SNPs for hg38
   table <- table[c("hgnc_id", "hgnc_symbol", "start_position", "end_position", "chromosome_name")]
   table <- subset(table, hgnc_symbol == exp_clump$exposure)
    table <- table[1,]
   merged <- merge(exp_clump, table, by.x = "exposure", by.y = "hgnc_symbol", all.x = T, all.y = F, no.dups = T)
   merged$pos.exposure<- as.numeric(merged$pos.exposure)
   merged$start_position <- as.numeric(merged$start_position)
   merged$end_position <- as.numeric(merged$end_position)
   merged$cistrans <- "trans"
   merged$cistrans <- with(merged,
                        ifelse(pos.exposure >= start_position - 500000 & pos.exposure <= end_position + 500000,
                               "cis", "trans"))
  subset <- subset(merged, cistrans == "cis")
if (dim(subset)[1]>0){
## do this for specific gene name
ids <- unique(subset$exposure)
clump<- ld_clump(
    dplyr::tibble(rsid=subset$SNP, pval=subset$pval.exposure, id=subset$exposure),
    plink_bin = genetics.binaRies::get_plink_binary(),
    bfile = "/user/work/ry20205/new_hg38/EUR"
)
clump$SNP <- clump$rsid
exp_dat <- merge(clump, subset, by="SNP")
if (dim(exp_dat)[1]>0){
out <- read.table("~/decode_analysis/data/all_format.txt", header=TRUE)
out$beta.outcome<-as.numeric(out$beta.outcome)
out$se.outcome<-as.numeric(out$se.outcome)
out$eaf.outcome<-as.numeric(out$eaf.outcome)
out$pval.outcome<-as.numeric(out$pval.outcome)
out$samplesize.outcome <- "24495" # 12,496(gbm=6191, nongbm=6305) and 18190 controls sample size of SNP should all "30686" gbm "24381" nongbm "24495"be total case controls
dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = out)
if (dim(dat)[1]>1){
steiger <- steiger_filtering(dat)
steiger_true <- subset(steiger, steiger_dir == TRUE) 
res <- generate_odds_ratios(mr(steiger_true))
res <- generate_odds_ratios(res$result)
out <- paste0("/user/work/ry20205/ukbb_data/results/", split_end, "_all_sig_mr_cis.txt")
write.table(res, file=out, row.names = F, quote = F, sep = "\t")} else {
res <- generate_odds_ratios(mr(dat))
out <- paste0("/user/work/ry20205/ukbb_data/results/", split_end, "_all_sig_mr_cis.txt")
write.table(res, file=out, row.names = F, quote = F, sep = "\t")}}}}
