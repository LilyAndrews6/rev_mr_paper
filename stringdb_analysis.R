#install packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("STRINGdb")
library(STRINGdb)

#pull human data from STRINGdb
string_db <- STRINGdb$new(version="12.0", species=9606)
string_db$ppi_enrichment(string_ids)

##for ukbb
ukbb_prot_list <- read.csv("ukbb_prot_list.txt", sep="")

ppi_list <- list()
end <- length(string_ids)
for (i in 1:100){
rand_prot <- sample(ukbb_prot_list$prot, end) #random selection of UKBB protein list
ppi <- string_db$ppi_enrichment(rand_prot) #find enrichment of randomly selected protein
ppi_list <- rbind(ppi_list, ppi$enrichment) #merge protein and enrichment data
}

total <- unlist(ppi_list)
mean <- mean(total)
print(mean)

##for deCODE
decode_prot_list <- read.csv("decode_prot_list.txt", sep="")

ppi_list <- list()
end <- length(string_ids)
for (i in 1:100){
rand_prot <- sample(decode_prot_list$prot, end) #random selection of deCODE protein list
ppi <- string_db$ppi_enrichment(rand_prot) #find enrichment of randomly selected protein
ppi_list <- rbind(ppi_list, ppi$enrichment) #merge protein and enrichment data
}
total <- unlist(ppi_list)
mean <- mean(total)
print(mean)

##for INTERVAL
interval_prot_list <- read.csv("interval_prot_list.txt", sep="")

ppi_list <- list()
end <- length(string_ids)
for (i in 1:100){
rand_prot <- sample(interval_prot_list$prot, end) #random selection of INTERVAL protein list
ppi <- string_db$ppi_enrichment(rand_prot) #find enrichment of randomly selected protein
ppi_list <- rbind(ppi_list, ppi$enrichment) #merge protein and enrichment data
}
total <- unlist(ppi_list)
mean <- mean(total)
print(mean)
