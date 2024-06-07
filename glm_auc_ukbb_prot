#load packages
library(pROC)

#read in UKBB prot data and revMR results of associated proteins
ukbb_raw_control<- merge(data_participant_control, data_prot_control)
ukbb_raw_case<- merge(data_participant_case, data_prot_case)
ukbb_raw_case_control <- rbind(ukbb_raw_case, ukbb_raw_control)


##create weighted score per protein
df <- data.frame()
for (i in 1:nrow(sort.ukbb_revmr_res)){
  print(ukbb_revmr_res$outcome[i]) #check which protein is running
  score <- as.numeric(ukbb_raw_case_control[,i+3])*as.numeric(ukbb_revmr_res$b[i]) #create score multiplying protein level data by beta generated from revMR
  lm_model <- glm(ukbb_raw_case_control$treatment ~ score, family = "binomial") #run generalised linear model
  model <- coef(summary(lm_model))[2,] #output results
  df <- rbind(df,model)
  df$prot[i] <- ukbb_revmr_res$outcome[i]
  roc_object <- roc(response=ukbb_raw_case_control$treatment, predictor=score) #generate area under the curve estimate
  rev_roc <- pROC::auc(roc_object)
  df$roc[i] <- rev_roc
}
colnames(df)[1:4] <- c("Estimate", "Std.Error", "z-value", "Pr(>|z|)")
write.table(df, file="protein_score.txt", row.names = F, quote=F, sep="\t") #export results
