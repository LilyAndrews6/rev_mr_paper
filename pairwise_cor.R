##create pairwise plot to compare protein betas betwen studies
pdf(file="pairwise_plot.pdf) #save as a pdf file
pairs(protein_betas[1:3], panel=function(x,y){ #protein_beta 1:3 are the columns b_ukbb, b_decode and b_interval
  points(x,y)
  abline(lm(y~x), col="red")})
dev.off()
