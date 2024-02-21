install.packages("qtl")
setwd('/Users/John/Desktop/NotreDame/Ferdig Lab/QTL scans/KHxMal_cross/')      
library(dplyr)
library(qtl)
##your input##
PHENOTYPES = c('Fitness_Elo')

cross_data <- read.cross(format="csv",file="NewMap for crt allele split scans ppq paper MapProgMapSNPscleaned_09-10-2023.csv",na.strings="NA",genotypes=c(0,1))

PHENOTYPES = c(c(colnames(cross_data[["pheno"]]))[2:22])#length(cross_data[["pheno"]])])#[4:8])
print(PHENOTYPES)
setwd('/Users/John/Desktop/')      
for (i in 1:length(PHENOTYPES)){
  phenotype_column=PHENOTYPES[i] 
  
###########################################################################
  out <- scanone(cross_data, method="hk", pheno.col=phenotype_column)
  operm <- scanone(cross_data, method="hk", pheno.col=phenotype_column, n.perm=1000, verbose=FALSE)
  
  cutofflim <- summary(operm, alpha=c(0.37, 0.05, 0.01))
  
  summary(out)
  summary(out, threshold=2)
  ##
  total_num_of_markers=totmar(cross_data)
  ##
  cutoff1 <- out
  cutoff1$lod <- rep(cutofflim[1],total_num_of_markers)#changed from 1109 for kv file to 1044 for kbs file
  
  cutoff2 <- out
  cutoff2$lod <- rep(cutofflim[2],total_num_of_markers)
  
  cutoff3 <- out
  cutoff3$lod <- rep(cutofflim[3],total_num_of_markers)
  
  png(file=paste(phenotype_column,".png",sep=""), height=8, width=11, units="in", res=300)
  plot(out, bandcol="gray90", col="black", bgrect="white", main=phenotype_column, ylab="LOD score", xlab="Chromosome number (cM)", ylim=c(0,8), cex.lab = 1.25, cex.axis = 1.4)
  plot(cutoff1, lty=1, add=TRUE, col="black")
  plot(cutoff2, lty=1, add=TRUE, col="black")
  plot(cutoff3, lty=1, add=TRUE, col="black")
  text(x=1000608,y=(cutoff1$lod[1]-0.05),"63%", pos=3, cex=1.15)
  text(x=1000608,y=(cutoff2$lod[1]-0.05),"95%", pos=3, cex=1.15)
  text(x=1000608,y=(cutoff3$lod[1]-0.05),"99%", pos=3, cex=1.15)
  dev.off()
  
  LOD=lodint(out,as.numeric(as.character(unlist(max(out)[1]))),expandtomarkers=TRUE)###highest LOD 
  write.csv(LOD, paste(phenotype_column, "1.5LODforGvis.csv", sep=""))

  write.csv(out, paste("LOD_Scores_", phenotype_column, ".csv", sep=""))
}


######################################################################################
#Calculate and plot residuals
#run the scanone for one phenotype before running this
#this script uses the "out" variable from the scanone, so you have to run the main scan
#for a single phenotype and then run this. 

#find marker with the highest LOD score
max(out1.c4)

#set this marker to be your covariate and then rerun the scanone function
g<- pull.geno(cross_data) [,"M443"]
out1.c4 <- scanone(cross_data, method="hk", pheno.col=phenotype_column, addcovar = g)


summary(out1.c4)
summary(out1.c4, threshold=2)

total_num_of_markers=totmar(cross_data)

cutoff1 <- out1.c4
cutoff1$lod <- rep(cutofflim[1],total_num_of_markers)#changed from 1109 for kv file to 1044 for kbs file

cutoff2 <- out1.c4
cutoff2$lod <- rep(cutofflim[2],total_num_of_markers)

cutoff3 <- out1.c4
cutoff3$lod <- rep(cutofflim[3],total_num_of_markers)

png(file=paste(phenotype_column,"_residual.png",sep=""), height=8, width=11, units="in", res=300)
plot(out1.c4, bandcol="gray90", col="black", bgrect="white", main=phenotype_column, ylab="LOD score", xlab="Chromosome number (cM)", ylim=c(0,8), cex.lab = 1.25, cex.axis = 1.4)
plot(cutoff1, lty=1, add=TRUE, col="black")
plot(cutoff2, lty=1, add=TRUE, col="black")
plot(cutoff3, lty=1, add=TRUE, col="black")
text(x=1000608,y=(cutoff1$lod[1]-0.05),"63%", pos=3, cex=1.15)
text(x=1000608,y=(cutoff2$lod[1]-0.05),"95%", pos=3, cex=1.15)
text(x=1000608,y=(cutoff3$lod[1]-0.05),"99%", pos=3, cex=1.15)
dev.off()

LOD=lodint(out1.c4,as.numeric(as.character(unlist(max(out)[1]))),expandtomarkers=TRUE)###highest LOD 
write.csv(LOD, paste(phenotype_column, "_residual_1.5LODforGvis.csv", sep=""))

write.csv(out, paste("LOD_Scores_residual_", phenotype_column, ".csv", sep=""))

####If you want to make a plot that directly compares the main scan to residual scan
plot(out, out1.c4, col=c("blue", "red"))

### New 2D scan
out2 <- scantwo(cross_data, method="hk", pheno.col="AUC_Raw")
operm2 <- scantwo(cross_data, method="hk", pheno.col="AUC_Raw", n.perm=5) #change 5 here to do a
summary(out2, perms=operm2, alpha=0.2, pvalues=TRUE)
dev.new()
png("2D RSA QTL vr1.png", height=8, width=11, units="in", res=220)
plot(out2,
     chr=c(5,6,7,9,10,11,14),
     main="PQP.AUC Pairwise Scan", col.scheme = "viridis")
dev.off()

plot(out2)
plot(operm2)