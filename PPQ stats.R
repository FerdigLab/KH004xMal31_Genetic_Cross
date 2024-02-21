### Linear Model on PPQ IC50 data
# load packages
library(dplyr)
library(ggplot2)


# Set working directory  ---> Update this
setwd("/Users/kbuttons/Documents/Research/P01/Data/Project 2/Mal31XKH004/2022/github")

# read in PPQ data
PPQ_rep_data <- read.delim("IC50data_PPQPheno.csv",header=TRUE,sep=",",as.is=TRUE)

# trim NA points (PSA was assayed on different days for some progeny, so subset seperately)
PPQ_rep_data_trimmed <- PPQ_rep_data[which(is.na(PPQ_rep_data[,"LIMITED.POINT.IC50"])==FALSE),]
PPQ_rep_data_trimmed_PSA <- PPQ_rep_data[which(is.na(PPQ_rep_data[,"PSA"])==FALSE),]



# Summarize
MeanByParasite <- PPQ_rep_data %>% group_by(Parasite) %>% summarize(mean.IC50 = mean(PPQ.IC50, na.rm=TRUE), sd.IC50 = sd(PPQ.IC50, na.rm=TRUE), 
                                                                    mean.AUC = mean(AUC, na.rm=TRUE), sd.AUC = sd(AUC, na.rm=TRUE), 
                                                                    mean.PSA = mean(PSA, na.rm=TRUE), sd.PSA = sd(PSA, na.rm=TRUE), n())

MeanByParasiteProg <- PPQ_rep_data %>% group_by(Parasite,Geno.Code.CRT.PMCN,CRT.geno,Plasmepsin.CNV,Plasmepsin.geno) %>% summarize(LIMITED.POINT.IC50 = mean(LIMITED.POINT.IC50, na.rm=TRUE), sd.LIMITED.POINT.IC50 = sd(LIMITED.POINT.IC50, na.rm=TRUE), 
                                                                                                                                   AUC = mean(AUC, na.rm=TRUE), sd.AUC = sd(AUC, na.rm=TRUE), 
                                                                                                                                   PSA = mean(PSA, na.rm=TRUE), sd.PSA = sd(PSA, na.rm=TRUE), 
                                                                                                                                   IC50 = mean(PPQ.IC50, na.rm=TRUE), sd.IC50 = sd(PPQ.IC50, na.rm=TRUE), n())

write.csv(MeanByParasiteProg,"MeanPPQByProgeny.csv")



# split samples by CRT allele
CRTK <- which(MeanByParasiteProg[,"CRT.geno"]=="K")
CRTM <- which(MeanByParasiteProg[,"CRT.geno"]=="M")

# AUC stats full model
AUC_full <- anova(lm("AUC~CRT.geno*Plasmepsin.CNV",data=MeanByParasiteProg)) #significant interaction and non-significant main effect

# AUC stats for main effects
AUC_crt <- anova(lm("AUC~CRT.geno",data=MeanByParasiteProg))
AUC_PMCN <- anova(lm("AUC~Plasmepsin.geno",data=MeanByParasiteProg))

# AUC stats PMCN by crt allele
AUC_KCRT <- anova(lm("AUC~Plasmepsin.CNV",data=MeanByParasiteProg[CRTK,]))
AUC_MCRT <- anova(lm("AUC~Plasmepsin.CNV",data=MeanByParasiteProg[CRTM,]))

AUC_cor <- c(AUC_KCRT$`Pr(>F)`[1],AUC_MCRT$`Pr(>F)`[1],AUC_full$`Pr(>F)`[1])


# Limited Point IC50 full model
LMIC50_full <- anova(lm("LIMITED.POINT.IC50~CRT.geno*Plasmepsin.CNV",data=MeanByParasiteProg))

# LPIC50 stats for main effects
LMIC50_crt <- anova(lm("LIMITED.POINT.IC50~CRT.geno",data=MeanByParasiteProg))
LMIC50_PMCN <- anova(lm("LIMITED.POINT.IC50~Plasmepsin.geno",data=MeanByParasiteProg))

# LPIC50 stats PMCN by crt allele
LMIC50_KCRT <- anova(lm("LIMITED.POINT.IC50~Plasmepsin.CNV",data=MeanByParasiteProg[CRTK,]))
LMIC50_MCRT <- anova(lm("LIMITED.POINT.IC50~Plasmepsin.CNV",data=MeanByParasiteProg[CRTM,]))

LMIC50_cor <- c(LMIC50_KCRT$`Pr(>F)`[1],LMIC50_MCRT$`Pr(>F)`[1],LMIC50_full$`Pr(>F)`[1])


# PSA stats full model
PSA_full <- anova(lm("PSA~CRT.geno*Plasmepsin.CNV",data=MeanByParasiteProg))

# AUC stats for main effects
PSA_CRT <- anova(lm("PSA~CRT.geno",data=MeanByParasiteProg))
PSA_PMCN <- anova(lm("PSA~Plasmepsin.geno",data=MeanByParasiteProg))

# AUC stats PMCN by crt allele
PSA_KCRT <- anova(lm("PSA~Plasmepsin.CNV",data=MeanByParasiteProg[CRTK,]))
PSA_MCRT <- anova(lm("PSA~Plasmepsin.CNV",data=MeanByParasiteProg[CRTM,]))

PSA_cor <- c(PSA_KCRT$`Pr(>F)`[1],PSA_MCRT$`Pr(>F)`[1],PSA_full$`Pr(>F)`[1])

