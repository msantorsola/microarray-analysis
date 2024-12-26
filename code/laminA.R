#!/usr/bin/env R

library(sscore)
options(width=60)

library(affy)


setwd('/gpfs/scratch/userexternal/msantors/GSE124409_LMNB1-KD_MDA-MB-231/laminA')


# import CEL files
celpath = "/gpfs/scratch/userexternal/msantors/GSE124409_LMNB1-KD_MDA-MB-231/laminA" 
data = ReadAffy(celfile.path=celpath)



#Retrieving pheno data using affy
ph = data@phenoData
ph


#sample names 

ph@data[ ,1] = c("preLA","ctrl","ctrl")
ph


#Normalization using MAS 5.0
data.mas5 <- mas5(data) 
class(data.mas5)
#[1] "ExpressionSet"
#attr(,"package")
#[1] "Biobase"

#expression matrix (probesets/genes in rows, chips in columns) 
exprSet.nologs = exprs(data.mas5)

# chip names
colnames(exprSet.nologs)

#log-transformation of the expression values 
exprSet = log(exprSet.nologs, 2)

#save expression matrix 
write.table(exprSet, file="prelamA_mas5_matrix_log.txt", quote=F, sep="\t")

#Absent/Present call by Affy A/P call algorithm on the CEL files
data.mas5calls = mas5calls(data)

#get the A/P calls
data.mas5calls.calls = exprs(data.mas5calls)

# save A/P calls
write.table(data.mas5calls.calls, file="prelamA_A_P_mas5calls.txt", quote=F, sep="\t")


summary(data.mas5calls.calls)


#Differential Expression Analysis

#Multichip comparisons

labels <- c(0,1,1) #preLAMA vs 2 ctrl
SScore.multi <- SScore(data,classlabel=c(0,1,1))

#S-Score values
sscores <- exprs(SScore.multi)
head(sscores)
#Class 0 vs 1
write.table(sscores, file="prelamA_DEG_all.txt", quote=F, sep="\t")
          
#Calculate p-values corresponding to rejection of the null hypothesisand acceptance of the alternative hypothesis of differential gene expression. 
#Cutoff valuesfor the S-Scores can be set to achieve the desired level of significance.  
#As an example,an absolute S-Score value of 3 (signifying 3 standard deviations from the mean, a typical cutoff value) would correspond to a p-value of 0.003.  
#Under this scenario, the significantgenes can be found as:

signif <- geneNames(data)[abs(sscores) >= 3]
write.table(signif, file="prelamA_DEG_sscore_3_sign.txt", quote=F, sep="\t")

