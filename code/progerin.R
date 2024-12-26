#!/usr/bin/env R


library(sscore)
options(width=60)

library(affy)

setwd('/gpfs/scratch/userexternal/msantors/GSE124409_LMNB1-KD_MDA-MB-231/progerin')
celpath = "/gpfs/scratch/userexternal/msantors/GSE124409_LMNB1-KD_MDA-MB-231/progerin"

# import CEL files data.p = ReadAffy(celfile.path=celpath)

#add annotation to phenoData
ph@data.p[ ,1] = c("ctrl","ctrl", "progerin")
ph

#Normalization using MAS 5.0
data.p.mas5 <- mas5(data.p) #normalized expression values for all probesets, along with other information. 
class(data.p.mas5)
#[1] "ExpressionSet"
#attr(,"package")
#[1] "Biobase"

exprSet.nologs = exprs(data.p.mas5)

colnames(exprSet.nologs)

#log-transform expression values 
exprSet = log(exprSet.nologs, 2)


#save_data
write.table(exprSet, file="progerin_mas5_matrix_log.txt", quote=F, sep="\t")

# Affy A/P call algorithm 
data.mas5calls = mas5calls(data.p)

# Get the  A/P calls
data.mas5calls.calls = exprs(data.mas5calls)

#save the calls as a matrix
write.table(data.mas5calls.calls, file="progerin_A_P_mas5calls.txt", quote=F, sep="\t")


labels <- c(1,1,0) #progerin vs 2 ctrl
SScore.multi <- SScore(data.p,classlabel=c(1,1,0))

#S-Score values
sscores <- exprs(SScore.multi)
head(sscores)
write.table(sscores, file="progerin_DEG_all.txt", quote=F, sep="\t", row.names=T)

#Calculate p-values 

signif <- geneNames(data.p)[abs(sscores) >= 3]
write.table(signif, file="progerin_DEG_sscore>3_sign.txt", quote=F, sep="\t")

