#!/usr/bin/env R

library(sscore)
options(width=60)


library(oligo)

GEO <- "GSE63059"
dir.geo <- paste(tempdir(),GEO,sep="/")
dir.create(dir.geo, showWarnings = FALSE)

list.cels <- c("GSM1539076","GSM1539077","GSM1539078","GSM1539079")

# pull the compressed .CEL files:
cels.get <- function(x)
		GEOquery::getGEOSuppFiles(GEO = x,
		makeDirectory = FALSE,
		baseDir = dir.geo,
		filter_regex = "*.CEL.gz")
		
lapply(list.cels,cels.get)
files.geo <- paste(dir.geo,list.files(path=dir.geo,pattern = ".gz"),sep="/")
# gunzip the compressed data files:
fun.gunzip <- function(x)
			R.utils::gunzip(filename = x,
			overwrite=TRUE,
			remove=FALSE)
lapply(files.geo,fun.gunzip)		

celpath1 <- "/var/folders/88/bzknwr5x3nd1tdlsgq3s0lnc0000gn/T/Rtmp9BL4OX/GSE63059/GSM1539076_NB4-VEC_HTA-2_0_.CEL"
celpath2 <- "/var/folders/88/bzknwr5x3nd1tdlsgq3s0lnc0000gn/T/Rtmp9BL4OX/GSE63059/GSM1539077_NB4-TETON.377_HTA-2_0_.CEL"

GCSs.single <- GCSscore(celFile1 = celpath1, celFile2 = celpath2)
class(GCSs.single)[1]

# convert GCSscore single-run from ExpressionSet to data.table:
GCSs.pinkd <-data.table::as.data.table(cbind(GCSs.single@featureData@data,GCSs.single@assayData[["exprs"]]))
# show all column names included in the output:
colnames(GCSs.pinkd)
write.table(GCSs.pinkd, "GSE63059_pinkd_sscore.csv", sep='\t', row.names=F, quote=F)
GCSs.pinkd_sign <- subset(GCSs.pinkd, abs(GCSs.pinkd$Sscore) >= 3)

####################
celpath1 <- "/var/folders/88/bzknwr5x3nd1tdlsgq3s0lnc0000gn/T//Rtmp9BL4OX/GSE63059/GSM1539078_NB4-DMSO_HTA-2_0_.CEL"
celpath2 <- "/var/folders/88/bzknwr5x3nd1tdlsgq3s0lnc0000gn/T//Rtmp9BL4OX/GSE63059/GSM1539079_NB4-ATRA_HTA-2_0_.CEL"

GCSs.single <- GCSscore(celFile1 = celpath1, celFile2 = celpath2)
class(GCSs.single)[1]

# convert GCSscore single-run from ExpressionSet to data.table:
GCSs.atra <-data.table::as.data.table(cbind(GCSs.single@featureData@data,GCSs.single@assayData[["exprs"]]))
# show all column names included in the output:
colnames(GCSs.atra)

GCSs.atra_sign <- subset(GCSs.atra, abs(GCSs.atra$Sscore) >= 3)
write.table(GCSs.atra, "GSE63059_atra_sscore.csv", sep='\t', row.names=F, quote=F)

#Two groups of samples
ph@data[ ,2] = c("control","control","control","mutant","mutant","mutant")
colnames(ph@data)[2]="source"
ph@data

#set sample groups
groups = ph@data$source

f = factor(groups,levels=c("control","mutant"))


design = model.matrix(~ 0 + f)
colnames(design) = c("control","mutant")

#calculate the mean expression levels using the lmFit() method. 
data.fit = lmFit(data.matrix,design)

data.fit$coefficients[1:10,]

#define groups to compare 
contrast.matrix = makeContrasts(mutant-control,levels=design)
data.fit.con = contrasts.fit(data.fit,contrast.matrix)

#run the t-test by using the eBayes()

data.fit.eb = eBayes(data.fit.con)

names(data.fit.eb)

#log fold changes
data.fit.eb$coefficients[1:10,]


#t-statistic 
data.fit.eb$t[1:10,]

#p-value
data.fit.eb$p.value[1:10,]


#Three groups of samples
ph@data[ ,2] = c("control","control","control","drug","drug","drug","exercise","exercise","exercise")
colnames(ph@data)[2]="source"
ph@data

f = factor(groups,levels = c("control","drug","exercise"))

#define design 
design = model.matrix(~ 0 + f) 
colnames(design)=c("control","drug","exercise")
design
data.fit = lmFit(data.rma,design)

#compare groups
contrast.matrix = makeContrasts(exercise-control,drug-exercise,drug-control,levels=design)
data.fit.con = contrasts.fit(data.fit,contrast.matrix)
data.fit.eb = eBayes(data.fit.con)


#Log fold changes 

data.fit.eb$coefficients[1:10,]
data.fit.eb$t[1:10,]
data.fit.eb$p.value[1:10,]
data.fit.eb$F[1:10]
data.fit.eb$F.p.value[1:10]


#Paired data
ph@data[ ,2] = c("before","treated","before","treated","before","treated")
colnames(ph@data)[2]="Treatment"
ph@data[ ,3] = c("patient1","patient1","patient2","patient2","patient3","patient3")
colnames(ph@data)[3]="Patient"
ph@data


groupsP = ph@data$Patient 
groupsT = ph@data$Treatment
fp = factor(groupsP,levels=c("patient1","patient2","patient3"))
ft = factor(groupsT,levels=c("before","treated"))


paired.design = model.matrix(~ fp + ft)
colnames(paired.design)=c("Intercept","Patient2vsPatient1","Patient3versusPatient1","AftervsBefore")
data.fit = lmFit(data.rma,paired.design)
data.fit$coefficients

data.fit.eb = eBayes(data.fit)
data.fit.eb$p.value[,1:20]


options(digits=2)
tab = topTable(data.fit.eb,coef=2,number=200,adjust.method=”BH”)

topgenes = tab[tab[, "adj.P.Val"] < 0.001, ]
dim(topgenes)


#Up- and down-regulated genes
topups = topgenes[topgenes[, "logFC"] > 1, ]
dim(topups)
topdowns = topgenes[topgenes[, "logFC"] < -1, ]
dim(topdowns)


DEresults = decideTests(data.fit.eb,method='global',adjust.method="BH",p.value=0.05,lfc=1)
DEresults[1:10,]


