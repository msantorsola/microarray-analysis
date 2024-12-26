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

