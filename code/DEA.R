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


