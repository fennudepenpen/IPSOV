setwd("C:/Users/shen/OneDrive/OV")

options(stringsAsFactors = F)

library(readxl)
covar <- read_excel("clinical_covariates.xlsx")
table(covar$platform)
table(covar$Study,covar$platform)

### U133A
covar1=covar[covar$platform=="[HG-U133A] Affymetrix Human Genome U133A Array",]
table(covar1$Study)

study = names(table(covar1$Study))

geneMatrix1 = NULL
for (i in 1:length(study)){
  load(paste0("./GEO/",study[i],".RData"))
  if (study[i] %in% c("GSE3149","GSE23554")) geneMatrix = log(geneMatrix,base = 2)
  if (study[i] %in% c("GSE26712")) geneMatrix = log(exp(geneMatrix),base = 2)
  
  geneMatrix1 = cbind(geneMatrix1,geneMatrix)
}
hist(geneMatrix1[1,])

p = intersect(covar1$ID,colnames(geneMatrix1))
geneMatrix = geneMatrix1[,match(covar1$ID,colnames(geneMatrix1))]
# adjust batch effects
require(sva)
geneMatrix = ComBat(geneMatrix,covar1$Study)

save(covar1,ano,geneMatrix,file="HG-U133A.RData")


# U133A2.0
covar1=covar[covar$platform=="[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array",]
table(covar1$Study)

# covar1 = subset(covar1,Study!="GSE19829")

study = names(table(covar1$Study))

geneMatrix1 = NULL
for (i in 1:5){
  if(study[i] != c("GSE19829")) load(paste0("./GEO/",study[i],".RData"))
  if(study[i] %in% c("GSE19829")) load(paste0("./GEO/GSE19829-GPL570.RData"))
  if (study[i] %in% c("GSE18520")) geneMatrix = log(geneMatrix+1,base = 2)
  if (study[i] %in% c("GSE30161")) geneMatrix = log(exp(geneMatrix),base = 2)
  
  geneMatrix1 = cbind(geneMatrix1,geneMatrix)
}

load("./GEO/GSE9891.RData")
geneMatrix1 = geneMatrix1[match(rownames(geneMatrix),rownames(geneMatrix1)),]
geneMatrix = cbind(geneMatrix1,geneMatrix)

hist(geneMatrix[1,])

p = intersect(covar1$ID,colnames(geneMatrix))
geneMatrix = geneMatrix[,match(covar1$ID,colnames(geneMatrix))]

require(sva)
geneMatrix = ComBat(geneMatrix,covar1$Study)

ano = ano[match(rownames(geneMatrix),ano$ID),]

save(covar1,ano,geneMatrix,file="[n=632]U133A2.0.RData")


# Agilent-014850 Whole Human Genome Microarray 4x44K G4112F (Probe Name version)
covar1=covar[covar$platform=="Agilent-014850 Whole Human Genome Microarray 4x44K G4112F (Probe Name version)",]
table(covar1$Study)

study = names(table(covar1$Study))

geneMatrix1 = NULL
for (i in 1:length(study)){
  load(paste0("./GEO/",study[i],".RData"))
  if (study[i] %in% c("GSE32062")) geneMatrix = geneMatrix[match(rownames(geneMatrix1),rownames(geneMatrix)),]
  
  geneMatrix1 = cbind(geneMatrix1,geneMatrix)
}
hist(geneMatrix1[1,])

p = intersect(covar1$ID,colnames(geneMatrix1))
geneMatrix = geneMatrix1[,match(covar1$ID,colnames(geneMatrix1))]
require(sva)
geneMatrix = ComBat(geneMatrix,covar1$Study)

ano = ano[match(rownames(geneMatrix),ano$ID),]

save(covar1,ano,geneMatrix,file="[n=637]Agilent-014850.RData")

# Operon human v3 ~35K 70-mer two-color oligonucleotide microarrays
covar1=covar[covar$platform=="Operon human v3 ~35K 70-mer two-color oligonucleotide microarrays.",]
load("./GEO/GSE13876.RData")
geneMatrix = geneMatrix[,match(covar1$ID,colnames(geneMatrix))]

save(covar1,ano,geneMatrix,file="[n=415]Operon human v3 ~35K.RData")

# ABI Human Genome Survey Microarray Version 2
covar1=covar[covar$platform=="ABI Human Genome Survey Microarray Version 2",]
load("./GEO/GSE49997.RData")
ano = ano[match(rownames(geneMatrix),ano$ID),]
geneMatrix = geneMatrix[,match(covar1$ID,colnames(geneMatrix))]

save(covar1,ano,geneMatrix,file="[n=204]ABI Human Genome.RData")


# TCGA
covariate_TCGA <- read_csv("covariate_TCGA.csv")
exp_TCGA <- read_delim("OV__gene_Array__tissueTypeAll__20180128102820.txt", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)
a1 = exp_TCGA$gene_id
exp_TCGA = as.matrix(exp_TCGA[,-1])
rownames(exp_TCGA) = a1


exp_TCGA = exp_TCGA[,substr(colnames(exp_TCGA),14,15)=="01"]

exp_TCGA = ComBat(exp_TCGA,substr(colnames(exp_TCGA),6,7))

p = intersect(covariate_TCGA$bcr_patient_barcode,substr(colnames(exp_TCGA),1,12))

covariate_TCGA = covariate_TCGA[match(p,covariate_TCGA$bcr_patient_barcode),]
colnames(exp_TCGA) = substr(colnames(exp_TCGA),1,12)
exp_TCGA = exp_TCGA[,match(p,colnames(exp_TCGA))]

save(exp_TCGA,covariate_TCGA,file="[n=539]TCGA_OV_clean.RData")

