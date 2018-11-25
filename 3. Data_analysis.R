setwd("F:/OneDrive/OV")  # home
options(stringsAsFactors = F)
look = function(x,n=5) print(x[1:n,1:n])

source("D:/Dropbox/code.R")
source("D:/Dropbox/Rcode/survivalROCCI.R")
source("D:/Dropbox/Rcode/enrichment.R")
source("D:/Dropbox/Rcode/ggsurvplot.R")
source("D:/Dropbox/Rcode/ggsurvplot_adjust.R")
require(survival)
library(readxl)
require(GSVA)


immune = read_excel("../Database/Immune genes/Immune_Geneappend3.xls")
pathway = names(table(immune$Category))

gset = list()
for (i in 1:length(pathway)){
  gset[[i]] = subset(immune,Category==pathway[i])$Symbol
  
}
names(gset) = pathway

load("./Data_clean/[n=520new]TCGA_OV_HU133A.RData")
covar1$Debulk = ifelse(covar1$residual_disease_largest_nodule%in%c("No Macroscopic disease"),1,
                       ifelse(covar1$residual_disease_largest_nodule%in%c("11-20 mm",">20 mm","1-10 mm"),0,NA))
coxph(Surv(os, death)~Debulk,data=covar1)


gene = intersect(ano$`Gene Symbol`,immune$Symbol)
probe = ano[ano$`Gene Symbol`%in%gene,"ID"]
geneMatrix = t(geneMatrix[rownames(geneMatrix)%in%probe,])

## choose probes
rlt_133A = NULL
for (i in 1:ncol(geneMatrix)){
  
  cur_cox = coxph(Surv(os, death)~geneMatrix[,i]+age+cstage+grade+Debulk, data = covar1)
  cur_rlt = c(gene=colnames(geneMatrix)[i],hr=exp(coef(cur_cox)[1]), p_univ=summary(cur_cox)$coef[1,ncol(summary(cur_cox)$coef)])  

  rlt_133A = rbind(rlt_133A,cur_rlt)
  
}
rlt_133A = as.data.frame(rlt_133A);rownames(rlt_133A) = 1:nrow(rlt_133A)

probe = rlt_133A[rlt_133A$p_univ<0.05,]$gene
geneMatrix = geneMatrix[,colnames(geneMatrix)%in%probe]


colnames(geneMatrix) = ano$`Gene Symbol`[match(colnames(geneMatrix),ano$ID)]


genelist = ano$`Gene Symbol`[match(probe,ano$ID)]
geneano = data.frame(genelist,probe)


## perform GSVA analysis
gsva_shen = function(dataset = geneMatrix,method = "ssgsea",gset=gset){
  if (method=="ssgsea") b = gsva(t(dataset), gset,method=method)
  if (method=="gsva")b = gsva(t(dataset), gset,method=method)$es.obs
  c = NULL
  for (i in 1:nrow(b)){
    model = coxph(Surv(os,death)~b[i,],data=covar1)
    c = rbind(c,c(rownames(b)[i],coef(summary(model))[1,2],coef(summary(model))[1,5]))
    
  }
  return(list(c,b))
}


rlt = NULL
score = NULL


load("./Data_clean/[n=520new]TCGA_OV_HU133A.RData")
geneMatrix = t(geneMatrix[rownames(geneMatrix)%in%probe,])
colnames(geneMatrix) = ano$`Gene Symbol`[match(colnames(geneMatrix),ano$ID)]
rlt1 = gsva_shen(geneMatrix)
rlt = cbind(rlt,rlt1[[1]][,2:3]);rownames(rlt) = rlt1[[1]][,1]
score = cbind(score,rlt1[[2]])

load("./Data_clean/[n=409]HG-U133A.RData")
geneMatrix = t(geneMatrix[rownames(geneMatrix)%in%probe,])
colnames(geneMatrix) = ano$`Gene Symbol`[match(colnames(geneMatrix),ano$ID)]
rlt1 = gsva_shen(geneMatrix)
rlt = cbind(rlt,rlt1[[1]][,2:3]);rownames(rlt) = rlt1[[1]][,1]
score = cbind(score,rlt1[[2]])


load("./Data_clean/[n=632]U133A2.0.RData")
geneMatrix = t(geneMatrix[rownames(geneMatrix)%in%probe,])
colnames(geneMatrix) = ano$`Gene Symbol`[match(colnames(geneMatrix),ano$ID)]
rlt1 = gsva_shen(geneMatrix)
rlt = cbind(rlt,rlt1[[1]][,2:3]);rownames(rlt) = rlt1[[1]][,1]
score = cbind(score,rlt1[[2]])

## choose genes with same name
gene = names(table(as.character(geneano$genelist)))
ident_gene = function(gene,geneMatrix){
  d = NULL;dc=NULL
  for (i in 1:length(gene)){
    d1 = as.matrix(geneMatrix[,colnames(geneMatrix)==gene[i]] )
    if (ncol(d1) == 0) next
    if (ncol(d1) == 1) d = cbind(d,d1)
    if (ncol(d1) > 1) {
      m = apply(d1,2,function(x) mean(x,na.rm = T))
      d1 = d1[,which.max(m)]
      d = cbind(d,d1)
    }
    dc = c(dc,gene[i])
  }
  colnames(d) = dc
  return(d)
}



load("./Data_clean/[n=637]Agilent-014850.RData")
rownames(geneMatrix) = ano$GENE_SYMBOL[match(rownames(geneMatrix),ano$ID)]
geneMatrix = t(geneMatrix[rownames(geneMatrix)%in%geneano$genelist,])
geneMatrix =  ident_gene(gene,geneMatrix)

rlt1 = gsva_shen(geneMatrix)
rlt = cbind(rlt,rlt1[[1]][,2:3]);rownames(rlt) = rlt1[[1]][,1]
score = cbind(score,rlt1[[2]])

load("./Data_clean/[n=415]Operon human v3 ~35K.RData")
rownames(geneMatrix) = ano$`Gene Symbol`[match(rownames(geneMatrix),ano$ID)]
geneMatrix = t(geneMatrix[rownames(geneMatrix)%in%geneano$genelist,])
geneMatrix =  ident_gene(gene,geneMatrix)

rlt1 = gsva_shen(geneMatrix)
rlt = cbind(rlt,rlt1[[1]][,2:3]);rownames(rlt) = rlt1[[1]][,1]
score = cbind(score,rlt1[[2]])

load("./Data_clean/[n=204]ABI Human Genome.RData")
rownames(geneMatrix) = ano$`Gene Symbol`[match(rownames(geneMatrix),ano$ID)]
geneMatrix = t(geneMatrix[rownames(geneMatrix)%in%geneano$genelist,])
geneMatrix =  ident_gene(gene,geneMatrix)
rlt1 = gsva_shen(geneMatrix)
# impute one immune category
rlt1[[2]] = rbind(rlt1[[2]][1:7,],rep(0,204),rlt1[[2]][8:14,])
score = cbind(score,rlt1[[2]])


load("covar_TCGAnewdebulk.RData")
covar = covar[covar$os>0,]
score = score[,match(covar$ID,colnames(score))]
covar  = covar[match(colnames(score),covar$ID),]

score = t(score)

c = NULL
for (i in 1:ncol(score)){
  model = coxph(Surv(os,death)~score[,i],data=covar)
  c = rbind(c,c(colnames(score)[i],coef(summary(model))[1,2],coef(summary(model))[1,5]))
  
}
c

model = coxph(Surv(os,death)~score,data=covar)
id_s = 1:ncol(score)

model = coxph(Surv(os,death)~score[,id_s],data=covar)


index = 0 ; for(i in 1:length(id_s)) index = index+coef(model)[i]*score[,id_s[i]]

covar$index = index

summary(coxph(Surv(os,death)~index,data=covar))
summary(coxph(Surv(os,death)~mcut(index,0.5),data=covar))

summary(coxph(Surv(os,death)~index+stage+grade+Debulk,data=covar))

table(covar$platform)
summary(coxph(Surv(os,death)~index,data=subset(covar,platform=="TCGA")))
summary(coxph(Surv(os,death)~index,data=subset(covar,platform=="[HG-U133A] Affymetrix Human Genome U133A Array")))
summary(coxph(Surv(os,death)~index,data=subset(covar,platform=="[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array")))
summary(coxph(Surv(os,death)~index,data=subset(covar,platform=="Agilent-014850 Whole Human Genome Microarray 4x44K G4112F (Probe Name version)")))
summary(coxph(Surv(os,death)~index,data=subset(covar,platform=="Operon human v3 ~35K 70-mer two-color oligonucleotide microarrays.")))
summary(coxph(Surv(os,death)~index,data=subset(covar,platform=="ABI Human Genome Survey Microarray Version 2")))



summary(coxph(Surv(os,death)~index+stage+grade+Debulk+as.factor(study),data=subset(covar,platform=="[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array")))

sapply(split(covar$index,covar$platform),median)

covar$index_d = ifelse(covar$index<=-0.23,1,2)

summary(coxph(Surv(os,death)~index_d,data=subset(covar,platform=="TCGA")))
summary(coxph(Surv(os,death)~index_d+age+stage+Debulk+grade+as.factor(study),data=subset(covar,platform=="TCGA")))

pdf("./plotnew/km-training.pdf",width=8,height = 8)
covar1 = subset(covar,platform=="TCGA")
ggsurvplot_adjust(data=covar1,time = "os",event = "death",marker = "index_d"
                  ,xmax = 90,timeby = 30,ystratalabs = c("Low-risk", "High-risk"),col = c(2,4),
                  xlab = "Overall survival (months)",covar = c("age","grade","stage","Debulk"))
dev.off()



summary(coxph(Surv(os,death)~index_d,data=subset(covar,platform=="[HG-U133A] Affymetrix Human Genome U133A Array")))
summary(coxph(Surv(os,death)~index_d+stage+Debulk+grade+as.factor(study),data=subset(covar,platform=="[HG-U133A] Affymetrix Human Genome U133A Array")))

pdf("./plotnew/km-v1.pdf",width=8,height = 8)
covar1 = subset(covar,platform=="[HG-U133A] Affymetrix Human Genome U133A Array")
ggsurvplot_adjust(data=covar1,time = "os",event = "death",marker = "index_d"
                  ,xmax = 180,timeby = 60,ystratalabs = c("Low-risk", "High-risk"),col = c(2,4),
                  xlab = "Overall survival (months)",covar = c("grade","stage","Debulk"))
dev.off()


summary(coxph(Surv(os,death)~index_d,data=subset(covar,platform=="[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array")))
summary(coxph(Surv(os,death)~index_d+stage+Debulk+grade+as.factor(study),data=subset(covar,platform=="[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array")))

pdf("./plotnew/km-v2.pdf",width=8,height = 8)
covar1 = subset(covar,platform=="[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array")
ggsurvplot_adjust(data=covar1,time = "os",event = "death",marker = "index_d"
                  ,xmax = 120,timeby = 24,ystratalabs = c("Low-risk", "High-risk"),col = c(2,4),
                  xlab = "Overall survival (months)",covar = c("grade","stage","Debulk"))
dev.off()

summary(coxph(Surv(os,death)~index_d,data=subset(covar,platform=="Agilent-014850 Whole Human Genome Microarray 4x44K G4112F (Probe Name version)")))
summary(coxph(Surv(os,death)~index_d+stage+Debulk+grade+as.factor(study),data=subset(covar,platform=="Agilent-014850 Whole Human Genome Microarray 4x44K G4112F (Probe Name version)")))

pdf("./plotnew/km-v3.pdf",width=8,height = 8)
covar1 = subset(covar,platform=="Agilent-014850 Whole Human Genome Microarray 4x44K G4112F (Probe Name version)")
ggsurvplot_adjust(data=covar1,time = "os",event = "death",marker = "index_d"
                  ,xmax = 120,timeby = 24,ystratalabs = c("Low-risk", "High-risk"),col = c(2,4),
                  xlab = "Overall survival (months)",covar = c("grade","stage","Debulk"))
dev.off()


summary(coxph(Surv(os,death)~index_d,data=subset(covar,platform=="Operon human v3 ~35K 70-mer two-color oligonucleotide microarrays.")))
summary(coxph(Surv(os,death)~index_d+age,data=subset(covar,platform=="Operon human v3 ~35K 70-mer two-color oligonucleotide microarrays.")))
pdf("./plotnew/km-v4.pdf",width=8,height = 8)
covar1 = subset(covar,platform=="Operon human v3 ~35K 70-mer two-color oligonucleotide microarrays.")
ggsurvplot_adjust(data=covar1,time = "os",event = "death",marker = "index_d"
                  ,xmax = 120,timeby = 24,ystratalabs = c("Low-risk", "High-risk"),col = c(2,4),
                  xlab = "Overall survival (months)",covar = c("age"))
dev.off()

summary(coxph(Surv(os,death)~index_d,data=subset(covar,platform=="ABI Human Genome Survey Microarray Version 2")))
summary(coxph(Surv(os,death)~index_d+age+stage,data=subset(covar,platform=="ABI Human Genome Survey Microarray Version 2")))
pdf("./plotnew/km-v5.pdf",width=8,height = 8)
covar1 = subset(covar,platform=="ABI Human Genome Survey Microarray Version 2");covar1$index_d = mcut(covar1$index,0.56)
ggsurvplot_adjust(data=covar1,time = "os",event = "death",marker = "index_d"
                  ,xmax = 120,timeby = 12,ystratalabs = c("Low-risk", "High-risk"),col = c(2,4),
                  xlab = "Overall survival (months)",covar = c("age"))
dev.off()




## RMS ratio
require(survRM2)

covar1$index_d1 = -covar1$index_d+2;covar1$index_d1[is.na(covar1$index_d1)]=1
table(covar1$index_d)
rmst2(covar1$os,covar1$death,covar1$index_d1,tau = NULL)


## meta analysis for overall effect
meta <- read_excel("./meta.xlsx")
require(meta)
m = metagen(TE = log(meta$HR),seTE = meta$se,sm = "HR",n.e = meta$N)
# study.results: show results for individual studies

pdf("./plotnew/meta.pdf")
forest(m,studlab = meta$Study,comb.random = F,lty.fixed=0,
       leftcols=c("studlab","n.e","effect","ci"),leftlabs = c("Study","N","HR","95% CI"),
       rightcols = F,
       print.tau2 = F,col.diamond.fixed = "red",col.square="black" )
dev.off()

## subgroup analysis
summary(coxph(Surv(os,death)~index,data=subset(covar,age<=60)))
summary(coxph(Surv(os,death)~index,data=subset(covar,ctype==1)))
summary(coxph(Surv(os,death)~index,data=subset(covar,stage>2)))
summary(coxph(Surv(os,death)~index,data=subset(covar,Debulk==1)))
summary(coxph(Surv(os,death)~index,data=subset(covar,grade<3)))
summary(coxph(Surv(os,death)~index,data=subset(covar)))

meta <- read_excel("./plotnew/Stratification.xlsx",sheet = "Strat")
require(meta)
m = metagen(TE = log(meta$HR),seTE = meta$se,sm = "HR",n.e = meta$N)
# study.results: show results for individual studies

pdf("./plotnew/Strat.pdf")
forest(m,studlab = meta$Subgroup,comb.random = F,lty.fixed=0,
       leftcols=c("studlab","n.e","effect","ci"),leftlabs = c("Study","N","HR","95% CI"),
       rightcols = F,
       print.tau2 = F,col.diamond.fixed = "red",col.square="black" )
dev.off()



#### sub pathway analysis
cate = "Cytokine_Receptors"
gset1 = list(gset[[cate]])

immune1 = immune[immune$Symbol%in%gene,]
gene1 = immune1[immune1$Category==cate,]$Symbol
probe1 = geneano[geneano$genelist%in%gene1,]$probe


score = NULL

load("./Data_clean/[n=520new]TCGA_OV_HU133A.RData")
geneMatrix = t(geneMatrix[rownames(geneMatrix)%in%probe1,])
colnames(geneMatrix) = ano$`Gene Symbol`[match(colnames(geneMatrix),ano$ID)]
rlt1 = gsva_shen(geneMatrix,gset = gset1);names(rlt1[[2]]) = rownames(geneMatrix)
score = c(score,rlt1[[2]])

load("./Data_clean/[n=409]HG-U133A.RData")
geneMatrix = t(geneMatrix[rownames(geneMatrix)%in%probe1,])
colnames(geneMatrix) = ano$`Gene Symbol`[match(colnames(geneMatrix),ano$ID)]
rlt1 = gsva_shen(geneMatrix,gset = gset1);names(rlt1[[2]]) = rownames(geneMatrix)
score = c(score,rlt1[[2]])


load("./Data_clean/[n=632]U133A2.0.RData")
geneMatrix = t(geneMatrix[rownames(geneMatrix)%in%probe1,])
colnames(geneMatrix) = ano$`Gene Symbol`[match(colnames(geneMatrix),ano$ID)]
rlt1 = gsva_shen(geneMatrix,gset = gset1);names(rlt1[[2]]) = rownames(geneMatrix)
score = c(score,rlt1[[2]])

load("./Data_clean/[n=637]Agilent-014850.RData")
rownames(geneMatrix) = ano$GENE_SYMBOL[match(rownames(geneMatrix),ano$ID)]
geneMatrix = t(geneMatrix[rownames(geneMatrix)%in%gene1,])
#geneMatrix =  ident_gene(gene,geneMatrix)

rlt1 = gsva_shen(geneMatrix,gset = gset1);names(rlt1[[2]]) = rownames(geneMatrix)
score = c(score,rlt1[[2]])

load("./Data_clean/[n=415]Operon human v3 ~35K.RData")
rownames(geneMatrix) = ano$`Gene Symbol`[match(rownames(geneMatrix),ano$ID)]
geneMatrix = t(geneMatrix[rownames(geneMatrix)%in%geneano$genelist,])
geneMatrix =  ident_gene(gene,geneMatrix)

rlt1 = gsva_shen(geneMatrix,gset = gset1);names(rlt1[[2]]) = rownames(geneMatrix)
score = c(score,rlt1[[2]])

load("./Data_clean/[n=204]ABI Human Genome.RData")
rownames(geneMatrix) = ano$`Gene Symbol`[match(rownames(geneMatrix),ano$ID)]
geneMatrix = t(geneMatrix[rownames(geneMatrix)%in%geneano$genelist,])
geneMatrix =  ident_gene(gene,geneMatrix)
rlt1 = gsva_shen(geneMatrix,gset = gset1);names(rlt1[[2]]) = rownames(geneMatrix)
score = c(score,rlt1[[2]])

score = score[match(covar$ID,names(score))]
summary(coxph(Surv(os,death)~score,data=covar))


covar$score = score

survfit(Surv(os,death)~mcut(score),data=covar)


# random select
hr=NULL;p=NULL
for (i in 1:10000){
  id = sample(1:nrow(covar),500)
  m = summary(coxph(Surv(os,death)~index,data=covar[id,]))
  hr[i] = m$concordance[1];p[i] = coef(m)[5]
}
hist(-log(p,10),nclass = 100,cex.axis=1.5,cex.lab=1.5,main="")
abline(v=-log(0.05,10),col="blue",lwd=2,lty=2)

hist(hr,nclass = 100,cex.axis=1.5,cex.lab=1.5,main="",xlab="C-index")
abline(v=0.624,col="blue",lwd=2,lty=2)

# pathway enrichment
kegg = keggenrich(gene)
a1 = goenrich(gene,type="BP")
a2 = goenrich(gene,type="CC",pvalueCutoff=0.5)
a3 = goenrich(gene,type="MF")


# Integrate IPSOV with clinical characteristics

summary(coxph(Surv(os,death)~index+stage+grade+Debulk,data=covar))
model = coef(coxph(Surv(os,death)~index+stage+grade+Debulk,data=covar))

index_clinical = model[1]*covar$index+model[2]*covar$stage+model[3]*covar$grade+model[4]*covar$Debulk
index_clinicalonly = model[2]*covar$stage+model[3]*covar$grade+model[4]*covar$Debulk
covar$index_clinical = index_clinical
covar$index_clinicalonly = index_clinicalonly

summary(coxph(Surv(os,death)~index_clinicalonly,data=covar))
summary(coxph(Surv(os,death)~index_clinical,data=covar))


library("compareC")
table(covar$platform)
covar1=subset(covar,platform=="Agilent-014850 Whole Human Genome Microarray 4x44K G4112F (Probe Name version)")
summary(coxph(Surv(os,death)~index_clinicalonly,data=covar1))
summary(coxph(Surv(os,death)~index_clinical,data=covar1))
covar2 = covar1[,c("os","death","index_clinical","index_clinicalonly")];covar2 = covar2[complete.cases(covar2),]
compareC(covar2$os,covar2$death,covar2$index_clinical,covar2$index_clinicalonly)

summary(coxph(Surv(os,death)~stage,data=covar))
summary(coxph(Surv(os,death)~grade,data=covar))
summary(coxph(Surv(os,death)~Debulk,data=covar))

### ano list
geneano1 = geneano[!duplicated(geneano$genelist),]
geneano2 = merge(geneano1,immune,by.x="genelist",by.y="Symbol")
dupgene = geneano2$genelist[duplicated(geneano2$genelist)]
geneano3 = geneano2[!duplicated(geneano2$genelist),]

for (i in 1:length(dupgene)){
  geneano4  = geneano2[geneano2$genelist %in% dupgene[i],]
  newcat = paste(as.character(geneano4$Category),collapse = ", ")
  geneano3$Category[geneano3$genelist==dupgene[i]]=newcat
  
}

write.csv(geneano3,"./plotnew/geneano.csv")

#### plot RMS
source("D:/Dropbox/Rcode/RMSCurve.R")
table(covar$platform)
covar1=subset(covar,platform=="[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array")

pdf("./plotnew/RMS3.pdf")
RMSCurve(os = covar1$os,death = covar1$death,marker=covar1$index,n.point=150,time=100)
par(new=T)
RMSCurve(os = covar1$os,death = covar1$death,marker=covar1$index_clinical,n.point=150,time=100,col = "red")
dev.off()

