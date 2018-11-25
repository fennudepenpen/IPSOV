setwd("./TN")
options(stringsAsFactors = F)

require(GEOquery)
require(stringr)

# downGSE function by SShen
# download verson 2
downGSE_v2 = function(gse,platform,multiplesets = F,outname = NA,outpheno = T){

  gse1 = substr(gse,1,2)
  if(!multiplesets) url = paste0("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE",gse1,"nnn/GSE",gse,"/matrix/GSE",gse,"_series_matrix.txt.gz")
  if(multiplesets) url = paste0("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE",gse1,"nnn/GSE",gse,"/matrix/GSE",gse,"-",platform,"_series_matrix.txt.gz")
  
  download.file(url,destfile=paste0("GSE",gse,"_series_matrix.txt.gz"))
  gse2 = parseGEO(paste0("GSE",gse,"_series_matrix.txt.gz"))
  
  # extract microarray file
  geneMatrix = exprs(gse2)
  
  # extract annotation file
  ano = fData(gse2)
  ano = ano[match(rownames(geneMatrix),ano$ID),]
  
  # extract pheno file
  pheno = pData(gse2)
  pheno1 = varLabels(gse2)
  
  index1 = str_extract(pheno1,"characteristics_")
  index2 = which(!is.na(index1))
  for(i in 1:length(index2)){
    char = pheno[,index2[i]]
    index3 = which(!is.na(str_extract(char,":")))
    title = unlist(lapply(char,function(x) strsplit(x,": ")[[1]][1]))
    if (length(table(title))>1) 
    {title = paste0(names(table(title)),collapse = "/")} else
     {title = names(table(title))}
    
    pheno[,index2[i]] = unlist(lapply(char,function(x) strsplit(x,": ")[[1]][2]))
    colnames(pheno)[index2[i]] = title
  }
  
  # convert characters to numeric values
  if (is.na(outname)) outname = paste0("GSE",gse)
  write.csv(pheno,paste0(outname,".csv"),row.names = F)
  pheno = read.csv(paste0(outname,".csv"))
  
  # output platform
  print(paste("platform :",gse2$platform_id[1]))
  
 # output files
  if(!multiplesets) save(geneMatrix,ano,pheno,file = paste0(outname,".RData"))
  if(multiplesets) save(geneMatrix,ano,pheno,file = paste0(outname,"-",platform,".RData"))
  
  # remove files
  file.remove(paste0("GSE",gse,"_series_matrix.txt.gz"))
  if(!outpheno) file.remove(paste0(outname,".csv"))

  
  
}

downGSE_v2(gse = "14764",multiplesets=F , platform ="GPL96",outname = NA,outpheno = F)

# loop function for these GEO datasets
# GSE14764
# GSE23554
# GSE26712
# GSE3149
# GSE18520
# GSE19829
# GSE26193
# GSE30161
# GSE63885
# GSE9891
# GSE49997
# GSE13876
# GSE17260
# GSE32062
# GSE53963
# GSE73614
