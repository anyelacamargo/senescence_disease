
# This script prepares data for QTL mapping
rm(list=ls()); # Delete files
cat("\014");


stackFiles = function(patternword)
{
  files <- list.files(pattern = patternword);
  myfiles = do.call(rbind, lapply(files, 
                                  function(x) read.table(x, stringsAsFactors = TRUE, header=T)));
  return(myfiles);
  
}



createQTLTable = function(mapdata)
{
  setwd('qtl');
  qtldata = stackFiles('qtls');
  imputedata = stackFiles('imputed');
  h = merge(qtldata, mapdata, by.x='peak.SNP', by.y='Marker');
  h = h[order(h$phenotype,h$chromosome, h$cM),];
  h$logP = round(h$logP,2);
  h$genomewide.pvalue = round(h$genomewide.pvalue,2);
  h$cM = round(h$cM,2);
  accession_list = as.vector(unique(imputedata$Accession));
  accession = data.frame(matrix(1, nrow=nrow(h), ncol= length(accession_list)));
  colnames(accession) = accession_list;
  qtlfull = data.frame(h, accession);
  
  for(trait in as.vector(h$phenotype))
  {
    for(marker in as.vector(h$peak.SNP))
    {
      
      for(a in accession_list)
      {
     
        m = intersect(which(imputedata$Phenotype == trait),
                      intersect(which(imputedata$SNP == marker), which(imputedata$Accession == a)));
        
        if(length(m) > 0)
        {
          i = which(qtlfull$peak.SNP == marker);
          j = which(qtlfull$phenotype == trait);
          i = intersect(i,j);
          qtlfull[i, a] = round(imputedata[m,'mean'],2)
        }
        
      }
      
    }
  }
  
return(qtlfull);
}

get_MAGIC = function(mapfile)
{
  #phenodata
  f = list()
  # map
  f$map_raw = read.table(mapfile, sep=',', header=TRUE);
  return(f);
}


f = get_MAGIC('wheat_geno_coordinates.csv'); 
t = createQTLTable(f$map_raw);
write.table(t[,c(2,1,9,10,7,8,12:19)], file='qtltable.csv', sep=',', row.names = F);
setwd('..');