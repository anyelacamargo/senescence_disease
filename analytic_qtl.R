
# This script prepares data for QTL mapping
rm(list=ls()); # Delete files
cat("\014");
setwd("M:/anyela/repo/senescence_disease");
source('generic.R')
source('qtl/manhattan.r');
source('qtl/qvalue.r');


createQTLTable = function(mapdata)
{
  setwd('qtl');
  qtldata = stackDataSets(stackFiles('qtls'));
  imputedata = stackDataSets(stackFiles('imputed'));
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


break();
f = get_MAGIC('wheat_geno_coordinates.csv'); 
t = createQTLTable(f$map_raw);
write.table(t[,c(2,1,9,10,7,8,12:19)], file='qtltable.csv', sep=',', row.names = F);
setwd('..');

setwd('qtl');
filelist = stackFiles('scan')
u = loadMarkerFile(filelist[c(1:3,19)], TRUE, '\t', 'from');


for(cname in colnames(u)[4:14])
{
  #tiff(paste('qtl_', cname,'.tiff'),  width = 1400, height = 1080,res=200);
  manhattan(u[,c('marker','chrom', 'pos',cname)], fdr.level = 0.05, cname)
  #dev.off();
  if(readline(cname) == 'q') { break()}
}
