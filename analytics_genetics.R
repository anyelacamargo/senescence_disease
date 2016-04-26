
# This script prepares data for QTL mapping
rm(list=ls()); # Delete files
cat("\014");
source('genohack.R');
#library(mpMap)

# Convert SNPs
convert.snp1 <- function(x) {
  #convert to {-1,0,1,NA}
  alleles <- c(0, 2); # 0=AA, 2=BB
  y <- rep(NA,length(x));
  #print(alleles);
  y[which(x==alleles[1])] <- 0;
  y[which(x==alleles[2])] <- 1;
  return(y)
}

# Convert SNPs
convert.snp2 <- function(x) {
  #print(x)
  #convert to {-1,0,1,NA}
  alleles <- c('AA', 'BB', 'AB'); # 0=AA, 2=BB
  y <- rep(NA,length(x));
  #print(alleles);
  y[which(x==alleles[1])] <- 0;
  y[which(x==alleles[2])] <- 1;
  y[which(x==alleles[3])] <- 0;
  #break();
  return(y)
}

searchPolym = function(data)
{
  
  profile = data;
  alleles = unique(profile);
  if(length(alleles) == 1) 
  {
    return(0);
  }
  else {return(1)};
}


process_founders = function(datafounder, mapdata, kw){

  
  colnames(datafounder);
  clf = ncol(datafounder); 
  rf = nrow(datafounder);# Cols dataframe
  # Search for markers with non-polymorphic alleles
  flag=as.numeric(sapply(rownames(datafounder), 
                         function(x) searchPolym(as.numeric(datafounder[x,9:clf]))));
  i = which(flag == 0);
  # List of non-poly alleles
  nonpolymarkers = datafounder$X13074[i];
  
  # Convert 
  datafounder = datafounder[-i, ];
  F <- apply(datafounder[,9:clf],2, convert.snp1);
  founders = data.frame(colnames(datafounder)[9:clf], t(F));
  colnames(founders) = c('genotype', as.vector(datafounder$X));
  mapfounder = merge(mapdata, datafounder, by.x = 'Marker', by.y=kw);
  return(mapfounder);
  
}

process_rils = function(data, map_raw)
{
  rils  = data;  
  colnames(rils);
  rils_raw_clean = merge(mapfounder[,1:2], rils, by.x='Marker', by.y='original.order');
  clr = ncol(rils_raw_clean); 
  rr = nrow(rils_raw_clean);# Cols dataframe
  R <- apply(rils_raw_clean[,4:clr],2, convert.snp2);
  v = as.vector(rils[2, 3:(clr-1)]);
  rils = data.frame(t(v), t(R));
  colnames(rils) = c('genotype', as.vector(rils_raw_clean$Marker));
  return(rils);
  
}

create_phenodata = function(phenotable, riltable)
{
  
  c = ncol(riltable);
  g = t(riltable[1:2,3:c]);
  k = merge(g, phenotable, by.x='2', by.y='genotype');
  k = k[,-1];
  colnames(k)[1] = 'SUBJECT.NAME';
  return(k);
}

create_pedfile = function()
{
  
  rildata = merge(phenodata[,c(1,17)], rils, by.x='genotype', by.y='genotype');
  write.table(rildata, file='rildisease.ped', sep=' ', quote = FALSE, row.names = FALSE, col.names = FALSE);
  write.table(rils, file='rildisease1.ped', sep=' ', quote = FALSE, row.names = FALSE, col.names = FALSE);
}


get_MAGIC = function(phenofile, foundersgenofile, rildgenofile, mapfile)
{
  #phenodata
  pheno_raw = read.table(phenofile, header=TRUE, sep=',');
  # foundersdata
  founders_raw = read.table(foundersgenofile, header=T, sep=',');
  #rilsdata
  rils_raw = read.table(rildgenofile, header=T, sep=',');
  # map
  map_raw = read.table(mapfile, sep=',', header=TRUE);
  # map x founder
  mapfounder = process_founders(founders_raw, map_raw, 'X13074');
    #create allele file
  phenodata = create_phenodata(pheno_raw, rils_raw);
  mapfounder = data.frame(mapfounder, chromosome= sapply(mapfounder$Chr, function(x) paste('chr', x, sep='')));
  mapfounder = mapfounder[,c(1,20,3,2, 6:19)];
  write.table(mapfounder[,c(1:3)],file='map.wheat.txt', sep='\t', quote = FALSE, row.names = FALSE);
  write.table(phenodata,file='wheat.phenotype', sep='\t', quote = FALSE, row.names = FALSE);
  g2a(mapfounder, ".wheat.alleles");
  return(mapfounder);
}


mapfounder = get_MAGIC('wheat_pheno.csv', 'founders.geno.csv','wheat_geno.csv','wheat_geno_coordinates.csv'); 


