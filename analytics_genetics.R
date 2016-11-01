
# This script prepares data for QTL mapping
rm(list=ls()); # Delete files
cat("\014");
source('genohack.R');
#ibrary('mpMap')
library(rrBLUP)

# Convert SNPs
convert.snp <- function(x) {
  #convert to {-1,0,1,NA}
  alleles <- c(0, 2, 1); # 0=AA, 2=BB
  y <- rep(NA,length(x));
  #print(alleles);
  y[which(x==alleles[1])] <- 'AA';
  y[which(x==alleles[2])] <- 'BB';
  y[which(x==alleles[3])] <- 'AB';
  return(y)
}


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


process_founders = function(datafounder, mapdata, kw)
{
  datafounder= founders_raw
  mapdata = map_raw
  f = list();
  datafounder = merge(mapdata[,1:2], datafounder, by.y='X13074', by.x='Marker');
  datafounder = datafounder[,c(3,4,2,6:8,1, 9:17)];
  colnames(datafounder)[3] = 'pos';
  colnames(datafounder)[7] = 'X13074';
  
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
  foundersV = data.frame(datafounder$X13074, F);
  colnames(foundersV)[1] = 'Marker';
  f$foundersH = data.frame(colnames(datafounder)[9:clf], t(F));
  colnames(f$foundersH) = c('genotype', as.vector(datafounder$X));
  f$mapfounder = merge(mapdata, foundersV, by.x = 'Marker', by.y='Marker');
  f$mapfounder = data.frame(f$mapfounder, 
                            chromosome= sapply(f$mapfounder$Chr, function(x) paste('chr', x, sep='')));
  return(f);
  
}

process_rils = function(dataril, datamap)
{
  f = list()
  rils = dataril
  #rils = rils_raw
  #datamap=map_raw
  
  
  rils_raw_clean = merge(datamap[,1:2], rils, by.x='Marker', by.y='original.order');
  clr = ncol(rils_raw_clean); 
  rr = nrow(rils_raw_clean);# Cols dataframe
  f$R <- apply(rils_raw_clean[,4:clr],2, convert.snp2);
  #v = as.vector(rils_raw_clean[2, 3:(clr-1)]);
  markername = as.vector(rils_raw_clean[1:rr, 'Marker']);
  genotypename = as.vector(rils[1, 3:(clr-1)]);
  genotypename = unlist(t(genotypename));
  
  f$rildataV = data.frame(markername, f$R);
  colnames(f$rildataV) = c('Marker', genotypename);
  rownames(f$rildataV) = 1:nrow(f$rildataV);
  f$rildataH = data.frame(genotypename, t(f$R));
  colnames(f$rildataH) = c('Marker', markername);
  rownames(f$rildataH) = 1:nrow(f$rildataH);
  return(f);
  
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


check_PStruct = function(data)
{
  
  A <- A.mat(data);
  eig.result <- eigen(A)
  lambda <- eig.result$values
  return(lambda);
 
}


get_MAGIC = function(phenofile, foundersgenofile, rildgenofile, mapfile)
{
  #phenodata
  f = list()
  write.table(f$phenodata, file='wheat.phenotype', sep='\t', quote = FALSE, row.names = FALSE);
  write.table(f$founders_raw,file='founders.geno2.txt', sep=' ', quote = FALSE, row.names = FALSE);
  #g2a(f$mapfounder, ".wheat.alleles");
  return(f);
}


writefiles = function(mapwheat, founder, ril, phenodata, founderraw, pedidata, mapdata)
{
  
  write.table(mapwheat,file='map.wheat.txt', sep='\t', 
              quote = FALSE, row.names = FALSE);
  write.table(founder, file='founderfile.txt', sep='\t');
  write.table(ril, file='finalfile.txt', sep='\t');
  
  write.table(phenodata, file='wheat.phenotype', sep='\t', quote = FALSE);
  write.table(founderraw, file='founders.geno2.txt', sep=' ', quote = FALSE, row.names = FALSE);
  write.table(pedidata, file='pedfile.txt', sep='\t');
  write.table(mapdata, file='mapfile.txt', sep='\t');
}

create_pedigree = function(peddata, rildata)
{
  
  y = rildata[1:2,3:722];
  y = t(y)
  rownames(y) = 1:nrow(y);
  yy = merge(y,peddata, by.x='2', by.y='id');
  colnames(yy)[1:2] = c('Marker', 'id' );
  return(yy);
}



transform_wheatgeno = function()
{
  
  wg2 = read.table('wheat_geno2.csv', sep=',', header=F);
  wg2a = wg2[, -c(1,2,4:11)];
  clr = ncol(rils_raw);
  clr1 = ncol(wg2a);
  n = as.vector(rils_raw[1, 3:(clr)]);
  n = unlist(t(n));
  
  genotypename = as.vector(rils_raw[2, 3:(clr)]);
  genotypename = unlist(t(genotypename));
  
  genotypename1 = as.vector(wg2a[1, 2:(clr1)]);
  genotypename1 = unlist(t(genotypename1));
  
  markername = as.character(wg2a[2:nrow(wg2a), 1]);
  
  
  i = match(genotypename1, genotypename);
  mn = c('',n[i]);
  colname = colnames(rils_raw)[3:(clr)]
  cn = c('',colname[i])
  
  R <- apply(wg2a[2:18602, 2:644],2, convert.snp);
  
  markers = cbind(markername,R);
  m = rbind(cn, mn, c('',genotypename1), markers);
  
  write.table(m, file='wheat_geno_alt.csv', sep=',', quote=F, row.names=F)
}


#f = get_MAGIC('wheat_pheno.csv', 'founders.geno.csv','wheat_geno.csv','wheat_geno_coordinates.csv'); 
phenofile = 'w8_pheno.csv';
foundersgenofile = 'C:/Anyela/repo/senescence_disease/founders.geno.csv';
rildgenofile = 'C:/Anyela/repo/senescence_disease/wheat_geno.csv';
mapfile = 'wheat_geno_coordinates.csv';


pheno_raw = read.table(phenofile, header=TRUE, sep=',');
# foundersdata
founders_raw = read.table(foundersgenofile, header=T, sep=',');
#rilsdata
rils_raw = read.table(rildgenofile, header=T, sep=',');
# map
map_raw = read.table(mapfile, sep=',', header=TRUE);
# pedigree
pedigree_raw = read.table('pedigree_table.csv',sep=',', header=TRUE);

parents = process_founders(founders_raw, map_raw, 'X13074');
rils = process_rils(rils_raw, map_raw);
phenodata = create_phenodata(pheno_raw, rils_raw);
peddata = create_pedigree(pedigree_raw, rils_raw);

peddata$id = as.numeric(peddata$id);
peddata$mother = as.numeric(peddata$mother);
peddata$father = as.numeric(peddata$father);
k = match(parents$mapfounder$Marker, colnames(rils$rildataH)[2:ncol(rils$rildataH)]);

writefiles(parents$mapfounder[,c('Marker', 'chromosome', 'cM' )], parents$mapfounder[,c(1,5:12)],
                      rils$rildataH[,c(1,k+1)], phenodata, founders_raw, peddata[,2:4], map_raw[,c(1,3,2)]);

break();
# Check for population structure
p = check_PStruct(t(rils$R));
png('PCA_PS.png', width = 1600, height = 1600,res=200);
plot(p/sum(p),ylab="Fraction Explained", xlab='ordered  eigenvalues')
dev.off()

# Create mpcrop object
#cross = read.mpcross(founderfile='founderfile1.txt', finalfile='finalfile.txt', 
                     #pedfile='pedfile.txt', mapfile='mapfile.txt', phenofile = 'wheat.phenotype');

