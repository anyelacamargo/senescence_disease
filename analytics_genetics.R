# This script does QTL mapping 
rm(list=ls()); # Delete files
cat("\014");
library(mpMap)


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
  y[which(x==alleles[3])] <- 2;
  #break();
  return(y)
}

# phenodata
phenodata = read.table('wheat_pheno.csv', header=TRUE, sep=',');
# Founders
founders_raw = read.table('founders.geno.csv', header=T, sep=',');
colnames(founders_raw);
clf = ncol(founders_raw); 
rf = nrow(founders_raw);# Cols dataframe
F <- apply(founders_raw[,9:clf],2, convert.snp1);
founders = data.frame(colnames(founders_raw)[9:clf], t(F));
colnames(founders) = c('genotype', as.vector(founders_raw$X));

# Rils
rils_raw = read.table('wheat_geno.csv', header=T, sep=',');
clr = ncol(rils_raw); 
rr = nrow(rils_raw);# Cols dataframe
R <- apply(rils_raw[3:rr,3:clr],2, convert.snp2);
v = as.vector(rils_raw[2, 3:clr]);
rils = data.frame(t(v), t(R));
colnames(rils) = c('genotype', as.vector(rils_raw$original.order[3:rr]));

#Map data
mapdata = read.table('wheat_geno_coordinates.csv', header=TRUE, sep=',');


nind=720;
pedm=cbind(1: (8+nind), rep(0, 8+nind), rep(0, 8+nind), c(rep(0, 8), rep(1, nind)));

write.table(rils[1:10,], file='rils1.csv', sep=',', row.names = FALSE);
write.table(founders, file='founders.csv', sep=',', row.names = FALSE );
write.table(ped, file='ped.csv', sep=',', row.names = FALSE);
write.table(rils, file='rils.csv', sep=',', row.names = FALSE);

break();
rils = read.table('rils.csv', header=T,sep=',');
founders = read.table('founders.csv', header=T,sep=',');
ped = read.table('ped.csv', header=T, sep=',');

o <- mpcross(finals=rils, founders=founders, fid=levels(founders$id),  id=levels(rils$id), 
                 pedigree=pedm);

mpprob(o, program='happy');
o1 <- read.mpcross(finalfile='rils.csv', founderfile='founders.csv', pedfile='ped.csv', 
                  phenofile='wheat_pheno.csv', mapfile = 'wheat_geno_coordinates.csv');




