
# This script does QTL mapping 
rm(list=ls()); # Delete files
cat("\014");
source('genohack.R');
#library(mpMap)

runmpMap = function()
{
  #mpMap data
  # 
  # nind=720;
  # pedm=cbind(1: (8+nind), rep(0, 8+nind), rep(0, 8+nind), c(rep(0, 8), rep(1, nind)));
  # 
  # write.table(rils[1:10,], file='rils1.csv', sep=',', row.names = FALSE);
  # write.table(founders, file='founders.csv', sep=',', row.names = FALSE );
  # write.table(ped, file='ped.csv', sep=',', row.names = FALSE);
  # write.table(rils, file='rils.csv', sep=',', row.names = FALSE);
  # 
  # break();
  # rils = read.table('rils.csv', header=T,sep=',');
  # founders = read.table('founders.csv', header=T,sep=',');
  # ped = read.table('ped.csv', header=T, sep=',');
  # mapdata = read.table('wheat_geno_coordinates.csv', header=TRUE, sep=',');
  # 
  # 
  # o <- mpcross(finals=rils, founders=founders, fid=levels(founders$id),  id=levels(rils$id), 
  #                  pedigree=pedm);
  # 
  # mpprob(o, program='happy');
  # o1 <- read.mpcross(finalfile='rils.csv', founderfile='founders.csv', pedfile='ped.csv', 
  #                   phenofile='wheat_pheno.csv', mapfile = 'wheat_geno_coordinates.csv');
  # 
  # 
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
  y[which(x==alleles[2])] <- 2;
  y[which(x==alleles[3])] <- 1;
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


# phenodata
phenodata = read.table('wheat_pheno.csv', header=TRUE, sep=',');
# Founders
founders_raw = read.table('founders.geno.csv', header=T, sep=',');
#rils
rils_raw = read.table('wheat_geno.csv', header=T, sep=',');
# map
map_raw = read.table('wheat_geno_coordinates.csv', sep=',', header=TRUE);


# Process founders
colnames(founders_raw);
clf = ncol(founders_raw); 
rf = nrow(founders_raw);# Cols dataframe

# Search for markers with non-polymorphic alleles
flag=as.numeric(sapply(rownames(founders_raw), 
            function(x) searchPolym(as.numeric(founders_raw[x,9:clf]))));
i = which(flag == 0);
# List of non-poly alleles
nonpolymarkers = founders_raw$X13074[i];



# Convert 
founders_raw = founders_raw[-i, ];
F <- apply(founders_raw[,9:clf],2, convert.snp1);
founders = data.frame(colnames(founders_raw)[9:clf], t(F));
colnames(founders) = c('genotype', as.vector(founders_raw$X));

# map x founder

mapfounder = merge(map_raw, founders_raw, by.x = 'Marker', by.y='X13074');

# Rils
# delete poly markers
colnames(rils_raw);
rils_raw_clean = merge(mapfounder[,1:2], rils_raw, by.x='Marker', by.y='original.order');
clr = ncol(rils_raw_clean); 
rr = nrow(rils_raw_clean);# Cols dataframe
R <- apply(rils_raw_clean[,4:clr],2, convert.snp2);
v = as.vector(rils_raw[2, 3:(clr-1)]);
rils = data.frame(t(v), t(R));
colnames(rils) = c('genotype', as.vector(rils_raw_clean$Marker));

rildata = merge(phenodata[,c(1,17)], rils, by.x='genotype', by.y='genotype');
write.table(rildata, file='rildisease.ped', sep=' ', quote = FALSE, row.names = FALSE, row.names = FALSE);
write.table(rils, file='rildisease1.ped', sep=' ', quote = FALSE, row.names = FALSE, row.names = FALSE);
g2a(mapfounder, ".alleles");
