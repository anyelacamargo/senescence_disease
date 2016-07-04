library(gdata);
library(ggplot2);
source('generic.r')


doodle = read.table('w8_doodle.csv', header=T, sep=',');
freqs = colnames(doodle)[3:dim(doodle)[2]];
freqs = as.numeric(sub("X", "", freqs));
colnames(doodle)[3:dim(doodle)[2]] = freqs;
groundtruth = read.table('w8_groundtruth.csv', header=TRUE, sep=',');

# Read table with frequencies
it = read.table('contable.csv', header=T, sep=',');
it = data.frame(it, hex=rgb(it[,2:4], max=255));
# Select intensities
intens = c(1,4,5,16,20,21,24,25,26,27,30,31,32,50,54,55,60,108,112,116);
#
intens = c(25,26,27,30,31,32,50,54,55,60,108,112,116);
# Select subdataset

idx = which(groundtruth$disease_score2 == 4);
idtagname_list = unique(groundtruth$barcode[idx]) 
# 
idtagname_list = c('W8-027112','W8-100112');

for(idtagname in idtagname_list)
{
  print(idtagname)
  ix = which(doodle$idtag == idtagname);
  scsub = doodle[ix,];
  h =  average_similar_measures(scsub, 'timestamp');
  hh = convert2numeric(h, colnames(h)[3:(ncol(h)-2)])
  hh = calibratedata(hh);
  cdatacal = get_subdata(hh, it, intens);
  
  l = which(groundtruth$barcode == idtagname);
  p = print_plot(cdatacal, 
                 as.numeric(as.character(groundtruth[l, 'flag.leaf.senescence'])), 
                 'Area', 'FLS DAS');
  pdf(paste(idtagname, '_dis.pdf',sep=''), paper='special');
  print(p);
  dev.off();
  #if(readline(idtagname) == 'q') { break; }
}

