rm(list=ls()); # Delete files
source('generic.r')

# Make sure you use ggplot for 2.14
# Read metadata



metafilename = 'w8_metadata.csv'; # Read metadata
metatable = read_data(metafilename, TRUE);
groundtruthfilename = 'w8_groundtruth.csv';
groundtruthtable = read_data(groundtruthfilename, TRUE);

cname_list = list("X" = 'genotype', "row" = 'location', "X1..earLength" = 'first.ear.length', 
         "X2nd.ear.L" = "second.ear.length", "X3rd.ear.L" = 'third.ear.length', 
         "X1.fg.lf.L" = "first.fg.lf.length","X1.ear.wt" = "first.ear.weight",
         'other.ear.wt' = "other.ear.weith", 'plant.wt' = 'plant.weight', 
         'stem.ht' = 'stem.height', 'top.internodeL' = 'top.internode.length',
         'GS65' = "gs65");

for(cname in names(cname_list))
{
  groundtruthtable = change_colname(groundtruthtable, cname, cname_list[[cname]]); 
}


rname_list = data.frame(rbind(c(colname ='tiller.count', oldname = '>', newname = 11),
                   c(colname ='tiller.count', oldname = '<', newname = 10),
                   c(colname = 'flag.leaf.senescence', oldname = '*** base rot ', newname = NA),
                   c(colname = 'flag.leaf.senescence', oldname = '171?', newname = NA),
                   c(colname = 'flag.leaf.senescence', oldname = '', newname = NA),
                   c(colname = 'plant.weight', oldname = '6*** contaminant', newname = NA),
                   c(colname = 'gs55', oldname = '?', newname = NA),
                   c(colname = 'gs39', oldname = '?', newname = NA),
                   c(colname = 'first.fg.lf.length', oldname = '?', newname = NA),
                   c(colname = 'first.fg.lf.length', oldname = '', newname = NA)));

# Change row values
for(i in 1:nrow(rname_list))
{
  groundtruthtable = change_rowvalue(groundtruthtable, as.vector(rname_list$colname[i]), 
                                     as.vector(rname_list$oldname[i]),
                                     as.vector(rname_list$newname[i])); 
}


# Merge metadata and groundtruth
data_all = merge(metatable,groundtruthtable, groundtruth, 
                 by.x='barcode', by.y='barcode');

# Create a list with parameters to be used in the plot
plotname_list = list("main" = '', 'x.title' = 'DAS Flag Leaf Senescence',
                     'y.title' = 'Disease score', 'aes.x' = 'flag.leaf.senescence', 
                     'aes.y' = 'disease_score2');


# Plot per trait
for(cname in colnames(data_all)[12:26])
{
  plotname_list[['aes.x']] = cname;
  plotname_list[['x.title']] = cname; 
  p = create_figure(data_all, plotname_list);
  tname = paste(strsplit(cname, '[.]')[[1]], collapse = '');
  print(tname)
  png(paste(tname, '.png', sep=''));
  print(p);
  dev.off();
}



## Analysis of color data