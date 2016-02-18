source('generic.r')
# Read metadata

metafilename <- 'metadata.csv';
metatable <- read_data(metafilename, TRUE);
groundtruthfilename <- 'groundtruth.csv';
groundtruthtable <- read_data(groundtruthfilename, TRUE);

# Merge metadata and groundtruth
data_all = merge(metatable,groundtruthtable, groundtruth, 
                 by.x='barcode', by.y='barcode');
# Average scores
data_all = data.frame(data_all, av_diseasescore=apply(data_all[17:18],1,mean));

# Get rid of non-sense
data_all$FlagLeafSenescence_DAS[data_all$FlagLeafSenescence_DAS == '' ] = NA;
data_all$FlagLeafSenescence_DAS[data_all$FlagLeafSenescence_DAS == '*** base rot '] = NA;

i = !is.na(data_all$FlagLeafSenescence_DAS);
data_all = data_all[i,];
#data_all$FlagLeafSenescence_DAS = as.numeric(data_all$FlagLeafSenescence_DAS);

ggplot(data=data_all, aes(x=disease_score2, y= FlagLeafSenescence_DAS)) + 
   geom_line() + 
   geom_point()

p <- ggplot(data_all, aes(FlagLeafSenescence_DAS, disease_score2))
p + geom_point()

create_figure= function(data)
{
  copydata = data;
  p = ggplot(copydata, aes(x=FlagLeafSenescence_DAS, y=disease_score2, group = 1)) +        
    geom_point(shape=1) +    # Use hollow circles
    geom_smooth(method=lm) +
    xlab("FLS DAS") + ylab("Disease score") +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1)) 
  png('s_d.png');
  print(p);
  dev.off();  
}

create_figure(data_all)
dev.off()
