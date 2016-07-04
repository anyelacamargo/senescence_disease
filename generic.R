library('ggplot2'); 
library(data.table);#
#library('RODBC');
#library('RPostgreSQL');
#library('RMySQL');



# Merge files
datatamerge_files <- function(p = 'txt')
{
 	filenames <- list.files(pattern = p);
 	c <- do.call("rbind", lapply(filenames, read.csv, header = F));
 	c <- read.table('output.csv', sep=',', header=F);
}


# Function to read data from external files
read_data = function(datafilename, h)
{
  idata <- read.table(datafilename, sep=',', header=h, na.strings = "NA", stringsAsFactors=TRUE);
  
  return(idata);
  
}

# Change colname
change_colname = function(data, oldp, newp)
{
  copydata  = data;
  setnames(copydata, oldp, newp)
  return(copydata);
  
}

# Change row value
change_rowvalue = function(data, fieldname, oldp, newp )
{
  copydata  = data;
  copydata[[fieldname]][copydata[[fieldname]] == oldp ] = newp;
  return(copydata);
}


## Create a plot, passes character strings. Fit a line

create_figure= function(data, attrib_list, data2)
{
  copydata = data;
  p = ggplot(na.omit(copydata), aes_string(x=attrib_list[['aes.x']], y=attrib_list[['aes.y']], 
                                           group = 1)) + 
    geom_point(shape=1) +    # Use hollow circles
    geom_smooth(method=lm) +
    geom_point(data=data2, aes_string(y=attrib_list[['aes.y']], x=attrib_list[['aes.x']]), 
               color='red', size = 3) +
    ggtitle(attrib_list[['x.title']]) +
    xlab('DAS') + ylab(attrib_list[['y.title']]) +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1))
    
  
 return(p) 
}

calibratedata = function(data)
{
  copydata = data;
  for(tdate in unique(copydata$time))
  {
    t = which(copydata$time == tdate);
    for(i in 3:(ncol(copydata)-2))
    {
      if(as.Date(tdate) <= as.Date("2015-04-01")) 
      { 
        copydata[t,i] = copydata[t,i]*0.27;
      }
      else 
      { 
        copydata[t,i] = copydata[t,i]*0.67; # 0.56
      }
    }
  }    
  return(copydata)
}


convert2numeric = function(data, col_list)
{
  copydata=data;
  for(cname in col_list)
  {
    copydata[[cname]] = as.numeric(as.character(copydata[[cname]]))
  }
  return(copydata);
}


# subset of data
average_similar_measures = function(data, featurename)
{
  copydata = data;
  copydata[[featurename]] = factor(copydata[[featurename]]);
  #print(plantid)
  # Arrange new dataset
  l <- split(copydata, copydata[[featurename]]);
  h <- do.call(rbind, lapply(1:length(l), function(x) takeMeans(l[x])));
  h = h[,c(-1,-2)];
  h = h[,c(dim(h)[2], 1:(dim(h)[2]-1))]
  h <- data.frame(idtag=unique(copydata$idtag), h);
  colnames(h)[3:(dim(h)[2])] <- colnames(scsub)[3:(dim(scsub)[2])];
  return(h);
}

# Create plot but adding a vertical line to indicate something
## TODO: MAKE It GENERIC
print_plot = function(data, ldate, ytrait, xtrait)
{
 
  #days before phenotyping
  
  dbp = seq(as.Date("2014-10-25"), as.Date("2015-04-27"), by="days")
  i = which(as.Date(unique(data$time)) == as.Date(dbp[ldate]));
  print(i)
  if(length(i) == 0) { i = which(as.Date(unique(data$time)) == as.Date(dbp[ldate+1]));}
  copydata = data;
  p = ggplot(copydata,aes(x = time, y = value, fill=intensity, group=intensity)) +
    geom_area() + geom_line(aes(ymax=value), position="stack") +
    geom_vline(xintercept = i, colour='red', size=2) +
    scale_fill_manual(values = as.character(unique(copydata$intensity))) + 
    labs(y = ytrait, x = xtrait ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          text = element_text(size=8)) + 
    ggtitle(idtagname) +
    #theme(axis.text.x = element_blank() ) + theme(legend.position="none")
    theme(axis.text.x = element_blank(),
        axis.line = element_line(colour = "black"), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) 
  return(p);
}


sum_intensities <- function(data)
{
  
  copydata = data;
  copydata = data.frame(copydata, sumint = sapply(rownames(copydata), 
                                                  function(x) sum(copydata[x,3:131])));
  return(copydata);
  
}


convert2hex <- function(data1, data2)
{
  copydata1 = data1;
  copydata2 = data2;
  hexlist = c();
  
  for(i in 1:length(copydata1))
  {
    idx = which(copydata2$i == copydata1[i])
    hexlist = append(hexlist, as.character(copydata2[idx,5]));
    
  }
  return(hexlist)
}


get_subdata = function(data, conver_table, intens)
{
  
  
  copydata = data;
  # I us want to use the same intensities # Vey crude way but it seems they're the most relevants.
  
  u = colnames(copydata) %in% intens; # Search for intens
  scc1= data.frame(idtag = copydata$idtag, 
                   time = copydata$time, copydata[,u]); # Create frame with sel intensi 
  freq = colnames(scc1)[3:dim(scc1)[2]];
  freq = as.numeric(sub("X", "", freq));
  colnames(scc1)[3:dim(scc1)[2]] = freq; # Delete X in columns
  #write.table(scc1, file='t.csv', sep=',', row.names=F);
  
  mm = as.matrix(scc1[, 3:dim(scc1)[2]]); # Deletes -4
  vector_data = unmatrix(mm,byrow=T); # Convert matrix into vector
  vector_data = as.numeric(vector_data); # Convert values to numeric
  idtag = as.character(unique(scc1$idtag)); # Unique idtags
  time = as.character(unique(scc1$time)) # unique time points
  hexlist = convert2hex(freq, conver_table); # 
  conveg = expand.grid(idtag=idtag, intensity=hexlist,time=time);
  cdata = data.frame(conveg, value=vector_data);
  #cdatacal = calibrate_data(cdata);
  return(cdata);
}

# Take row means

takeMeans <- function(x) 
{
  x1 = rapply(x, mean);
  x1['time'] = names(x); 
  return(x1)
}
# subset of data
get_subset = function(data)
{
  copydata = data;
  #print(plantid)
  # Arrange new dataset
  l <- split(copydata, copydata$timestamp);
  h <- do.call(rbind, lapply(1:length(l), function(x) takeMeans(l[x])));
  h = h[,c(-1,-2)];
  h = h[,c(dim(h)[2], 1:(dim(h)[2]-1))]
  h <- data.frame(idtag=unique(copydata$idtag), h);
  colnames(h)[3:(dim(h)[2])] <- colnames(scsub)[3:(dim(scsub)[2])];
  return(h);
}

### DB

connect2DB = function(query_list, db_list)
{
  drv <- dbDriver(db_list[['dbdriver']]);
  dbhost <- db_list[['dbip']]; 
  dbport <- db_list[['dbport']];   
  con <- dbConnect(drv, host=dbhost, port=dbport, dbname=db_list[['dbname']],
                   user=db_list[['dbuser']], password=db_list[['dbpass']]);
  rs = dbSendQuery(con, query_list);
  data <- fetch(rs, n = -1);
  dbDisconnect(con);
  return(data);
}




#Export data in weka format

create_query_temp = function(tdate)
{
  s = '%';
  st = sprintf("SELECT logdate, s3_1, s3_2,s3_3,s3_4,s3_5,s3_6,s3_7 
               FROM phenomics.glasshouse_sensors where logdate > '%s' and compartent = '5' ",tdate);
  return(gsub("\r?\n|\r", " ", st));
}

# split timestamp into date and time_ Lemnatec
process_timestamp = function(data, timename, sp, timeaction)
{
  copydata = data;
  d <- data.frame(copydata, date=sapply(copydata[[timename]],  
                                        function(x) strsplit(as.character(x), sp)[[1]][1]));
  d <- data.frame(d, time=sapply(d[[timename]],  
                                 function(x) strsplit(strsplit(as.character(x), ' ')[[1]][2],':')[[1]][1]));
  i = which(colnames(d) == 'time');
  colnames(d)[i] = timeaction;
  
  return(d);
}

# Merge two datasets
merge_data <- function(d1, d2, genoname, gname)
{

  i = which(colnames(d2) == genoname);
  colnames(d2)[i] = 'genotypename';

  print(colnames(d2));
	d <- data.frame(d1, genotypename=sapply(d1$id_tag, function(x)
	  d2[which(d2$barcode == as.character(x)), colnames(d2)[i]]));
	d <- data.frame(d, water.treatment=sapply(d$id_tag, function(x)
	  d2[which(d2$barcode == as.character(x)), 'control']));
	d <- data.frame(d, rep=sapply(d$id_tag, function(x)
	  d2[which(d2$barcode == as.character(x)), 'rep']));
	return(d);
}



## TBC

# read_data <- function(datafilename, h)
# {
# 	idata <- read.table(datafilename, sep=',', header=h);
# 	
# 	return(idata);
# 
# }
# 
# 
# splitbyidtag <- function(data)
# {
#   copydata = data;
#   t = as.character(unique(copydata$id_tag));
#   y = split(d, t)
# }
#  
# 
# 
# # This function has to be changed all the time as the metatable files is always different

# 
# 
# merge_dataheight <- function(datatable, metatable)
# {
# 	d <- data.frame(datatable, time=sapply(datatable$time_stamp, 
#                                          function(x) strsplit(as.character(x), '_')[[1]][2])); # 3 for w1
#      	d <- data.frame(d, id_tag=sapply(d$time_stamp, 
#                                         function(x) strsplit(as.character(x), '_')[[1]][1])); # 1 for w1
# 	d <- data.frame(d, genotypename=sapply(d$id_tag, 
#                                          function(x) metatable[which(metatable$id_tag == as.character(x)), 'line']));
# 	d <- data.frame(d, photoperiod=sapply(d$id_tag, 
#                                         function(x) metatable[which(metatable$id_tag == as.character(x)), 'photoperiod']));
# 	d <- data.frame(d, temperature=sapply(d$id_tag, 
#                                         function(x) metatable[which(metatable$id_tag == as.character(x)), 'temp']));
# 	d <- data.frame(d, water.treatment=sapply(d$id_tag, 
#                                             function(x) metatable[which(metatable$id_tag == as.character(x)), 'water.treatment']));
# 	d <- data.frame(d, rep=sapply(d$id_tag, 
#                                 function(x) metatable[which(metatable$id_tag == as.character(x)), 'rep']));
# 	return(d);
# }
# 
# 
# unify_data <- function(cdata)
# {
# 	nd <- c();
# 	t <- c();
# 	t <- as.character(unique(cdata$time));
# 	for(i in t) 
# 	{ 
# 		if(length(which(cdata$time == i)) < 200)
# 		{
# 			nd <- append(nd, i); 
# 		}
# 		
# 	}
# 	return(nd);
# }
# 
# search_uncompletetime <- function(cdata, nd)
# {
# 	inx <- c();
# 	for(d in nd)
# 	{
# 		i <- which(cdata$time == d);	
# 		inx <- append(inx, i);
# 
# 	}
# 	return(inx);
# 
# }
# 
# 
# merge_meta <- function(metatable)
# {
# 	d <- data.frame(metatable, id_tag=sapply(metatable$code, function(x) 
#     if(x < 10000) paste('W1_0',x, sep='') else paste('W1_',x, sep='')));
# 
# }
# 
# plot_data <- function(data)
# {
#   copydata = data;
# 	j <- which(as.Date(copydata$time) > as.Date("2013-09-01"));
#   copydata <- copydata[c(-j),];
# 	#par(mfrow=(c(3,1)));
# 	for(wt in as.numeric(unique(copydata$water.treatment)))
# 	{
# 		fname <- paste('watertreatment_', wt, '.png',sep='');
# 		png(fname, width = 880, height = 880, pointsize = 22);
# 		i <- which(copydata$water.treatment == wt);
# 		copydataw <- copydata[i,];
# 		#print(dim(cdataw)) 
# 		with(copydataw, interaction.plot(time, as.character(genotypename), 
#     log(areacal), cex.axis=0.5, legend=F, las=2, lwd=2, pch=c(1:33), lty=1:33, col=rainbow(33),
#         main=paste('Water treatment:  ', wt, sep='')));
# 		legend('topright', legend = unique(as.character(copydataw$genotypename)), bty="n", col=rainbow(32), pch=c(1:32), lty=1:33, title = 'Genotypes', cex=0.5);
# 		#if (readline(wt) == "q")    {       break;    }
# 		dev.off();
# 	}
# 
# }
# 
# 
# dummy_func <- function()
# {
# 
# 	filenames <- list.files(path='Z:/lemnatec_data/rotations/MS/0', pattern = "png");
# 	t <- as.matrix(filenames);
# 	t <- data.frame(t, name=sapply(t, function(x) strsplit(as.character(x),'.png')[[1]][1]));
# 	h <- match(as.character(t$name), as.character(datatable$filename));
# 	f <- is.na(h);
# 	y <- t[f,];
# 	write.table(y, file='missingfiles.csv', sep=',', quote=F, row.names=F);
# 
# }
# 
# 
# plotggplot <- function()
# {
# 
#   ggplot(data = cdata, aes(x = areacal, y = time, colour = line, group=line)) +  
#     stat_summary(fun.y=mean, geom="point")
# 
# }
# 
# 
# plotsdummy <- function(cdata)
# {
# 
# 	with(cdata[c(-t),], interaction.plot(time, genotypename, heightpixel, las=2, cex.axis=0.7, col=rainbow(20), lwd=2));
# 	demo1.aov <- aov(heightmm ~ genotypename * time + Error(id_tag), data = cdata)
# 	xyplot(heightmm ~ time, data = cdata, groups = id_tag,  type = "o", panel = panel.superpose)
# 	xyplot(heightmm ~ time | genotypename, data = cdata, groups = id_tag, type = "o", panel = panel.superpose);
# 	time.quad <- lme(heightmm ~ genotypename * time + time2,   random = list(id = pdDiag(~ time)), cdata) 
# 	summary(time.quad)
# 	 with(cdata,all(table(time,genotypename, water.treatment, temperature)>=1));
# 	time.quad <- lme(heightmm ~ genotypename * water.treatment,   random=~1|id_tag, cdata) 
# 
# 
# }
# 
# 
# growthcode <- function()
# {
# 
# 	gd <- read.table('growthstageguide.csv', header=T, sep=',');
# 	gd <- data.frame(gd, code=sapply(rownames(gd), function(x) paste('G', gd[x,1], gd[x,2], sep='')));
# 	return(gd);
# 
# }
# 
# 
# simulategd <- function(cdata, gd)
# {
# 	t <- as.character(unique(cdata$time));
# 	i <- which(cdata$time == t[length(t)-1]);
# 	cdata <- data.frame(cdata, gdcode=sapply(row.names(cdata), function(x) if(is.na(match(x, i))){ as.character(gd[sample(22:30,1),'code']) } else {as.character(gd[sample(31:33,1),'code'])}));
# 	return(cdata);
# 	
# 
# }
# 
# simulateftl <- function(metatable)
# {
# 
# 	v <- rep(0,200);
# 	for(type in as.character(unique(metatable$photoperiod)))
# 	{
# 		i <- which(metatable$photoperiod == type);
# 		if(type == 'Control')
# 		{
# 			s <- 30:44;
# 			v[i] <- sample(s, length(s))
# 		} 	
# 		else if(type == 'Ppd Late')
# 		{
# 			s <- 45:49;
# 			v[i] <- sample(s, length(s))
# 
# 		}
# 		else
# 		{
# 			s <- 12:30;
# 			v[i] <- sample(s, length(s))
# 
# 		}
# 
# 	}
# 
# 	 with(m, bwplot(ft ~ line | photoperiod,  xlab=list(las=2), main='Days to Flowering', scales=list(x=list(rot=90))))
# 
# }
# 
# 
# growthrate <- function(data, character_list)
# {
#  time_list <- as.numeric(as.character(unique(data$time)));
#  ecotype_list <- unique(data$ecotype);
#  time_list_adj <- c(17,22,25,28,32);
#  datacopy <- data;
# 
#  for(ecotype in ecotype_list[1:19])
#  {
#    for(time in 2:length(time_list))
#    {
#      sample_list <- as.character(unique(datacopy[intersect(which(datacopy$ecotype == ecotype), which(datacopy$time == time_list[time-1])),'id']));
#      for(sample in sample_list)
#      {
#        for(character in character_list)
#        { 
#          np <- intersect(intersect(which(datacopy$ecotype == ecotype), which(datacopy$time == time_list[time-1])), which(datacopy$id == sample));
#          nn <- intersect(intersect(which(datacopy$ecotype == ecotype), which(datacopy$time == time_list[time])), which(datacopy$id == sample));
#          if(!is.na(data[np, character]) & !is.na(data[nn, character]))
#          {
#           
#            datacopy[nn,character] <- log(data[nn, character]) - log(data[np, character]) / ( time_list_adj[time] - time_list_adj[time-1]);
#          }
#        }
#      }
#    }
#  }
#  return(datacopy);
# }
# 
# 
# wue <- function(data)
# {
#  datacopy <- data;
#  i <- union(which(datacopy$writer_label != 'W2_sv0'), which(as.Date(datacopy$time) == as.Date('2013-09-11')));
#  datacopy <- datacopy[c(-i),];
#  datacopy$time <- as.numeric(datacopy$time);
#  time_list <- unique(datacopy$time);
#  line_list <- as.character(unique(datacopy$id_tag));
# 
# 
#  for(line in line_list)
#  {
#    for(time in 2:(length(time_list)-1))
#    {
#      #sample_list <- as.character(unique(datacopy[intersect(which(datacopy$id_tag == line), which(datacopy$time == time_list[time-1])),'id']));
#      #for(sample in sample_list)
#      #{
# 		print(line)	
# 		print(time)
# 		#print(sample)
#          np <- intersect(which(datacopy$id_tag == line), which(datacopy$time == time_list[time-1]));
#          nn <- intersect(which(datacopy$id_tag == line), which(datacopy$time == time_list[time])); 
# 		print(nn);  
# 		print(np);  
# 	    	print(datacopy$weight_before[np]);    
# 		print(datacopy$weight_before[nn]);
# 		print(datacopy$water_amount[np]);
# 		print(datacopy$water_amount[nn]);           
#          datacopy[nn, 'wue'] <- (datacopy[nn[1], 'weight_after'] - datacopy[nn[1], 'weight_before']) / (datacopy[nn[1], 'water_amount'] - datacopy[np[1], 'water_amount']);
# 	   print(datacopy[nn, 'wue']);
# 	  
# 
#      #}
#    }
#  }
#  return(datacopy);
# }
# 
# # ------------- Read metadata
# 
# 
# tidy_cols = function(data, clist)
# {
#   copydata = data;
#   ttc = clist['X1']
#   i = match(ttc, colnames(copydata));
#   colnames(copydata)[i] =  clist['X2'];
#   return(copydata);
#   
# }
# 
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)

  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }

  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )

  # Rename the "mean" column
  datac <- rename(datac, c("mean" = measurevar))

  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult

  return(datac)
}
# 
# # Interaction plot water vs time
# plot_ave_area= function(data)
# {
#   copydata = data;
#   for(ename in unique(copydata$genotypename))
#   {    
#     print(ename)
#     fname = paste("Intplot_AreaPerGenotype", ename, ".png", sep="");
#     
#     i = which(copydata$genotype == ename);
#     sdata = copydata[i, ];
#     pd <- position_dodge(0.1) # move them .05 to the left and right
#     tgc <- summarySE(sdata, measurevar="area", groupvars=c('date','water.treatment','genotypename'));
#     #png(fname);
#     p = ggplot(data = tgc,
#            aes(x = date, y = area, colour = water.treatment, 
#                group=interaction(water.treatment, genotypename))) +
#       ggtitle(ename) +
#       geom_path(alpha = 0.5) +
#       geom_errorbar(aes(ymin=area-se, ymax=area+se), width=.1, position=pd) +
#       stat_summary(fun.y=mean, geom="point")+
#       stat_summary(fun.y=mean, geom="line") +
#       theme(axis.text.x=element_text(angle=90),
#             axis.text.x = element_text(size=20));
#     print(p);
#     if(readline(fname) == 'q') { break();}
#     #dev.off();
#   }    
# }
# 
# 
# 
# plot_int_data = function(data, genotypename_list)
# {
#   copydata = data;
#   for(ename in genotypename_list)
#   {    
#     print(ename)
#     fname = paste("Intplot_AreaPerGenotype", ename, ".png", sep="");
#     
#     i = which(copydata$genotype == ename);
#     sdata = copydata[i, ];
#     with(sdata, interaction.plot(date, water.treatment, area, col=rainbow(2),
#                      las=2, cex.axis=0.6, lwd=2, main=ename));
#     if(readline(fname) == 'q') { break();}
#   }
# }
# 
# 
# plot_area= function(data)
# {
#   copydata = data;
#   for(ename in unique(copydata$genotypename))
#   {    
#     fname = paste("Intplot_AreaPerPlant", ename, ".png", sep="");
#     
#     i = which(copydata$genotype == ename);
#     sdata = copydata[i, ];
#     pd <- position_dodge(0.1) # move them .05 to the left and right
#     tgc <- summarySE(sdata, measurevar="area", groupvars=c('water.treatment', 'genotypename', 'date', 'id_tag'));
#     png(fname);
#     p = ggplot(data = sdata,
#                aes(x = date, y = area, colour = water.treatment, shape = id_tag,
#                    group=interaction(water.treatment, id_tag))) +
#       ggtitle(ename) +
#       geom_path(alpha = 0.5) +
#       geom_errorbar(aes(ymin=area-se, ymax=area+se), width=.1, position=pd) +
#       stat_summary(fun.y=mean, geom="point")+
#       stat_summary(fun.y=mean, geom="line") +
#       theme(axis.text.x=element_text(angle=90),
#             axis.text.x = element_text(size=20));
#     print(p);
#     #if(readline(fname) == 'q') { break();}
#     dev.off();
#   }    
# }
# 
# # Plot water use efficiency
plot_wue = function(data)
{
  copydata = data;
  for(ename in unique(copydata$genotypename))
  {
    fname = paste("Intplot_AreaPerGenotype_WUE", ename, ".png", sep="");

    i = which(copydata$genotype == ename);
    sdata = copydata[i, ];
    pd <- position_dodge(0.1) # move them .05 to the left and right
    tgc <- summarySE(sdata, measurevar="wue",
                     groupvars=c('water.treatment', 'genotypename', 'date', 'temp'));
    png(fname);
    mn = min(sdata$temp);
    mx = max(sdata$temp);
    p = ggplot(data = tgc,
        aes(x = date, y = wue, colour=temp, shape=water.treatment,
        group=interaction(water.treatment, genotypename))) +
        ggtitle(ename) +
        geom_path(alpha = 0.5) +
        geom_errorbar(aes(ymin=wue-se, ymax=wue+se), width=.1, position=pd) +
        stat_summary(fun.y=mean, geom="point", size = 3)+
        stat_summary(fun.y=mean, geom="line") +
        theme(axis.text.x=element_text(angle=90),
        axis.text.x = element_text(size=20)) +
        scale_colour_gradient(limits=c(mn, mx), low="yellow", high="red");;
        print(p);
    #if(readline(fname) == 'q') { break();}
    dev.off();
  }
}
# 
# # Plot temperature
# plot_temp = function(data)
# {
#   copydata = data;
#   for(ename in unique(copydata$genotypename))
#   {    
#     fname = paste("Intplot_ArePerGenotype_temp", ename, ".png", sep="");
#     
#     i = which(copydata$genotype == ename);
#     sdata = copydata[i, ];
#     pd <- position_dodge(0.1) # move them .05 to the left and right
#     tgc <- summarySE(sdata, measurevar="area", groupvars=c('water.treatment', 
#                      'genotypename', 'date', 'temp'));
#     png(fname);
#     mn = min(sdata$temp);
#     mx = max(sdata$temp);
#     p = ggplot(data = tgc, aes(x=date, y=area, colour=temp, shape=water.treatment, group=
#       interaction(water.treatment, temp))) +
#       ggtitle(ename) + geom_path(alpha = 0.5) + 
#       geom_errorbar(aes(ymin=area-se, ymax=area+se), width=.1, position=pd) +
#       stat_summary(fun.y=mean, geom="point", size=4) +
#       stat_summary(fun.y=mean, geom="line") +
#       theme(axis.text.x=element_text(angle=90),
#             axis.text.x = element_text(size=20)) +
#       scale_colour_gradient(limits=c(mn, mx), low="yellow", high="red");
#     print(p);
#     #if(readline(fname) == 'q') { break();}
#     dev.off();
#   }      
# }
# 
# 
# #average data
# avgdata <- function(data, fields, variables)
# {
#   print(fields)
#   print(variables)
#   copydata <- data;
#   copydata = aggregate(cbind(fields) ~ variables, data = copydata, FUN= "mean" )
#   return(copydata) 
# }
# 
# 
# # Split temp data by date by hour
# splittime = function(data)
# {
#   
#   copydata = data;
#   d = strsplit(copydata, ' ')[[1]][1];
#   h = strsplit(strsplit(copydata, ' ')[[1]][2], ':')[[1]][1];
#   return(paste(d,h));
# }
# 

PlotCorrelation = function(data)
{
  dta = data;
  # Plot correlations
  dta.r <- cor(dta, use = 'complete.obs') # get correlations
  dta.col <- dmat.color(dta.r) # get colors
  dta.o <- order.single(dta.r);
  cpairs(dta, dta.o, panel.colors=dta.col, gap=.5, cex.labels=1.2); 
  
}

pcacomplete <- function(data)
{
  copydata <- data;
  res <- PCA(copydata, quali.sup=c(1:5), graph=F, scale.unit=TRUE);
  #png(fname, width = 2000, height = 880, pointsize = 22);
  plot(res, choix = "var", cex=0.9, label='var',  col.var = rainbow(10), lwd=2, 
       font.lab=2, font=2, ylab='hh');
  #dev.off();
  plotloadings(res$var$contrib, time);
  #lo <- sweep(res$var$coord,2,sqrt(res$eig[1:ncol(res$var$coord),1]),FUN="/");
  #plotloadings(lo, time);
  res <- PCA(copydata[6:25], graph=F);
  u = dimdesc(res, axes=c(1,2));
  return(res);
  
}

