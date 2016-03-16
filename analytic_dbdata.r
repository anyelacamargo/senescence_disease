library('ggplot2'); # install package like that install.packages('ggplot2')
library('RODBC');
library('RPostgreSQL');
library('RMySQL');


# Set working directory

datatamerge_files <- function()
{
	filenames <- list.files(pattern = "txt");
	#c <- do.call("rbind", lapply(filenames, read.csv, header = F));
	c <- read.table('output.csv', sep=',', header=F);
}


connect2DB = function(y)
{
  drv <- dbDriver("PostgreSQL");
  dbhost <- '144.124.105.250'; 
  dbport <- '5432';   
  con <- dbConnect(drv, host=dbhost, port=dbport, dbname='production',
                   user='postgres', password='LemnaTec');
  rs = dbSendQuery(con, y);
  data <- fetch(rs, n = -1);
  dbDisconnect(con);
  return(data);
}

connect2DBmysql = function(y)
{
  drv <- dbDriver("MySQL");
  dbhost <- 'venom.ibers.aber.ac.uk'; 
  dbport <- '3306';   
  con <- dbConnect(drv, host=dbhost, dbname='phenomics',
                   user='dummy', password='google2ooo');
  rs = dbSendQuery(con, y);
  data <- fetch(rs, n = -1);
  dbDisconnect(con);
  return(data);
}

read_data <- function(datafilename, h)
{
	idata <- read.table(datafilename, sep=',', header=h);
	
	return(idata);

}


splitbyidtag <- function(data)
{
  copydata = data;
  t = as.character(unique(copydata$id_tag));
  y = split(d, t)
}
 
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

# This function has to be changed all the time as the metatable files is always different
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


merge_dataheight <- function(datatable, metatable)
{
	d <- data.frame(datatable, time=sapply(datatable$time_stamp, 
                                         function(x) strsplit(as.character(x), '_')[[1]][2])); # 3 for w1
     	d <- data.frame(d, id_tag=sapply(d$time_stamp, 
                                        function(x) strsplit(as.character(x), '_')[[1]][1])); # 1 for w1
	d <- data.frame(d, genotypename=sapply(d$id_tag, 
                                         function(x) metatable[which(metatable$id_tag == as.character(x)), 'line']));
	d <- data.frame(d, photoperiod=sapply(d$id_tag, 
                                        function(x) metatable[which(metatable$id_tag == as.character(x)), 'photoperiod']));
	d <- data.frame(d, temperature=sapply(d$id_tag, 
                                        function(x) metatable[which(metatable$id_tag == as.character(x)), 'temp']));
	d <- data.frame(d, water.treatment=sapply(d$id_tag, 
                                            function(x) metatable[which(metatable$id_tag == as.character(x)), 'water.treatment']));
	d <- data.frame(d, rep=sapply(d$id_tag, 
                                function(x) metatable[which(metatable$id_tag == as.character(x)), 'rep']));
	return(d);
}


unify_data <- function(cdata)
{
	nd <- c();
	t <- c();
	t <- as.character(unique(cdata$time));
	for(i in t) 
	{ 
		if(length(which(cdata$time == i)) < 200)
		{
			nd <- append(nd, i); 
		}
		
	}
	return(nd);
}

search_uncompletetime <- function(cdata, nd)
{
	inx <- c();
	for(d in nd)
	{
		i <- which(cdata$time == d);	
		inx <- append(inx, i);

	}
	return(inx);

}


merge_meta <- function(metatable)
{
	d <- data.frame(metatable, id_tag=sapply(metatable$code, function(x) 
    if(x < 10000) paste('W1_0',x, sep='') else paste('W1_',x, sep='')));

}

plot_data <- function(data)
{
  copydata = data;
	j <- which(as.Date(copydata$time) > as.Date("2013-09-01"));
  copydata <- copydata[c(-j),];
	#par(mfrow=(c(3,1)));
	for(wt in as.numeric(unique(copydata$water.treatment)))
	{
		fname <- paste('watertreatment_', wt, '.png',sep='');
		png(fname, width = 880, height = 880, pointsize = 22);
		i <- which(copydata$water.treatment == wt);
		copydataw <- copydata[i,];
		#print(dim(cdataw)) 
		with(copydataw, interaction.plot(time, as.character(genotypename), 
    log(areacal), cex.axis=0.5, legend=F, las=2, lwd=2, pch=c(1:33), lty=1:33, col=rainbow(33),
        main=paste('Water treatment:  ', wt, sep='')));
		legend('topright', legend = unique(as.character(copydataw$genotypename)), bty="n", col=rainbow(32), pch=c(1:32), lty=1:33, title = 'Genotypes', cex=0.5);
		#if (readline(wt) == "q")    {       break;    }
		dev.off();
	}

}


dummy_func <- function()
{

	filenames <- list.files(path='Z:/lemnatec_data/rotations/MS/0', pattern = "png");
	t <- as.matrix(filenames);
	t <- data.frame(t, name=sapply(t, function(x) strsplit(as.character(x),'.png')[[1]][1]));
	h <- match(as.character(t$name), as.character(datatable$filename));
	f <- is.na(h);
	y <- t[f,];
	write.table(y, file='missingfiles.csv', sep=',', quote=F, row.names=F);

}


plotggplot <- function()
{

  ggplot(data = cdata, aes(x = areacal, y = time, colour = line, group=line)) +  
    stat_summary(fun.y=mean, geom="point")

}


plotsdummy <- function(cdata)
{

	with(cdata[c(-t),], interaction.plot(time, genotypename, heightpixel, las=2, cex.axis=0.7, col=rainbow(20), lwd=2));
	demo1.aov <- aov(heightmm ~ genotypename * time + Error(id_tag), data = cdata)
	xyplot(heightmm ~ time, data = cdata, groups = id_tag,  type = "o", panel = panel.superpose)
	xyplot(heightmm ~ time | genotypename, data = cdata, groups = id_tag, type = "o", panel = panel.superpose);
	time.quad <- lme(heightmm ~ genotypename * time + time2,   random = list(id = pdDiag(~ time)), cdata) 
	summary(time.quad)
	 with(cdata,all(table(time,genotypename, water.treatment, temperature)>=1));
	time.quad <- lme(heightmm ~ genotypename * water.treatment,   random=~1|id_tag, cdata) 


}


growthcode <- function()
{

	gd <- read.table('growthstageguide.csv', header=T, sep=',');
	gd <- data.frame(gd, code=sapply(rownames(gd), function(x) paste('G', gd[x,1], gd[x,2], sep='')));
	return(gd);

}


simulategd <- function(cdata, gd)
{
	t <- as.character(unique(cdata$time));
	i <- which(cdata$time == t[length(t)-1]);
	cdata <- data.frame(cdata, gdcode=sapply(row.names(cdata), function(x) if(is.na(match(x, i))){ as.character(gd[sample(22:30,1),'code']) } else {as.character(gd[sample(31:33,1),'code'])}));
	return(cdata);
	

}

simulateftl <- function(metatable)
{

	v <- rep(0,200);
	for(type in as.character(unique(metatable$photoperiod)))
	{
		i <- which(metatable$photoperiod == type);
		if(type == 'Control')
		{
			s <- 30:44;
			v[i] <- sample(s, length(s))
		} 	
		else if(type == 'Ppd Late')
		{
			s <- 45:49;
			v[i] <- sample(s, length(s))

		}
		else
		{
			s <- 12:30;
			v[i] <- sample(s, length(s))

		}

	}

	 with(m, bwplot(ft ~ line | photoperiod,  xlab=list(las=2), main='Days to Flowering', scales=list(x=list(rot=90))))

}


growthrate <- function(data, character_list)
{
 time_list <- as.numeric(as.character(unique(data$time)));
 ecotype_list <- unique(data$ecotype);
 time_list_adj <- c(17,22,25,28,32);
 datacopy <- data;

 for(ecotype in ecotype_list[1:19])
 {
   for(time in 2:length(time_list))
   {
     sample_list <- as.character(unique(datacopy[intersect(which(datacopy$ecotype == ecotype), which(datacopy$time == time_list[time-1])),'id']));
     for(sample in sample_list)
     {
       for(character in character_list)
       { 
         np <- intersect(intersect(which(datacopy$ecotype == ecotype), which(datacopy$time == time_list[time-1])), which(datacopy$id == sample));
         nn <- intersect(intersect(which(datacopy$ecotype == ecotype), which(datacopy$time == time_list[time])), which(datacopy$id == sample));
         if(!is.na(data[np, character]) & !is.na(data[nn, character]))
         {
          
           datacopy[nn,character] <- log(data[nn, character]) - log(data[np, character]) / ( time_list_adj[time] - time_list_adj[time-1]);
         }
       }
     }
   }
 }
 return(datacopy);
}


wue <- function(data)
{
 datacopy <- data;
 i <- union(which(datacopy$writer_label != 'W2_sv0'), which(as.Date(datacopy$time) == as.Date('2013-09-11')));
 datacopy <- datacopy[c(-i),];
 datacopy$time <- as.numeric(datacopy$time);
 time_list <- unique(datacopy$time);
 line_list <- as.character(unique(datacopy$id_tag));


 for(line in line_list)
 {
   for(time in 2:(length(time_list)-1))
   {
     #sample_list <- as.character(unique(datacopy[intersect(which(datacopy$id_tag == line), which(datacopy$time == time_list[time-1])),'id']));
     #for(sample in sample_list)
     #{
		print(line)	
		print(time)
		#print(sample)
         np <- intersect(which(datacopy$id_tag == line), which(datacopy$time == time_list[time-1]));
         nn <- intersect(which(datacopy$id_tag == line), which(datacopy$time == time_list[time])); 
		print(nn);  
		print(np);  
	    	print(datacopy$weight_before[np]);    
		print(datacopy$weight_before[nn]);
		print(datacopy$water_amount[np]);
		print(datacopy$water_amount[nn]);           
         datacopy[nn, 'wue'] <- (datacopy[nn[1], 'weight_after'] - datacopy[nn[1], 'weight_before']) / (datacopy[nn[1], 'water_amount'] - datacopy[np[1], 'water_amount']);
	   print(datacopy[nn, 'wue']);
	  

     #}
   }
 }
 return(datacopy);
}

# ------------- Read metadata


tidy_cols = function(data, clist)
{
  copydata = data;
  ttc = clist['X1']
  i = match(ttc, colnames(copydata));
  colnames(copydata)[i] =  clist['X2'];
  return(copydata);
  
}

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

# Interaction plot water vs time
plot_ave_area= function(data)
{
  copydata = data;
  for(ename in unique(copydata$genotypename))
  {    
    print(ename)
    fname = paste("Intplot_AreaPerGenotype", ename, ".png", sep="");
    
    i = which(copydata$genotype == ename);
    sdata = copydata[i, ];
    pd <- position_dodge(0.1) # move them .05 to the left and right
    tgc <- summarySE(sdata, measurevar="area", groupvars=c('date','water.treatment','genotypename'));
    #png(fname);
    p = ggplot(data = tgc,
           aes(x = date, y = area, colour = water.treatment, 
               group=interaction(water.treatment, genotypename))) +
      ggtitle(ename) +
      geom_path(alpha = 0.5) +
      geom_errorbar(aes(ymin=area-se, ymax=area+se), width=.1, position=pd) +
      stat_summary(fun.y=mean, geom="point")+
      stat_summary(fun.y=mean, geom="line") +
      theme(axis.text.x=element_text(angle=90),
            axis.text.x = element_text(size=20));
    print(p);
    if(readline(fname) == 'q') { break();}
    #dev.off();
  }    
}



plot_int_data = function(data, genotypename_list)
{
  copydata = data;
  for(ename in genotypename_list)
  {    
    print(ename)
    fname = paste("Intplot_AreaPerGenotype", ename, ".png", sep="");
    
    i = which(copydata$genotype == ename);
    sdata = copydata[i, ];
    with(sdata, interaction.plot(date, water.treatment, area, col=rainbow(2),
                     las=2, cex.axis=0.6, lwd=2, main=ename));
    if(readline(fname) == 'q') { break();}
  }
}


plot_area= function(data)
{
  copydata = data;
  for(ename in unique(copydata$genotypename))
  {    
    fname = paste("Intplot_AreaPerPlant", ename, ".png", sep="");
    
    i = which(copydata$genotype == ename);
    sdata = copydata[i, ];
    pd <- position_dodge(0.1) # move them .05 to the left and right
    tgc <- summarySE(sdata, measurevar="area", groupvars=c('water.treatment', 'genotypename', 'date', 'id_tag'));
    png(fname);
    p = ggplot(data = sdata,
               aes(x = date, y = area, colour = water.treatment, shape = id_tag,
                   group=interaction(water.treatment, id_tag))) +
      ggtitle(ename) +
      geom_path(alpha = 0.5) +
      geom_errorbar(aes(ymin=area-se, ymax=area+se), width=.1, position=pd) +
      stat_summary(fun.y=mean, geom="point")+
      stat_summary(fun.y=mean, geom="line") +
      theme(axis.text.x=element_text(angle=90),
            axis.text.x = element_text(size=20));
    print(p);
    #if(readline(fname) == 'q') { break();}
    dev.off();
  }    
}

# Plot water use efficiency
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

# Plot temperature
plot_temp = function(data)
{
  copydata = data;
  for(ename in unique(copydata$genotypename))
  {    
    fname = paste("Intplot_ArePerGenotype_temp", ename, ".png", sep="");
    
    i = which(copydata$genotype == ename);
    sdata = copydata[i, ];
    pd <- position_dodge(0.1) # move them .05 to the left and right
    tgc <- summarySE(sdata, measurevar="area", groupvars=c('water.treatment', 
                     'genotypename', 'date', 'temp'));
    png(fname);
    mn = min(sdata$temp);
    mx = max(sdata$temp);
    p = ggplot(data = tgc, aes(x=date, y=area, colour=temp, shape=water.treatment, group=
      interaction(water.treatment, temp))) +
      ggtitle(ename) + geom_path(alpha = 0.5) + 
      geom_errorbar(aes(ymin=area-se, ymax=area+se), width=.1, position=pd) +
      stat_summary(fun.y=mean, geom="point", size=4) +
      stat_summary(fun.y=mean, geom="line") +
      theme(axis.text.x=element_text(angle=90),
            axis.text.x = element_text(size=20)) +
      scale_colour_gradient(limits=c(mn, mx), low="yellow", high="red");
    print(p);
    #if(readline(fname) == 'q') { break();}
    dev.off();
  }      
}


#average data
avgdata <- function(data, fields, variables)
{
  print(fields)
  print(variables)
  copydata <- data;
  copydata = aggregate(cbind(fields) ~ variables, data = copydata, FUN= "mean" )
  return(copydata) 
}


# Split temp data by date by hour
splittime = function(data)
{
  
  copydata = data;
  d = strsplit(copydata, ' ')[[1]][1];
  h = strsplit(strsplit(copydata, ' ')[[1]][2], ':')[[1]][1];
  return(paste(d,h));
}


create_query_image = function(expname, angule)
{
  si = '%'

  st = sprintf("SELECT r2_obj.analysis_id, snapshot.id_tag, snapshot.time_stamp, r2_obj.area, 
         r2_obj.writer_label
                  FROM 
                  public.r2_obj, public.snapshot, public.analysis
                  WHERE 
                  analysis.snapshot_id = snapshot.id AND
                  analysis.id = r2_obj.analysis_id AND
                  snapshot.id_tag LIKE  '%s_%s'and writer_label like '%s_sv%s%s'
                  ORDER BY
                  snapshot.time_stamp ASC;", expname,si,expname,angule, si);
  
  return(gsub("\r?\n|\r", " ", st));
  
}

create_query_water = function(expname)
{
  si = '%';
  st = sprintf("SELECT snapshot.id_tag, snapshot.time_stamp, weight_before, weight_after,
snapshot.water_amount
  FROM public.snapshot WHERE snapshot.id_tag LIKE '%s_%s' AND snapshot.water_amount <> '-1'
  ORDER BY snapshot.time_stamp ASC;",expname,si);
  
  return(gsub("\r?\n|\r", " ", st));
  
}

#Export data in weka format

create_query_temp = function(tdate)
{
  s = '%';
  st = sprintf("SELECT logdate, s3_1, s3_2,s3_3,s3_4,s3_5,s3_6,s3_7 
               FROM phenomics.glasshouse_sensors where logdate > '%s' and compartent = '5' ",tdate);
  return(gsub("\r?\n|\r", " ", st));
}

# ------- Call code ---------------------


# Read metadata
metafilename <- 'metadata.csv';
metatable <- read_data(metafilename, TRUE);# Get metadata
expname = strsplit(as.vector(metatable$barcode[1]), '-')[[1]][1] # Get exp code


#Query sql for image data

#i_000 = create_query_image(expname, '0');
#i_045 = create_query_image(expname, '45');
i_090 = create_query_image(expname, '90');

# Query for water data
w = create_query_water(expname);
# Query for env data
t = create_query_temp('2015-01-17'); # search for temp sata from that date onwards


# Connect to db, query and get results into dataframe
#image_000_raw <- connect2DB(i_000);
#image_045_raw <- connect2DB(i_045);
image_090_raw <- connect2DB(i_090);
water_data_raw = connect2DB(w);
temp_data_raw = connect2DBmysql(t);


# Merge sideviews
#data_all = merge(image_000_raw, image_090_raw, by.x='analysis_id', by.y='analysis_id');
#delname = c("id_tag.x","time_stamp.x", "writer_label.x", "writer_label.y"); 
data_all = image_090_raw;
delname = c("writer_label"); 
i = match(delname, colnames(data_all));
data_all = data_all[, c(-i)];


# Tidy time_spam
#i = match(c("area.x","id_tag.y", "time_stamp.y","area.y"), colnames(data_all));
i = match(c("area"), colnames(data_all));

# Rename cols
colnames(data_all)[i] = c('area_090');

# Change timestap to leave data only
data_all = process_timestamp(data_all, 'time_stamp',' ', 'timeimaged'); 
# Delete col
data_all = data_all[, c(-(which(colnames(data_all) == 'time_stamp')))];

# Tidy water data
water_data = process_timestamp(water_data_raw, 'time_stamp',' ', 'timewatered');
water_data = water_data[, c(-(which(colnames(water_data) == 'time_stamp')))];
water_data$water_amount[water_data$water_amount == 0] = 0.01;

# Tridy temp data
# Average by date by hour
temp_data_raw1 = data.frame(temp_data_raw, datehour = 
                sapply(temp_data_raw$logdate, function(x) splittime(x)));

temp_data = process_timestamp(temp_data_raw, 'logdate', ' ', 'timetemp');
# Merge area data with water data and temp
u = merge(data_all, water_data, by.x=c('id_tag', 'date'), by.y=c('id_tag', 'date'));

cdata <- merge_data(u, metatable, 'genotype'); # Change column name
colnames(cdata)[4] = 'area';
areanames = c('area');
# Average data per factor
cdata = aggregate(cbind(area, water_amount, 
                        weight_before, weight_after) ~ 
        date + timewatered + timeimaged + id_tag + genotypename + water.treatment +  rep, 
        data = cdata, FUN= "mean");

cdata$area = cdata$area * 0.02;  # Calibrate pixels to mm

cdatawue = cdata;
cdatawue$water_amount[cdatawue$water_amount == 0.01] = NA;
cdatawue = cdatawue[!is.na(cdatawue$water_amount),];
cdatawue = data.frame(cdatawue, wue = sapply(row.names(cdatawue), 
          function(x) cdatawue[x,'area']/cdatawue[x,'water_amount'])); # area/water_amount
m = which(as.Date(cdatawue$date) <= as.Date('2015-05-09'));
cdatawue = cdatawue[-m,]; # delete rows whose time is <= 2015-05-11
break();
# Plots

plot_int_data(cdata, unique(cdata$genotypename));
# Write final data to file
write.table(cdata, file='data_pheno.csv', sep=',', quote=F, row.names=FALSE);


