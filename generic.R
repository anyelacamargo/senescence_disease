library('ggplot2'); 
library(data.table);#
#library('RODBC');
#library('RPostgreSQL');
#library('RMySQL');



# Function to read data from external files
read_data = function(datafilename, h)
{
  idata <- read.table(datafilename, sep=',', header=h);
  
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

create_figure= function(data, attrib_list)
{
  copydata = data;
  p = ggplot(na.omit(copydata), aes_string(x=attrib_list[['aes.x']], y=attrib_list[['aes.y']], group = 1)) +        
    geom_point(shape=1) +    # Use hollow circles
    geom_smooth(method=lm) +
    xlab(attrib_list[['x.title']]) + ylab(attrib_list[['y.title']]) +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1)) 
 return(p) 
}


# Create plot but adding a vertical line to indicate something
## TODO: MAKE It GENERIC
print_plot = function(data, ldate)
{
  print(ldate)
  #days before phenotyping
  dbp = 100;
  copydata = data;
  p = ggplot(copydata,aes(x = time, y = areacal, fill=intensity, group=intensity)) +
    geom_area() + geom_line(aes(ymax=areacal), position="stack") +
    geom_vline(xintercept = (ldate - dbp), colour='red', size=2) +
    scale_fill_manual(values = as.character(unique(copydata$intensity))) + 
    #theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle(idtagname) +
    theme(axis.text.x = element_blank() ) + theme(legend.position="none") + 
    theme(axis.line = element_line(colour = "black"), 
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

calibration_eq <- function(d)
{
  
  if(as.Date(d$time) <= as.Date("2015-04-01")) 
  { 
    o = d$value*0.27;
  }
  else 
  { 
    o = d$value*0.67;#1.5
  }
  
  return(o);
}

#Calibrate dataset
calibrate_data <- function(data)
{
  copydata <- data;
  copydata <- data.frame(copydata, areacal=sapply(row.names(copydata), 
                                                  function(x) calibration_eq(copydata[x,])));
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


#scsub = scsub[order(as.Date(cdata$time, format="%y%m%d")),]
## 
# vc = c();
# for(i in 3:dim(sc)[2])
# {
#   vc = append(vc,mean(sc[,i]))
# }
# v = which(vc <= 10);
# v = v+2;
# scc1 = sc[,c(-3, -123:-127, -v)];

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
  cdatacal = calibrate_data(cdata);
  return(cdatacal);
}

# Take row means

takeMeans <- function(x) 
{
  x1 = rapply(x, mean);
  x1['time'] = names(x); 
  return(x1)
}

# Plot the thing;


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
