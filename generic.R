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


# Create a plot, passes character strings. Fit a line
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