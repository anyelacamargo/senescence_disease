library('ggplot2'); # install package like that install.packages('ggplot2')
#library('RODBC');
#library('RPostgreSQL');
#library('RMySQL');



# Function to read data from external files
read_data <- function(datafilename, h)
{
  idata <- read.table(datafilename, sep=',', header=h);
  
  return(idata);
  
}

