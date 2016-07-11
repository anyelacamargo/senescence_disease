rm(list=ls()); # Delete files
cat('\014') 
source('generic.r');
library(gclus);
library(ggplot2)
library(FactoMineR);
library(VIM);
library('mice');
library('psych');


calculatePvalue = function(data)
{
  l = list();
  pop2Mean = data;
  l$ccpop2Mean = cor(pop2Mean)
  numcol = ncol(pop2Mean);
  l$pvaluep2Mean = l$ccpop2Mean;
  
  for(i in 1:numcol)
  {
    for(j in 1:numcol)
    {
      l$pvaluep2Mean[i,j] <- cor.test(pop2Mean[,i],pop2Mean[,2],method= "spearman")$p.value;
    }
  }
  
  l$ccpop2Mean[l$pvaluep2Mean > 0.05] = NA;
  return(l);
}






writeHeatMatp <- function(data, pval, ftitle, colna, labnames, ks)
{
  copydata <- data
  heatmap.2(as.matrix(copydata), 
            key=T, keysize=ks,key.ylab=NA,density.info="none",
            col=bluered(256),na.color=colna,sepcolor='black',
            labRow=labnames,labCol=labnames,
            cellnote = round(as.matrix(pval),2), notecex= 0.6, cexRow= 0.9, cexCol = 0.9,
            trace="none",  notecol="black", main=ftitle, 
            Rowv=FALSE, 
            Colv=FALSE,
            #dendrogram=c("col"), 
            symbreaks=T, scale="none" 
  );
  
}

imputeData = function(data)
{
  copydata = data;
  tempdata <- mice(copydata,m=5,maxit=50,meth='pmm',seed=500)
  #summary(tempData)
  #xyplot(tempData,Ozone ~ Wind+Temp+Solar.R,pch=18,cex=1)
  print(densityplot(tempdata))
  #print(stripplot(tempdata, pch = 20, cex = 1.2))
  return(tempdata);
}

replaceImputedValue = function(data, trait_list)
{
  print(trait_list)
  d = list();
  d$copydata = data;
  d$impdata = imputeData(d$copydata[,trait_list]);
  
  for(trait in trait_list)
  {
    d$copydata[row.names(d$impdata$imp[[trait]]),trait] = d$impdata$imp[[trait]]$`1`;
  }
  
  return(d);
}

# Make sure you use ggplot for 2.14


fitmodel = function(data)
{
  copydata = data;
  f = aov(FLS ~ GS39 + GS55 + GS65 +d1 + d2 + d3 +PW+SH+TIL+TN+FEL+SEL+TEL+FFLL+FEW+OEW+SM, data = copydata);
  print(summary(f));
}

convert2Numeric = function(data, trait_list)
{
  copydata = data;
  # convert data to numeric
  for(trait in trait_list)
  {
    copydata[,trait] = as.numeric(as.character(copydata[,trait]));
  }
  return(copydata);
}

plotforDisease = function(data)
{
  copydata = data;
  plotname_list = list('main' = '', 'x.title' = 'DAS',
                       'y.title' = 'Disease score', 'aes.x' = 'FLS', 
                       'aes.y' = 'SM');
  
  # Plot per trait
  for(cname in colnames(copydata)[12:29])
  {
    
    plotname_list[['aes.x']] = cname;
    plotname_list[['x.title']] = cname; 
    a = aggregate(copydata[[cname]] ~ SM, data=copydata, FUN='mean');
    colnames(a)[2] = cname;
    p = create_figure(copydata, plotname_list, a);
    tname = paste(strsplit(cname, '[.]')[[1]], collapse = '');
    png(paste(tname, '.png', sep=''));
    print(p);
    dev.off();
    #if(readline(cname) == 'q') {break();}
  }
}

plotDiagrams = function(v1, v2, dataset, trait, genotypes, xl, yl, traitname)
{
  cd = c("blue", "red","gray", 'green4', 'black', 'yellow', 'purple', 'brown', 'deepskyblue1', 
         'coral1', 'green', 'magenta1',
         rep('pink', 208));
  p = qplot(v1, v2, data=dataset, 
            #geom=c("point", "smooth"), 
            #method="lm", formula=y~x, 
            colour= genotypes, size = trait,
            xlab=xl, ylab=yl) + theme_bw() + geom_point(position=position_dodge(width=0.5)) + 
    theme(panel.grid.minor = element_blank()) + scale_colour_manual(values= cd) +
    labs(size=traitname);
  return(p);
}


plotHist = function(data, traitlist,coor)
{
  copydata = data;
  par(mfrow=c(coor))
  for(trait in traitlist)
  {
    hist(dtaM[[trait]], col='red', xlab=trait, ylab='frequency', main='');
  }
}

replaceNumber = function(data, trait_list)
{
  copydata = data;
  for(cname in trait_list)
  {
    
    i = which(copydata[,cname] == 0);
    copydata[i,cname] = NA;
  }
  return(copydata);
}

averageValues = function(data)
{
  copydata = data;
  
  copydata = aggregate(cbind(TC, GS39, GS55, GS65, FLS, PW, SH, TIL, TN, FEL, SEL, TEL, FFLL, 
                      FEW, OEW, SM) ~ genotype + type, data = copydata, FUN= "mean", na.action = na.pass);
  return(copydata);
}

createRatios = function(data)
{
  copydata = data;
  copydata = data.frame(copydata,  d1 = sapply(rownames(copydata), 
                                 function(x) copydata[x,'GS55'] - copydata[x,'GS39']));
  copydata = data.frame(copydata,  d2 = sapply(rownames(copydata), 
                                 function(x) copydata[x,'GS65'] - copydata[x,'GS55']));
  copydata = data.frame(copydata,  d3 = sapply(rownames(copydata), 
                                 function(x) copydata[x,'FLS'] - copydata[x,'GS55']));
  return(copydata);
}


plotPCA = function(data)
{
  copydata = data;
  rownames(copydata) = c(as.character(dtaM$genotype[1:12]), as.character(13:nrow(copydata)));
 
  i = copydata[,c(2, 4:ncol(copydata))];
  res <- PCA(copydata[,c(2, 4:ncol(copydata))], graph=T, quali.sup=c(1), scale.unit=TRUE);
  pdf('PCARes.pdf');
  plot.PCA(res, axes=c(1, 2), choix="ind")
  plot(res, choix = "var", cex=0.9, label='var',  col.var = rainbow(10), 
       lwd=2, font.lab=2, font=2, ylab='hh');
  dev.off();
}

preproccesdata = function(data)
{
  copydata = data;
  
  # Noise
  rname_list = data.frame(rbind(c(colname ='TC', oldname = '>', newname = 11),
                                c(colname ='TC', oldname = '<', newname = 10),
                                c(colname = 'FLS', oldname = '*** base rot ', newname = NA),
                                c(colname = 'FLS', oldname = '171?', newname = NA),
                                c(colname = 'FLS', oldname = '', newname = NA),
                                c(colname = 'PW', oldname = '6*** contaminant', newname = NA),
                                c(colname = 'GS55', oldname = '?', newname = NA),
                                c(colname = 'GS39', oldname = '?', newname = NA),
                                c(colname = 'FFLL', oldname = '?', newname = NA),
                                c(colname = 'FFLL', oldname = '', newname = NA)));
  
  # Change row values
  for(i in 1:nrow(rname_list))
  {
    copydata = change_rowvalue(copydata, as.vector(rname_list$colname[i]), 
                                       as.vector(rname_list$oldname[i]),
                                       as.vector(rname_list$newname[i])); 
  }
  return(copydata);
}

changeCols = function(data)
{
  copydata = data;
  cname_list = list('X' = 'genotype', 'row' = 'location', 'X1..earLength' = 'FEL', 
                    'X2nd.ear.L' = 'SEL', 'X3rd.ear.L' = 'TEL', 
                    'X1.fg.lf.L' = 'FFLL','X1.ear.wt' = 'FEW',
                    'other.ear.wt' = 'OEW', 'plant.wt' = 'PW', 
                    'stem.ht' = 'SH', 'top.internodeL' = 'TIL',
                    'gs39' = 'GS39', 'gs55' = 'GS55', 'flag.leaf.senescence' = 'FLS',
                    'tiller.count' = 'TC', 'tiller.no' = 'TN');
  
  for(cname in names(cname_list))
  {
    copydata = change_colname(copydata, cname, cname_list[[cname]]); 
  }
  
  return(copydata);
}

deleteSample = function(data, sample_list)
{
  copydata = data;
  i = c();
  for(samplename in sample_list)
  {
    i = append(i, which(copydata$barcode == samplename));
  }
  
  copydata = copydata[-i, ];
  return(copydata);
}


countNA = function()
{
  p = length(dta$GS39[is.na(dta$GS39)]) + 
  length(dta$GS39[is.na(dta$GS55)]) +
  length(dta$GS39[is.na(dta$GS65)]) +
         length(dta$GS39[is.na(dta$FLS)]) +
                length(dta$GS39[is.na(dta$FFLL)]);
  
  i = length(dta$GS39) * 5;
  print((p * 100) / i);
  
}

metafilename = 'w8_metadata.csv'; # Read metadata
metatable = read_data(metafilename, TRUE);
groundtruthfilename = 'w8_groundtruth.csv';
groundtruthtable = read_data(groundtruthfilename, TRUE);

# preproccesdata
groundtruthtable = changeCols(groundtruthtable);
groundtruthtable = preproccesdata(groundtruthtable);
groundtruthtable = data.frame(groundtruthtable, 
                              SM = round(apply(groundtruthtable[,c('disease_score1', 'disease_score2')], 1, mean)));
groundtruthtable = data.frame(groundtruthtable, 
                              NANTH = sapply(groundtruthtable$GS55, function(x) if(is.na(x)) { 0 } else { 1}  ));

# Merge metadata and groundtruth
data_all = merge(metatable, groundtruthtable, by.x=cbind('barcode', 'genotype'), 
                 by.y=cbind('barcode','genotype'));

data_all = convert2Numeric(data_all, colnames(data_all)[12:ncol(data_all)]);

dta = data_all;
dta = deleteSample(dta, c('W8-095111', 'W8-002111'));
dta = replaceNumber(dta, colnames(dta)[12:c(ncol(dta)-3)]);
#break();
dta_imp = replaceImputedValue(dta, c('GS39', 'GS55', 'GS65', 'FLS', 'FFLL'));
#dta_imp = replaceImputedValue(dta, colnames(dta)[12:22]);
densityplot(dta_imp$impdata)

png('imputeddata.png', width = 1600, height = 1600,res=200);
stripplot(dta_imp$impdata, pch = 20, cex = 1.2)
dev.off();
dta_avg = averageValues(dta_imp$copydata);


dtaM = createRatios(dta_avg);
write.table(dtaM[,c(1,19)], file='wheat_phenoS.csv', sep=',', row.names = F);

# Plot correlations
tiff('corplot.tiff',  width = 1980, height = 1580,res=200);
pairs.panels(dtaM[,3:21]) 
dev.off();

pv = calculatePvalue(dtaM[,4:21]);
writeHeatMatp(pv$pvaluep2Mean, pv$pvaluep2Mean, 'Mixed pop (r)', '#EEE5DE',colnames(dtaM[,4:21]), 1);
writeHeatMatp(pv$ccpop2Mean, pv$ccpop2Mean, 'Mixed pop (r)', '#EEE5DE',colnames(dtaM[,4:21]), 1);

# Plot scatter plots
genotypes = c(as.character(dtaM$genotype[1:12]), rep('RIL', nrow(dtaM)-12));
FFLL= c(dtaM$FFLL[1:12]/100, dtaM$FFLL[13:nrow(dtaM)]/100);


tiff('d1vsFLS.tiff', width = 1080, height = 1080,res=200);
p = plotDiagrams(dtaM$d1, dtaM$FLS, dtaM, FFLL, genotypes, 'd1', 'FLS', 'FFLL')
print(p);
dev.off();

FLS= c(dtaM$FLS[1:12]/100, dtaM$FLS[13:nrow(dtaM)]/100);
tiff('d2vsd3.tiff',  width = 1080, height = 1080,res=200);
p = plotDiagrams(dtaM$d2, dtaM$d3, dtaM, FLS, genotypes, 'd2', 'd3', 'FLS')
print(p);
dev.off();

genotypes = c(as.character(dtaM$genotype[1:12]), rep('RIL', nrow(dtaM)-12));
SM= c(dtaM$SM[1:12]/10, dtaM$SM[13:nrow(dtaM)]/10);
tiff('d3vsFLS.tiff',  width = 1080, height = 1080,res=200);
p = plotDiagrams(dtaM$d3, dtaM$FLS, dtaM, SM, genotypes, 'd3', 'FLS', 'SM');
print(p);
dev.off();


# fit model


fitmodel(dtaM)
  
# PLOT PCA
plotPCA(dtaM);

# Plot histograms
png('histrait.tiff', width = 1600, height = 1600,res=200);
plotHist(dtaM, colnames(dtaM)[3:21], c(4,5));
dev.off();

pdf('histrait.pdf', width = 600, height = 600);
plotHist(dtaM, colnames(dtaM)[3:21], c(4,5));
dev.off();

png('selectedtraits.tiff', width = 1600, height = 1600,res=200);
plotHist(dtaM, colnames(dtaM)[c(5,7,9, 21)], c(2,2));
dev.off();


#clusplot(dtaM[,3:21], fit$cluster, color=TRUE, shade=TRUE,          labels=2, lines=0)
barcol = c('blue', 'blue','blue', 'red', 'blue','blue', 'blue','blue', 'blue',
           'red','blue', 'blue');

png('barplotParentsFLS.png', width = 1000, height = 600,res=200);
barplot(dtaM$FLS[1:12], names.arg = dtaM$genotype[1:12], las=2, 
        cex.names = 0.7, col=barcol, ylab = 'FLS DAS', 
        ylim = c(0,300))
dev.off();
dev.off()
