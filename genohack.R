alleleLine <- function(a, n)
{
  v <- a / sum(a);
  return(sprintf("allele %s %s", n, paste(as.character(v), collapse = " ")));
}


alleleBlock <- function(s, markerName, cM)
{
  b= c();
  alleles <- unique(s);
  if(length(alleles) == 1) 
  {
    b <- sprintf("# marker %s: no polymorphic", markerName);
    return(b);
  }
  else
  {
    b <- sprintf("marker %s %d %f", markerName, length(alleles), cM);
    for (m in alleles)
    {
      an = as.numeric(s == m);
      b <- c(b, alleleLine(an, m));
      
    }
  }
  return(b);
}


g2a <- function(gfname, aBasename)
{
  #gAll <- read.table(gfname, header = TRUE, stringsAsFactors = FALSE, sep=',');
  gAll = gfname;
  chrNameList <- unique(gAll$pos);
  founderNames <- colnames(gAll)[12:ncol(gAll)];
  aList <- list();
  for (chrName in chrNameList)
  {
    g <- gAll[gAll$pos == chrName, ];
    a <- sprintf("markers %d strains %d", nrow(g), length(founderNames));
    a <- c(a, sprintf("strain_names %s", paste(founderNames, collapse = " ")));
    for (i in 1:nrow(g))
    {
      markerName <- g[i, "X"];
      s <- as.numeric(g[i, founderNames]);
      cM = as.numeric(g[i, 'cM'])
      names(s) <- founderNames;
      a <- c(a, alleleBlock(s, markerName, cM));
    }
    afname <- sprintf("chr%s%s", chrName, aBasename);
    write(a, afname);
    aList[[chrName]] <- a;
  }
  return(invisible(aList));
}


